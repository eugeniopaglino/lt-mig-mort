library(mgcv)
library(data.table)
library(MASS) # Required for mvrnorm

#' Single-Step Relational Mortality Model with Posterior Simulations
#' 
#' @param mx_data data.table with columns x, Dx, Ex
#' @param mx_std_data data.table with columns x, mx_std
#' @param dof Spline basis dimension
#' @param basis Spline basis type
#' @param skip_open Logical. If TRUE, treats highest age Dx as NA.
#' @param skip_ages List of numeric vectors defining age bounds to exclude.
#' @param ages Numeric vector for prediction grid.
#' @param n_sim Integer. Number of posterior simulations to draw.
#' @return A list containing `estimates` (data.table) and `simulations` (data.table).
singularize_single_step <- function(
    mx_data, 
    mx_std_data,
    dof = 10,
    basis = 'tp',
    skip_open = FALSE,
    skip_ages = NULL,
    ages = c(0.2, 1.5:111.5),
    n_sim = 1000
) {
  
  # 1. Prepare Data
  data <- merge(mx_data, mx_std_data, by = c('x'), all.x = TRUE)
  data[, log_mx_std := log(mx_std)]
  
  if (skip_open) {
    data[x == max(x), Dx := NA]
  }
  
  if (!is.null(skip_ages)) {
    if (!is.list(skip_ages)) skip_ages <- list(skip_ages)
    for (age_range in skip_ages) {
      if (is.numeric(age_range) && length(age_range) == 2) {
        data[x %between% age_range, Dx := NA]
      }
    }
  }
  
  # 2. Fit the Single-Step Model
  # mgcv drops NA response rows safely during fitting but keeps the continuous grid intact.
  gam_fit <- gam(
    Dx ~ 1 + log_mx_std + s(x, k = dof, bs = basis) + offset(log(Ex)), 
    family = nb(),
    data = data,
    method='REML',
    select=T
  )
  
  out_list <- list(model_fit = gam_fit)
  
  # 3. Prepare Data for Prediction
  new_data <- data.table(x = ages)
  new_data <- merge(new_data, mx_std_data, by = 'x')
  new_data[, log_mx_std := log(mx_std)]
  
  # Satisfy mgcv's formula parser and mathematically nullify the offset
  # Setting Ex = 1 means log(Ex) = 0, so the prediction represents pure log(mx)
  new_data[, Ex := 1]
  
  # 4. Extract Point Estimates
  # type = "lpmatrix" returns the evaluated basis functions (the design matrix Xp).
  # Crucially, Xp does *not* include the offset. Therefore, Xp %*% beta yields
  # pure log(mx), completely ignoring log(Ex).
  Xp <- predict(gam_fit, newdata = new_data, type = "lpmatrix")
  
  # Point estimate of log(mx) and exponentiate to mx
  eta_hat <- Xp %*% coef(gam_fit)
  new_data[, mx_s := as.numeric(exp(eta_hat))]
  
  out_list$estimates <- new_data
  
  # 5. Generate Posterior Simulations
  if (n_sim > 0) {
    # Draw coefficients from the empirical Bayesian posterior
    # gam_fit$Vp is the posterior covariance matrix of the coefficients
    beta_sims <- MASS::mvrnorm(n_sim, mu = coef(gam_fit), Sigma = gam_fit$Vp)
    
    # Calculate simulated linear predictors: Xp (Ages x Coefs) %*% t(beta_sims) (Coefs x Sims)
    # The result is a matrix of size: (Number of Ages) x (n_sim)
    eta_sims <- Xp %*% t(beta_sims)
    
    # Convert to pure mortality rates
    rate_sims <- exp(eta_sims)
    
    # Format as a data.table and bind with the age column
    sim_dt <- as.data.table(rate_sims)
    setnames(sim_dt, paste0("sim_", 1:n_sim))
    out_list$simulations <- cbind(new_data[, .(x)], sim_dt)
  }
  
  return(out_list)
}

#' Fit Single-Step Model and Return Full Abridged Life Table Summaries
#' 
#' @param mx_data data.table with columns x, Dx, Ex
#' @param mx_std_data data.table with columns x, mx_std
#' @param dof Spline basis dimension
#' @param basis Spline basis type
#' @param skip_open Logical. If TRUE, treats highest age Dx as NA.
#' @param skip_ages List of numeric vectors defining age bounds to exclude.
#' @param ages Numeric vector for prediction grid.
#' @param n_sim Integer. Number of posterior simulations to draw.
#' @param max_abridged_age Integer. The starting age of your open-ended interval.
#' @return A list containing the raw `model_fit`, an `abridged_summary` data.table
#' and a `single_year_summary` data.table
fit_abridged_sim_lt <- function(
    mx_data, 
    mx_std_data,
    dof = 10,
    basis = 'tp',
    skip_open = FALSE,
    skip_ages = NULL,
    ages = c(0.2, 1.5:111.5),
    n_sim = 1000,
    max_abridged_age = 85
) {
  
  # 1. Run the base model to get single-year simulations
  model_out <- singularize_single_step(
    mx_data = mx_data, mx_std_data = mx_std_data, dof = dof, 
    basis = basis, skip_open = skip_open, skip_ages = skip_ages, 
    ages = ages, n_sim = n_sim
  )
  
  simulations_dt <- model_out$simulations
  sim_cols <- grep("sim_", names(simulations_dt), value = TRUE)
  
  N_ages <- nrow(simulations_dt)
  exact_ages <- 0:(N_ages - 1)
  
  # 2. Define setups for both summaries
  # Abridged
  abridged_starts <- c(0, 1, seq(5, max_abridged_age, by = 5))
  age_groups <- cut(exact_ages, breaks = c(abridged_starts, Inf), right = FALSE, labels = abridged_starts)
  
  metrics <- c("nMx", "nax", "nqx", "npx", "ndx", "lx", "nLx", "Tx", "ex")
  single_metrics <- c("Mx", "lx", "ex") # You can add others if needed
  
  # 3. Setup matrices
  n_sims <- length(sim_cols)
  
  # Matrices for Abridged
  abridged_mats <- setNames(lapply(metrics, function(m) matrix(NA, nrow = length(abridged_starts), ncol = n_sims)), metrics)
  # Matrices for Single-Year (Up to 110)
  single_mats <- setNames(lapply(single_metrics, function(m) matrix(NA, nrow = length(exact_ages), ncol = n_sims)), single_metrics)
  
  # 4. Loop through simulations
  for (i in seq_along(sim_cols)) {
    mx_vec <- simulations_dt[[sim_cols[i]]]
    
    # Build single year LT
    slt <- lt_single(Mx = mx_vec, x = exact_ages)
    setDT(slt)
    
    # -- A. Handle Abridged --
    slt[, age_group := age_groups]
    agg <- slt[, .(nLx = sum(Lx), ndx = sum(dx), lx = lx[1], Tx = Tx[1], ex = ex[1], n = .N), by = age_group]
    
    template <- data.table(age_group = factor(abridged_starts, levels = abridged_starts))
    agg <- merge(agg, template, by = "age_group", all.y = TRUE)
    setorder(agg, age_group)
    
    agg[, nMx := ndx / nLx]
    agg[, nqx := ndx / lx]
    agg[.N, nqx := 1] 
    agg[, npx := 1 - nqx]
    agg[, nax := fifelse(ndx > 0, (nLx - n * (lx - ndx)) / ndx, n / 2)]
    agg[.N, nax := ex]
    
    for(m in metrics) abridged_mats[[m]][, i] <- agg[[m]]
    
    # -- B. Handle Single-Year --
    for(m in single_metrics) single_mats[[m]][, i] <- slt[[m]]
  }
  
  # 5. Summarize functions
  summarize_mats <- function(mat_list, age_vec) {
    summary_list <- lapply(names(mat_list), function(m) {
      mat <- mat_list[[m]]
      dt_stats <- data.table(
        mean = apply(mat, 1, mean, na.rm = TRUE),
        l95  = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
        u95  = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
      )
      setnames(dt_stats, c(paste0(m, "_mean"), paste0(m, "_l95"), paste0(m, "_u95")))
      return(dt_stats)
    })
    cbind(data.table(age = age_vec), do.call(cbind, summary_list))
  }
  
  return(list(
    model_fit = model_out$model_fit,
    abridged_summary = summarize_mats(abridged_mats, abridged_starts),
    single_year_summary = summarize_mats(single_mats, exact_ages)
  ))
}