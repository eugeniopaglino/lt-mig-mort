setwd('C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt_anneliese')
library(here)
library(ggplot2)
library(data.table)
library(gt)

rm(list=ls())

i_am('R/test_count_models.R')

# Source the new count-based models 
source(here('R','mortality_models_counts.R'))
source(here('R','lt_functions.R'))

safe_gtsave <- function(gt_object, filename, ...) {
  # Check if file exists
  if (file.exists(filename)) {
    # Try to delete it
    tryCatch({
      file.remove(filename)
    }, warning = function(w) {
      stop("Could not overwrite file. It might be open in Word or locked by another process.")
    })
  }
  # Proceed to save
  gtsave(gt_object, filename, ...)
}

# Custom plotting theme for manuscript
theme_manuscript <- function() {
  theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# Reddish color palette inspired by the reference plot
fill_colors <- c("95% CI" = "#e2b8b8")
line_colors <- c("Model Mean" = "darkred", "Empirical Data" = "black", "Single-Year" = "black", "Abridged" = "black")

# Helper functions: Calculate midpoints for plotting interval rates
add_midpoints_abridged <- function(dt) {
  dt[, x_mid := age + 2.5]
  dt[age == 0, x_mid := 0.2]
  dt[age == 1, x_mid := 2.5]
  # For the open interval, cap it reasonably for plotting
  dt[age == max(age), x_mid := max(age) + 4.5] 
  return(dt)
}

add_midpoints_single <- function(dt) {
  dt[, x_mid := age + 0.5]
  dt[age == 0, x_mid := 0.2]
  # For the open interval, cap it reasonably for plotting
  dt[age == max(age), x_mid := max(age) + 1.5] 
  return(dt)
}

# =========================================================================
# TEST 1: BASIC-TESTING WITH HMD DATA
# =========================================================================

# Abridged data for test
exposures_abridged_data <- fread(here('data','exposures_abridged_data.txt'))
exposures_abridged_data[,c('female','male'):=NULL]
setnames(exposures_abridged_data,c('total'),c('Ex'))

deaths_abridged_data <- fread(here('data','deaths_abridged_data.txt'))
deaths_abridged_data[,c('female','male'):=NULL]
setnames(deaths_abridged_data,c('total'),c('Dx'))

lt_abridged_data <- merge(exposures_abridged_data,deaths_abridged_data,by=c('year','age'))

# Cleaning
lt_abridged_data[, x_start := as.integer(stringr::str_extract(age,'\\d+'))]
lt_abridged_data[, x_mid := x_start + 2.5]
lt_abridged_data[x_start == 0, x_mid := 0.2]
lt_abridged_data[x_start == 1, x_mid := 2.5]
lt_abridged_data[x_start == 110, x_mid := 111.5]

lt_abridged_data <- lt_abridged_data[, .(x = x_mid, Dx, Ex)]
lt_abridged_data[, mx := Dx / Ex]
setorder(lt_abridged_data, x)

# Standard
single_standard <- fread(here('outputs','data','single_standard.csv'))
# Select one standard, for females, males, or total also see the create_standard.R 
# script for a general framework to create a suitable standard from a set of mortality
# rates. Here we select total.
single_standard <- single_standard[sex=='Total']
single_standard[,mx_std:=exp(log_mx_std)]

# 1. Fit Single-Step Model and Extract Results using the new wrapper
model_output <- fit_abridged_sim_lt(
  mx_data = lt_abridged_data, 
  mx_std_data = single_standard,
  dof = 20,
  n_sim = 500,
  max_abridged_age = 85
)

my_lt_singularized <- model_output$single_year_summary
my_lt_singularized <- add_midpoints_single(my_lt_singularized)

my_lt_abridged <- model_output$abridged_summary
my_lt_abridged <- add_midpoints_abridged(my_lt_abridged)

# Plot: Age-Specific Mortality Rates
p_mx_basic <- ggplot() +
  geom_ribbon(data = my_lt_singularized, aes(x = x_mid, ymin = log(Mx_l95), ymax = log(Mx_u95), fill = "95% CI")) +
  geom_line(data = my_lt_singularized, aes(x = x_mid, y = log(Mx_mean), color = "Model Mean"), linewidth = 1) +
  geom_point(data = lt_abridged_data, aes(x = x, y = log(mx), color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Log Mortality Rate (log mx)") +
  theme_manuscript()

print(p_mx_basic)
ggsave(here('outputs', 'figures', 'p_mx_basic.png'), plot = p_mx_basic, width = 8, height = 6, dpi = 300, bg = "white")


# Baseline tables for comparison
my_lt_abridged <- lt_abridged(nMx = lt_abridged_data$mx)
setDT(my_lt_abridged)

lt_single_data <- fread(here('data','lt_single_data.txt'))
my_lt_single <- lt_single(Mx = lt_single_data$mx)
setDT(my_lt_single)

# Plot: Survival Function
p_lx_basic <- ggplot() +
  geom_ribbon(data = my_lt_singularized, aes(x = age, ymin = lx_l95, ymax = lx_u95, fill = "95% CI")) +
  geom_line(data = my_lt_single, aes(x = x, y = lx, color = "Single-Year"), linetype = "dashed") +
  geom_line(data = my_lt_singularized, aes(x = age, y = lx_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_lt_abridged, aes(x = x, y = lx, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Survival (lx)") +
  theme_manuscript()

print(p_lx_basic)
ggsave(here('outputs', 'figures', 'p_lx_basic.png'), plot = p_lx_basic, width = 8, height = 6, dpi = 300, bg = "white")

# Plot: Life Expectancy
p_ex_basic <- ggplot() +
  geom_ribbon(data = my_lt_singularized, aes(x = age, ymin = ex_l95, ymax = ex_u95, fill = "95% CI")) +
  geom_line(data = my_lt_single, aes(x = x, y = ex, color = "Single-Year"), linetype = "dashed") +
  geom_line(data = my_lt_singularized, aes(x = age, y = ex_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_lt_abridged, aes(x = x, y = ex, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Life Expectancy (ex)") +
  theme_manuscript()

print(p_ex_basic)
ggsave(here('outputs', 'figures', 'p_ex_basic.png'), plot = p_ex_basic, width = 8, height = 6, dpi = 300, bg = "white")


# =========================================================================
# TEST 2: STRESS-TESTING WITH PROBLEMATIC DATA
# =========================================================================

dt_prob <- fread(here('outputs','data','problematic_abridged_mx.csv'))
dt_prob[age==10, mx:=0] # Injecting a structural zero

dt_prob[, x_mid := age + 2.5]
dt_prob[age == 0, x_mid := 0.2]
dt_prob[age == 1, x_mid := 2.5]
dt_prob[age == 90, x_mid := 95.5] 

fake_exposures <- copy(lt_abridged_data)[,x:=fifelse(x>=90,95.5,x)][,.(Ex=sum(Ex)),by=.(x)]
fake_exposures[,Ex:=Ex/10000] # Radically shrinking exposures to simulate sparse data
dt_prob <- merge(dt_prob, fake_exposures[,.(x_mid=x, Ex)], by=c('x_mid'), all.x=T)
dt_prob[, Dx := floor(mx*Ex)]

prob_abridged <- dt_prob[, .(x = x_mid, Dx, Ex)]
prob_abridged[, mx := Dx / Ex]

# Run the wrapper model on sparse/problematic data
prob_model_output <- fit_abridged_sim_lt(
  mx_data = prob_abridged,
  mx_std_data = single_standard,
  dof = 15,
  skip_open = TRUE,
  skip_ages = list(c(0, 5), c(80, 100)),
  n_sim = 500,
  max_abridged_age = 85
)

my_prob_lt_singularized <- prob_model_output$single_year_summary
my_prob_lt_singularized <- add_midpoints_single(my_prob_lt_singularized)

# Plot: Problematic Data Mortality Rates
p_mx_prob <- ggplot() +
  geom_ribbon(data = my_prob_lt_singularized, aes(x = x_mid, ymin = log(Mx_l95), ymax = log(Mx_u95), fill = "95% CI")) +
  geom_line(data = my_prob_lt_singularized, aes(x = x_mid, y = log(Mx_mean), color = "Model Mean"), linewidth = 1) +
  geom_point(data = prob_abridged[mx > 0], aes(x = x, y = log(mx), color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = c("Model Mean" = "darkred", "Abridged" = "black")) +
  labs(x = "Age", y = "Log Mortality Rate (log mx)") +
  theme_manuscript()

print(p_mx_prob)
ggsave(here('outputs', 'figures', 'p_mx_prob.png'), plot = p_mx_prob, width = 8, height = 6, dpi = 300, bg = "white")

# Baseline comparison table
my_prob_lt_abridged <- lt_abridged(nMx = prob_abridged$mx)
setDT(my_prob_lt_abridged)

# Plot: Problematic Data Survival Function
p_lx_prob <- ggplot() +
  geom_ribbon(data = my_prob_lt_singularized, aes(x = age, ymin = lx_l95, ymax = lx_u95, fill = "95% CI")) +
  geom_line(data = my_prob_lt_singularized, aes(x = age, y = lx_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_prob_lt_abridged, aes(x = x, y = lx, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = c("Model Mean" = "darkred", "Abridged" = "black")) +
  labs(x = "Age", y = "Survival (lx)") +
  theme_manuscript()

print(p_lx_prob)
ggsave(here('outputs', 'figures', 'p_lx_prob.png'), plot = p_lx_prob, width = 8, height = 6, dpi = 300, bg = "white")


# Plot: Problematic Data Life Expectancy
p_ex_prob <- ggplot() +
  geom_ribbon(data = my_prob_lt_singularized, aes(x = age, ymin = ex_l95, ymax = ex_u95, fill = "95% CI")) +
  geom_line(data = my_prob_lt_singularized, aes(x = age, y = ex_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_prob_lt_abridged, aes(x = x, y = ex, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = c("Model Mean" = "darkred", "Abridged" = "black")) +
  labs(x = "Age", y = "Life Expectancy (ex)") +
  theme_manuscript()

print(p_ex_prob)
ggsave(here('outputs', 'figures', 'p_ex_prob.png'), plot = p_ex_prob, width = 8, height = 6, dpi = 300, bg = "white")


# =========================================================================
# TEST 3: APPLICATION TO RUSSIA 2000 DATA
# =========================================================================

# 1. Load and clean abridged data for Russia
exposures_abridged_data_rus <- fread(here('data','exposures_abridged_data_RUS.txt'))
exposures_abridged_data_rus[, c('female','male') := NULL]
setnames(exposures_abridged_data_rus, c('total'), c('Ex'))

deaths_abridged_data_rus <- fread(here('data','deaths_abridged_data_RUS.txt'))
deaths_abridged_data_rus[, c('female','male') := NULL]
setnames(deaths_abridged_data_rus, c('total'), c('Dx'))

lt_abridged_data_rus <- merge(exposures_abridged_data_rus, deaths_abridged_data_rus, by=c('year','age'))

lt_abridged_data_rus[, x_start := as.integer(stringr::str_extract(age,'\\d+'))]
lt_abridged_data_rus[, x_mid := x_start + 2.5]
lt_abridged_data_rus[x_start == 0, x_mid := 0.2]
lt_abridged_data_rus[x_start == 1, x_mid := 2.5]
lt_abridged_data_rus[x_start == 110, x_mid := 111.5]

lt_abridged_data_rus <- lt_abridged_data_rus[, .(x = x_mid, Dx, Ex)]
lt_abridged_data_rus[, mx := Dx / Ex]
setorder(lt_abridged_data_rus, x)

# 2. Fit Single-Step Model
model_output_rus <- fit_abridged_sim_lt(
  mx_data = lt_abridged_data_rus, 
  mx_std_data = single_standard,
  n_sim = 500,
  max_abridged_age = 85
)

my_lt_singularized_rus <- model_output_rus$single_year_summary
my_lt_singularized_rus <- add_midpoints_single(my_lt_singularized_rus)

# Plot: Russia Mortality Rates
p_mx_rus <- ggplot() +
  geom_ribbon(data = my_lt_singularized_rus, aes(x = x_mid, ymin = log(Mx_l95), ymax = log(Mx_u95), fill = "95% CI")) +
  geom_line(data = my_lt_singularized_rus, aes(x = x_mid, y = log(Mx_mean), color = "Model Mean"), linewidth = 1) +
  geom_point(data = lt_abridged_data_rus, aes(x = x, y = log(mx), color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Log Mortality Rate (log mx)") +
  theme_manuscript()

print(p_mx_rus)
ggsave(here('outputs', 'figures', 'p_mx_rus.png'), plot = p_mx_rus, width = 8, height = 6, dpi = 300, bg = "white")

# Baseline tables for Russia (Abridged)
my_lt_abridged_rus <- lt_abridged(nMx = lt_abridged_data_rus$mx)
setDT(my_lt_abridged_rus)

# Baseline tables for Russia (True Single-Year from HMD)
lt_single_data_rus <- fread(here('data','lt_single_data_RUS.txt'))
my_lt_single_rus <- lt_single(Mx = lt_single_data_rus$mx)
setDT(my_lt_single_rus)

# Plot: Russia Survival Function
p_lx_rus <- ggplot() +
  geom_ribbon(data = my_lt_singularized_rus, aes(x = age, ymin = lx_l95, ymax = lx_u95, fill = "95% CI")) +
  geom_line(data = my_lt_single_rus, aes(x = x, y = lx, color = "Single-Year"), linetype = "dashed") +
  geom_line(data = my_lt_singularized_rus, aes(x = age, y = lx_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_lt_abridged_rus, aes(x = x, y = lx, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Survival (lx)") +
  theme_manuscript()

print(p_lx_rus)
ggsave(here('outputs', 'figures', 'p_lx_rus.png'), plot = p_lx_rus, width = 8, height = 6, dpi = 300, bg = "white")

# Plot: Russia Life Expectancy
p_ex_rus <- ggplot() +
  geom_ribbon(data = my_lt_singularized_rus, aes(x = age, ymin = ex_l95, ymax = ex_u95, fill = "95% CI")) +
  geom_line(data = my_lt_single_rus, aes(x = x, y = ex, color = "Single-Year"), linetype = "dashed") +
  geom_line(data = my_lt_singularized_rus, aes(x = age, y = ex_mean, color = "Model Mean"), linewidth = 1) +
  geom_point(data = my_lt_abridged_rus, aes(x = x, y = ex, color = "Abridged"), size = 1.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  labs(x = "Age", y = "Life Expectancy (ex)") +
  theme_manuscript()

print(p_ex_rus)
ggsave(here('outputs', 'figures', 'p_ex_rus.png'), plot = p_ex_rus, width = 8, height = 6, dpi = 300, bg = "white")

# =========================================================================
# 3. GENERATING PUBLICATION TABLES WITH gt
# =========================================================================

# Define exact starting ages to display in the manuscript tables
target_ages <- c(0, 10, 20, 40, 60, 80, 90)

# Create a standardized 'exact_age' key
my_lt_singularized[, exact_age := age]
my_lt_abridged[, exact_age := floor(x)]
my_lt_single[, exact_age := floor(x)]

my_prob_lt_singularized[, exact_age := age]
my_prob_lt_abridged[, exact_age := floor(x)]

my_lt_singularized_rus[, exact_age := age]
my_lt_abridged_rus[, exact_age := floor(x)]
my_lt_single_rus[, exact_age := floor(x)]

# -------------------------------------------------------------------------
# TABLE 1: Basic Test Summary
# -------------------------------------------------------------------------
basic_summary <- merge(
  my_lt_abridged[exact_age %in% target_ages, .(exact_age, lx_abridged = lx, ex_abridged = ex)],
  my_lt_single[exact_age %in% target_ages, .(exact_age, lx_single = lx, ex_single = ex)],
  by = "exact_age", all = TRUE
)

basic_summary <- merge(
  basic_summary,
  my_lt_singularized[exact_age %in% target_ages, .(
    exact_age,
    lx_sing_mean = lx_mean, lx_sing_l95 = lx_l95, lx_sing_u95 = lx_u95,
    ex_sing_mean = ex_mean, ex_sing_l95 = ex_l95, ex_sing_u95 = ex_u95
  )],
  by = "exact_age", all = TRUE
)

basic_summary[, lx_abridged_fmt := sprintf("%.1f", lx_abridged * 100)]
basic_summary[, lx_single_fmt := sprintf("%.1f", lx_single * 100)]
basic_summary[, lx_sing_fmt := sprintf("%.1f [%.1f, %.1f]", 
                                       lx_sing_mean * 100, 
                                       lx_sing_l95 * 100, 
                                       lx_sing_u95 * 100)]

basic_summary[, ex_abridged_fmt := sprintf("%.2f", ex_abridged)]
basic_summary[, ex_single_fmt := sprintf("%.2f", ex_single)]
basic_summary[, ex_sing_fmt := sprintf("%.2f [%.2f, %.2f]", 
                                       ex_sing_mean, ex_sing_l95, ex_sing_u95)]

setnames(basic_summary, "exact_age", "Age")

basic_gt <- basic_summary[, .(Age, lx_abridged_fmt, lx_single_fmt, lx_sing_fmt, ex_abridged_fmt, ex_single_fmt, ex_sing_fmt)] |>
  gt() |>
  tab_spanner(
    label = "Survivors (%)",
    columns = c(lx_abridged_fmt, lx_single_fmt, lx_sing_fmt)
  ) |>
  tab_spanner(
    label = "Life Expectancy (ex)",
    columns = c(ex_abridged_fmt, ex_single_fmt, ex_sing_fmt)
  ) |>
  cols_label(
    lx_abridged_fmt = "Abridged",
    lx_single_fmt = "Single-Year",
    lx_sing_fmt = "Modeled [95% CI]",
    ex_abridged_fmt = "Abridged",
    ex_single_fmt = "Single-Year",
    ex_sing_fmt = "Modeled [95% CI]"
  ) |>
  cols_align(align = "center", columns = -Age) |>
  sub_missing(columns = everything(), missing_text = "—") |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

safe_gtsave(basic_gt, here('outputs', 'tables', 'Table_1_Basic_Test.docx'))


# -------------------------------------------------------------------------
# TABLE 2: Stress Test Summary
# -------------------------------------------------------------------------
stress_summary <- merge(
  my_prob_lt_abridged[exact_age %in% target_ages, .(exact_age, lx_abridged = lx, ex_abridged = ex)],
  my_prob_lt_singularized[exact_age %in% target_ages, .(
    exact_age,
    lx_sing_mean = lx_mean, lx_sing_l95 = lx_l95, lx_sing_u95 = lx_u95,
    ex_sing_mean = ex_mean, ex_sing_l95 = ex_l95, ex_sing_u95 = ex_u95
  )],
  by = "exact_age", all = TRUE
)

stress_summary[, lx_abridged_fmt := sprintf("%.1f", lx_abridged * 100)]
stress_summary[, lx_sing_fmt := sprintf("%.1f [%.1f, %.1f]", 
                                        lx_sing_mean * 100, 
                                        lx_sing_l95 * 100, 
                                        lx_sing_u95 * 100)]

stress_summary[, ex_abridged_fmt := sprintf("%.2f", ex_abridged)]
stress_summary[, ex_sing_fmt := sprintf("%.2f [%.2f, %.2f]", 
                                        ex_sing_mean, ex_sing_l95, ex_sing_u95)]
setnames(stress_summary, "exact_age", "Age")

stress_gt <- stress_summary[, .(Age, lx_abridged_fmt, lx_sing_fmt, ex_abridged_fmt, ex_sing_fmt)] |>
  gt() |>
  tab_spanner(
    label = "Survivors (%)",
    columns = c(lx_abridged_fmt, lx_sing_fmt)
  ) |>
  tab_spanner(
    label = "Life Expectancy (ex)",
    columns = c(ex_abridged_fmt, ex_sing_fmt)
  ) |>
  cols_label(
    lx_abridged_fmt = "Abridged (Problematic)",
    lx_sing_fmt = "Modeled [95% CI]",
    ex_abridged_fmt = "Abridged (Problematic)",
    ex_sing_fmt = "Modeled [95% CI]"
  ) |>
  cols_align(align = "center", columns = -Age) |>
  sub_missing(columns = everything(), missing_text = "—") |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

safe_gtsave(stress_gt, here('outputs', 'tables', 'Table_2_Stress_Test.docx'))

# -------------------------------------------------------------------------
# TABLE 3: Russia 2000 Summary
# -------------------------------------------------------------------------
rus_summary <- merge(
  my_lt_abridged_rus[exact_age %in% target_ages, .(exact_age, lx_abridged = lx, ex_abridged = ex)],
  my_lt_single_rus[exact_age %in% target_ages, .(exact_age, lx_single = lx, ex_single = ex)],
  by = "exact_age", all = TRUE
)

rus_summary <- merge(
  rus_summary,
  my_lt_singularized_rus[exact_age %in% target_ages, .(
    exact_age,
    lx_sing_mean = lx_mean, lx_sing_l95 = lx_l95, lx_sing_u95 = lx_u95,
    ex_sing_mean = ex_mean, ex_sing_l95 = ex_l95, ex_sing_u95 = ex_u95
  )],
  by = "exact_age", all = TRUE
)

rus_summary[, lx_abridged_fmt := sprintf("%.1f", lx_abridged * 100)]
rus_summary[, lx_single_fmt := sprintf("%.1f", lx_single * 100)]
rus_summary[, lx_sing_fmt := sprintf("%.1f [%.1f, %.1f]", 
                                     lx_sing_mean * 100, lx_sing_l95 * 100, lx_sing_u95 * 100)]
rus_summary[, ex_abridged_fmt := sprintf("%.2f", ex_abridged)]
rus_summary[, ex_single_fmt := sprintf("%.2f", ex_single)]
rus_summary[, ex_sing_fmt := sprintf("%.2f [%.2f, %.2f]", 
                                     ex_sing_mean, ex_sing_l95, ex_sing_u95)]
setnames(rus_summary, "exact_age", "Age")

rus_gt <- rus_summary[, .(Age, lx_abridged_fmt, lx_single_fmt, lx_sing_fmt, ex_abridged_fmt, ex_single_fmt, ex_sing_fmt)] |>
  gt() |>
  tab_header(
    title = "Table 3. Application: Russia 2000 Data",
    subtitle = "Abridged vs. Single (True) vs. Singularized Model (with 95% Credible Intervals)"
  ) |>
  tab_spanner(
    label = "Survivors (%)",
    columns = c(lx_abridged_fmt, lx_single_fmt, lx_sing_fmt)
  ) |>
  tab_spanner(
    label = "Life Expectancy (ex)",
    columns = c(ex_abridged_fmt, ex_single_fmt, ex_sing_fmt)
  ) |>
  cols_label(
    lx_abridged_fmt = "Abridged",
    lx_single_fmt = "Single-Year",
    lx_sing_fmt = "Modeled [95% CI]",
    ex_abridged_fmt = "Abridged",
    ex_single_fmt = "Single-Year",
    ex_sing_fmt = "Modeled [95% CI]"
  ) |>
  cols_align(align = "center", columns = -Age) |>
  sub_missing(columns = everything(), missing_text = "—") |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

safe_gtsave(rus_gt, here('outputs', 'tables', 'Table_3_Russia_2000.docx'))