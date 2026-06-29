# Loading necessary packages
setwd("C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt_anneliese")
library(here)
library(mgcv)
library(tidyverse)
library(data.table)

# Restart R if necessary
rm(list=ls())

i_am('R/create_standard.R')

in_dir <- here('data')
out_dir <- here('outputs','data')

# Ensure the output directories exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here('outputs', 'figures'), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD AND PREP HMD DATA (Pooled 2000-2019)
# ==============================================================================

HMD_data <- fread(here('data', 'Mx_1x1.txt'))

# Filter for the requested pooled period
target_years <- 2000:2019
HMD_data <- HMD_data[Year %in% target_years]

HMD_data <- melt(
  HMD_data,
  id.vars = c('Year', 'Age'),
  measure.vars = c('Male', 'Female','Total'), # Pooled into sex-specific standards
  variable.name = 'sex',
  value.name = 'Mx'
)

setnames(HMD_data, c('Year', 'Age'), c('year', 'x'))

# Clean age formats
HMD_data[x == '110+', x := "110"]
HMD_data[, x := as.double(x)]

# Create continuous age mid-points
HMD_data[x < 110 & x != 0, x := x + 0.5]
HMD_data[x == 0, x := 0.2]
HMD_data[x == 110, x := x + 1.5]

# Convert to raw log mortality rates
HMD_data[, nMx_raw := log(Mx)]

# ==============================================================================
# 2. GAM SMOOTHING (Pooled by Sex)
# ==============================================================================

# Grouping by sex only for pooled estimation
unique_groups <- unique(HMD_data[, .(sex)])
gam_fits <- vector(mode = "list", length = nrow(unique_groups))

cat("Fitting GAM models (Pooled by Sex for 2000-2019)...\n")
for (i in 1:nrow(unique_groups)) {
  grp_data <- HMD_data[sex == unique_groups$sex[i]]
  
  # Fit model on the pooled data
  gam_fits[[i]] <- gam(
    nMx_raw ~ s(x, k = 30, bs = 'tp'),
    data = grp_data,
    method = 'REML'
  )
  cat("Finished fitting:", unique_groups$sex[i], "\n")
}

# ==============================================================================
# 3. PREDICT SMOOTHED DATA
# ==============================================================================

HMD_new_data <- expand.grid(
  sex = c('Male', 'Female','Total'),
  x = c(0.2, 1.5:111.5)
)
setDT(HMD_new_data)

cat("Predicting Smoothed HMD Standard rates...\n")
for (i in 1:nrow(unique_groups)) {
  target_sex <- unique_groups$sex[i]
  HMD_new_data[sex == target_sex, nMx_std := predict(gam_fits[[i]], newdata = HMD_new_data[sex == target_sex])]
}

# ==============================================================================
# 4. KANNISTO EXTRAPOLATION & BLENDING
# ==============================================================================

HMD_new_data[, Mx_smooth := exp(nMx_std)]
HMD_new_data[, Mx_kannisto := as.numeric(NA)]

cat("Fitting Kannisto extrapolations (Ages 80-105)...\n")
for (target_sex in c('Male', 'Female','Total')) {
  # Fit data: strictly ages 80 to 105
  fit_dt <- HMD_new_data[sex == target_sex & x >= 80 & x <= 105]
  fit_dt[, Mx_cap := pmin(pmax(Mx_smooth, 1e-6), 0.999)]
  
  mod <- lm(qlogis(Mx_cap) ~ I(x - 80), data = fit_dt)
  
  pred_dt <- HMD_new_data[sex == target_sex & x >= 80]
  preds <- plogis(predict(mod, newdata = pred_dt))
  
  HMD_new_data[sex == target_sex & x >= 80, Mx_kannisto := preds]
}

# Blend: Use GAM smoothed Mx for < 100, Kannisto for >= 105.
HMD_new_data[, weight_kannisto := fcase(
  x < 100, 0,
  x >= 100 & x < 105, (x - 99.5) / 6,
  default = 1
)]

HMD_new_data[, Mx_final := fcase(
  weight_kannisto == 0, Mx_smooth,
  weight_kannisto == 1, Mx_kannisto,
  default = (1 - weight_kannisto) * Mx_smooth + weight_kannisto * Mx_kannisto
)]

HMD_new_data[, log_mx_std := log(Mx_final)]
HMD_new_data[, c("Mx_smooth", "Mx_kannisto", "weight_kannisto") := NULL]

# ==============================================================================
# 5. DIAGNOSTIC PLOTS
# ==============================================================================

# Final Standard Check
pdf(here('outputs', 'figures', 'standard_diagnostic.pdf'), width = 8, height = 6)

p_final <- ggplot(HMD_new_data, aes(x = x, y = nMx_std, color = sex)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Final Smoothed & Blended HMD Standard Mortality Rates",
    subtitle = "Pooled: 2000-2019",
    x = "Age", y = "Log(nMx Standard)", color = "Sex"
  ) +
  theme(legend.position = "bottom")

print(p_final)
dev.off()

# ==============================================================================
# 6. EXPORT STANDARD
# ==============================================================================

HMD_new_data_export <- HMD_new_data[, .(sex, x, log_mx_std)]
setorder(HMD_new_data_export,sex,x)

fwrite(HMD_new_data_export, here(out_dir, 'single_standard.csv'))