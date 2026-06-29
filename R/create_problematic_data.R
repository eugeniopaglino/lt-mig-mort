setwd('C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt_utilities')
library(here)
library(tinyplot)
library(data.table)

rm(list=ls())

i_am('create_problematic_data.R')

# 1. Create the table and base mortality
dt <- data.table(age = c(0, 1, seq(5, 90, by = 5)))
dt[, mx := 0.0002 * exp(0.08 * age)]

# 2. Add some problematic issues
dt[age == 0,  mx := 0.001]                 # Under-reported infant mortality
dt[age %between% c(20, 30), mx := mx * 3]    # Young adult crisis hump
dt[age == 70, mx := 0.1]                    # Age heaping spike
dt[age >= 85, mx := mx * 0.5]                # Implausible mortality deceleration

# 3. Save to CSV 
fwrite(dt, "problematic_abridged_mx.csv")

# Quick visual check
tinyplot(log(mx) ~ age, data = dt, type = "b", 
         main = "Problematic Abridged Mortality Profile",
         grid = TRUE, frame = FALSE)
