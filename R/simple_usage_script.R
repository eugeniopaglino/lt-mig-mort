setwd('C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt-mig-mort')
library(here)
library(data.table)

rm(list=ls())

# 1. Setup & Load Utilities
# Ensure these scripts are in your R folder
source(here('R', 'lt_functions.R'))
source(here('R', 'mortality_models_counts.R'))

# 2. Prepare your Data
# Ensure your data table has columns: x (age midpoint), Dx (deaths), Ex (exposures)
my_data <- fread(here('data', 'your_abridged_deaths_and_exposures.csv')) 
# Drop extra columns, change names, and create midpoints
setnames(my_data,c('Age','nDx','nEx'),c('x','Dx','Ex'))
my_data[,Year:=NULL]

my_data[, x_mid := x + 2.5]
my_data[x == 0, x_mid := 0.2]
my_data[x == 1, x_mid := 2.5]
# For the open-ended interval, you can use last_age + 1/nMx
my_data[x == max(x), x_mid := max(x) + (1/(Dx/Ex))]

# leave x as the age column, representing the representative age for
# each age interval
my_data[,x:=NULL]
setnames(my_data,'x_mid','x')

# Ensure your standard has columns: x (age), mx_std (mortality rate)
my_standard <- fread(here('outputs', 'data', 'single_standard.csv'))
# Select appropriate standard (by sex)
my_standard <- my_standard[sex == 'Total', .(x, mx_std = exp(log_mx_std))]

# 3. Fit the Model
# n_sim = 500 is usually enough for a quick diagnostic
model_results <- fit_abridged_sim_lt(
  mx_data = my_data, 
  mx_std_data = my_standard,
  skip_open = TRUE, # To skip open-ended interval (generally recommended)
  # skip_ages = list(c(0, 5), c(80, 100)), # Use to skip certain age intervals
  n_sim = 500,
  max_abridged_age = 85
)

# 4. Extract exclusively the Abridged Summary
# This contains the mean, l95, and u95 for each life table column
abridged_output <- model_results$abridged_summary

# 5. View or Export
print(abridged_output)

# Optional: Export to CSV
fwrite(abridged_output, here('outputs', 'data', 'abridged_summary_results.csv'))