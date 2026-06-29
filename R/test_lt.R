setwd('C:/Users/epaglino/OneDrive - University of Helsinki/projects/lt_utilities')
library(here)
library(data.table)

i_am('test_lt.R')

source(here('lt_functions.R'))

# Test for abridged table
lt_abridged_data <- fread(here('lt_abridged_data.txt'))
my_lt_abridged <- lt_abridged(nMx=lt_abridged_data$mx)

# We are close, the differences are because the HMD uses different nax values
plot(
  lt_abridged_data$ex,my_lt_abridged$ex,
  xlab = "Original Life Expectancy",
  ylab = "Calculated Life Expectancy",
  main = "Comparing Life Expectancies"
  )
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

# Test for single-year table
lt_single_data <- fread(here('lt_single_data.txt'))
my_lt_single <- lt_single(Mx=lt_single_data$mx)

# These are identical because the ax values are exactly the same
plot(
  lt_single_data$ex,my_lt_single$ex,
  xlab = "Original Life Expectancy",
  ylab = "Calculated Life Expectancy",
  main = "Comparing Life Expectancies"
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

