
# Setup as package
install.packages("devtools"); library(devtools)

# Get it setup first by creating a package
#devtools::create("tidyflo")

# Setup license
#usethis::use_gpl3_license()

# Create the documentation
devtools::document()

# Now lets make accesible the functions
devtools::load_all()

# Check loaing
devtools::check()

# Create test folder
usethis::use_testthat()


# Example data that are used in tests
fcs_dat_1 <- data.frame(
  x = c(1, 2, 3, 4, 5),
  y = c(2, 4, 6, 8, 10)
)
usethis::use_data(fcs_dat_1, overwrite = TRUE)

devtools::test()