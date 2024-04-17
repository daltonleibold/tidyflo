
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
