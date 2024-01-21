# List of required packages
required_packages <- c(
  "maps", "ggplot2", "ggrepel", "lubridate", "zoo", "tidyverse", "docstring", "geosphere",
  "corrplot", "lmodel2", "tseries", "forecast", "magrittr", "testthat"
)

# Check if each package is already installed, and install if not
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
}
