test_data_raw <- readxl::read_excel("data_raw/test_data.xlsx")
d <- as.numeric(test_data_raw$`EoS status` == 0)
rho <- as.numeric(test_data_raw$`EoS status` != 2)
set.seed(100)
# epsilon <- 0.001 * (runif(nrow(test_data_raw)) - 0.5)
# tau <- test_data_raw$`Time to Death` + epsilon
tau <- test_data_raw$`Time to Death`
p <- test_data_raw$`Time to DP`
d_p <- as.numeric(test_data_raw$`DP Status` == 0)
usethis::use_data(tau, d, rho, p, d_p, overwrite = TRUE)


