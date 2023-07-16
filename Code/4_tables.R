source("0_packages.R")

hrs <- import("hrs_full_analytic.rds")

#get case/observation numbers
cases <- n_distinct(hrs$hhidpn)
obs <- nrow(hrs)


