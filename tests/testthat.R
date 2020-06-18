Sys.setenv("R_TESTS" = "")

library(testthat)
library(cmcR)

test_check("cmcR")
