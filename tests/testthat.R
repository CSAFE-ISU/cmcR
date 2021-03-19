library(testthat)
library(cmcR)

test_check("cmcR",reporter = SummaryReporter)

# test_package("cmcR",reporter = SummaryReporter,filter="preProcess", invert=TRUE)
