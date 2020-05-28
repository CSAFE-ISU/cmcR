tmp1 <- x3ptools::read_x3p(tmpfile1)

test_that("preProcess_ functions work as expected", {
  testthat::expect_equal(tmp1$header.info$incrementY,2.5e-05)
  testthat::expect_equal(dim(tmp1$surface.matrix),c(145,144))

  testthat::expect_error(preProcess_ransac("not a matrix"))
  testthat::expect_error(preProcess_levelBF(tmp1$surface.matrix,
                                            useResiduals = logical(0)))
  testthat::expect_warning(
    testthat::expect_error(preProcess_cropWS(tmp1$surface.matrix,
                                           croppingThresh = numeric(0))))
  testthat::expect_error(preProcess_removeFPCircle(tmp1$surface.matrix,
                                                   aggregation_function = function(x){NULL}))
  testthat::expect_error(preProcess_gaussFilter(tmp1$surface.matrix,
                                                res = numeric(0)))
})
