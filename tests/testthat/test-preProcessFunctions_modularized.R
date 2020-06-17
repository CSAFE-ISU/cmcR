#TODO: - Write test to check that wavelength argument length agrees with
#filterType argument

tmp1 <- x3ptools::read_x3p(tmpfile1)

tmp2 <- tmp1$surface.matrix %>%
  magrittr::add(1) %>%
  cmcR::preProcess_ransac() %>%
  cmcR::preProcess_levelBF(useResiduals = FALSE)

testthat::test_that("preProcess_ functions work as expected", {
  testthat::expect_equal(tmp1$header.info$incrementY,2.5e-05)
  testthat::expect_equal(dim(tmp1$surface.matrix),c(145,144))

  testthat::expect_equal(as.numeric(quantile(tmp2,c(0,.5,1),na.rm = TRUE)),c(-2.143726e-05,7.115510e-07,2.005418e-05))

  testthat::expect_error(preProcess_ransac("not a matrix"))
  testthat::expect_error(preProcess_levelBF(tmp1$surface.matrix,
                                            useResiduals = logical(0)))

  testthat::expect_warning(testthat::expect_error(preProcess_cropWS(tmp1$surface.matrix,
                                             croppingThresh = numeric(0))))
  testthat::expect_error(preProcess_removeFPCircle(tmp1$surface.matrix,
                                                   aggregation_function = function(x){NULL}))
  testthat::expect_error(preProcess_gaussFilter(tmp1$surface.matrix,
                                                res = numeric(0)))
})
