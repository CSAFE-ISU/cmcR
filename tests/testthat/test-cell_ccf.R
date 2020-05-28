
tmp1 <- x3ptools::read_x3p(tmpfile1)
tmp2 <- x3ptools::read_x3p(tmpfile2)

comp <- cmcR::cellCCF_bothDirections(tmp1,
                                     tmp2,
                                     thetas = 0,
                                     cellNumHoriz = 6,
                                     regionToCellProp = 4,
                                     centerCell = "individualCell",
                                     scaleCell = "individualCell")

randInd <- rbinom(1,size = 1,prob = .5) + 1

testthat::test_that("cellCCF works as expected", {
  testthat::expect_equal(nrow(comp[[randInd]]$ccfResults$`0`),23)
  testthat::expect_false(any(is.na(comp[[randInd]]$ccfResults$`0`$fft.ccf)))
  testthat::expect_equal(comp[[randInd]]$params$theta,0)

  testthat::expect_error(cmcR::cellCCF(tmp1,"not a matrix"))
  testthat::expect_equal(cmcR::cellCCF(tmp1,
                                       tmp2,
                                       thetas = numeric(0)) %>%
                           .$ccfResults %>%
                           length(),
                         0)
})
