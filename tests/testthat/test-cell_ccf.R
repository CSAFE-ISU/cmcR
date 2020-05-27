suppressWarnings({
  tmp1 <- x3ptools::read_x3p(tmpfile1)
  tmp2 <- x3ptools::read_x3p(tmpfile2)
})

comp <- cellCCF_bothDirections(tmp1,
                               tmp2,
                               thetas = 0,
                               cellNumHoriz = 6,
                               regionToCellProp = 4,
                               centerCell = "individualCell",
                               scaleCell = "individualCell")

test_that("cellCCF works as expected", {
  testthat::expect_equal(nrow(comp[[rbinom(1,size = 1,prob = .5)]]$ccfResults$`0`),36)
  testthat::expect_false(any(is.na(comp[[rbinom(1,size = 1,prob = .5)]]$ccfResults$`0`$ccf)))
  testthat::expect_equal(comp[[rbinom(1,size = 1,prob = .5)]]$params$theta,0)

  testthat::expect_error(cellCCF(tmp1,"not a matrix"))
  testthat::expect_equal(cellCCF(tmp1,
                                 tmp2,
                                 thetas = numeric(0)) %>%
                           .$ccfResults %>%
                           length(),
                         0)
})
