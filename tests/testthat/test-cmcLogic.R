
tmp1 <- x3ptools::read_x3p(tmpfile1)
tmp2 <- x3ptools::read_x3p(tmpfile2)

comp <- cmcR::cellCCF_bothDirections(tmp1,
                               tmp2,
                               thetas = c(-24,-21,-18,18,21,24),
                               cellNumHoriz = 6,
                               regionToCellProp = 4,
                               centerCell = "individualCell",
                               scaleCell = "individualCell")

cmcs <- cmcR::cmcFilter_improved(cellCCF_bothDirections_output = comp,
                           ccf_thresh = .4,
                           dx_thresh = 20,
                           theta_thresh = 6)

testthat::test_that("cmcFilter_ functions work as expected", {
  testthat::expect_equal(nrow(cmcs$initialCMCs$comparison_1to2),20)
  testthat::expect_equal(nrow(cmcs$initialCMCs$comparison_2to1),20)
  testthat::expect_equal(nrow(cmcs$highCMCs),24)

  testthat::expect_error(cmcFilter_improved(cellCCF_bothDirections_output = comp,ccf_thresh = numeric(0)))
  testthat::expect_error(cmcFilter_improved(cellCCF_bothDirections_output = comp,dx_thresh = numeric(0)))
  testthat::expect_error(cmcFilter_improved(cellCCF_bothDirections_output = comp,theta_thresh = numeric(0)))
  testthat::expect_error(cmcFilter_improved(cellCCF_bothDirections_output = comp,consensus_function = function(x){NULL}))
})
