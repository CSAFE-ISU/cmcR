
tmp1 <- x3ptools::read_x3p(tmpfile1)
tmp2 <- x3ptools::read_x3p(tmpfile2)

cmcs1 <- cmcR::congruentMatchingCells(tmp1,
                                      tmp2,
                                      thetas = c(-24,-21,-18,18,21,24),
                                      cellNumHoriz = 6,
                                      regionToCellProp = 4,
                                      centerCell = "individualCell",
                                      scaleCell = "individualCell",
                                      ccf_thresh = .4,
                                      dx_thresh = 20,
                                      theta_thresh = 6)

comp <- cellCCF_bothDirections(tmp1,
                               tmp2,
                               thetas = c(-24,-21,-18,18,21,24),
                               cellNumHoriz = 6,
                               regionToCellProp = 4,
                               centerCell = "individualCell",
                               scaleCell = "individualCell")

cmcs2 <- cmcFilter_improved(cellCCF_bothDirections_output = comp,
                            ccf_thresh = .4,
                            dx_thresh = 20,
                            theta_thresh = 6)

test_that("congruentMatchingCells function works as expected", {
  testthat::expect_identical(cmcs1,cmcs2)
})
