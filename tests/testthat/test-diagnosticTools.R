
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

cellRegPairs <- cmcR::getCellRegionPairs(tmp1,
                                         tmp2,
                                         comp$comparison_1to2$ccfResults %>%
                                           cmcR::topResultsPerCell(),
                                         cellCCF_params = comp$comparison_1to2$params)

p1 <- cmcR::ccfMapPlot(cellRegPairs[[1]][[1]][[1]],
                       cellRegPairs[[1]][[1]][[2]],
                       returnGrob = TRUE)

p2 <- cmcR::cmcPlot(x3p1 = tmp1,x3p2 = tmp2,
                    cellCCF_bothDirections_output = comp,
                    cmcFilter_improved_output = cmcs,
                    type = "list",
                    x3pNames = c("x3p1","x3p2"))

p3 <- cmcR::cmcPerThetaBarPlot(comp)

testthat::test_that("diagnosticTools functions work as expected", {
  testthat::expect_equal(cellRegPairs %>%
                           purrr::flatten() %>%
                           length(),
                         comp$comparison_1to2$ccfResults %>%
                           cmcR::topResultsPerCell() %>%
                           nrow())

  testthat::expect_s3_class(p1,class = c("gtable","gTree","grob","gDesc"))

  testthat::expect_true(all(unlist(purrr::map(purrr::flatten(p2), ~ class(.) == c("gg","ggplot")))))
  testthat::expect_equal(names(p2),c("initialCMC","highCMC"))

  testthat::expect_s3_class(p3,class = c("gg","ggplot"))
})
