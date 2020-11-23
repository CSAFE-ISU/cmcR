`%>%` <- dplyr::`%>%`

x3p1 <- x3ptools::read_x3p(tmpfile1)
x3p2 <- x3ptools::read_x3p(tmpfile2)

x3p1 <- x3p1 %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_crop(region = "interior",
                        radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                               tau = .5,
                               method = "fn") %>%
  cmcR::preProcess_gaussFilter(wavelength = c(16,500),
                               filtertype = "bp") %>%
  x3ptools::sample_x3p()

x3p2 <- x3p2 %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_crop(region = "interior",
                        radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                               tau = .5,
                               method = "fn") %>%
  cmcR::preProcess_gaussFilter(wavelength = c(16,500),
                               filtertype = "bp") %>%
  x3ptools::sample_x3p()

cellTibble <- cmcR::comparison_allTogether(x3p1,x3p2,
                                           theta = -24,
                                           numCells = 64,
                                           maxMissingProp = .85) %>%
  dplyr::mutate(originalMethodClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                           x = x,
                                                           y = y,
                                                           theta = theta,
                                                           corr = pairwiseCompCor,
                                                           xThresh = 20,
                                                           corrThresh = .5,
                                                           thetaThresh = 3),
                highCMCClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                    x = x,
                                                    y = y,
                                                    theta = theta,
                                                    corr = pairwiseCompCor,
                                                    xThresh = 20,
                                                    corrThresh = .5,
                                                    thetaThresh = 3,
                                                    tau = 1))

cellTibble_rev <- cmcR::comparison_allTogether(x3p2,x3p1,
                                               theta = 24,
                                               numCells = 64,
                                               maxMissingProp = .85) %>%
  dplyr::mutate(originalMethodClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                           x = x,
                                                           y = y,
                                                           theta = theta,
                                                           corr = pairwiseCompCor,
                                                           xThresh = 20,
                                                           corrThresh = .5,
                                                           thetaThresh = 3),
                highCMCClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                    x = x,
                                                    y = y,
                                                    theta = theta,
                                                    corr = pairwiseCompCor,
                                                    xThresh = 20,
                                                    corrThresh = .5,
                                                    thetaThresh = 3,
                                                    tau = 1))


x3pPlt <- cmcR::x3pListPlot(list("name1" = x3p1,
                                 "name2" = x3p2),
                            type = "list")

cmcPlt <- cmcR::cmcPlot(x3p1,
                        x3p2,
                        cellTibble,
                        cellTibble_rev,
                        x3pNames = c("name1","name2"),
                        corColName = "pairwiseCompCor")

cmcPlt_list <- cmcR::cmcPlot(x3p1,
                             x3p2,
                             cellTibble,
                             cellTibble_rev,
                             x3pNames = c("name1","name2"),
                             corColName = "pairwiseCompCor",type = "list")

testthat::test_that("diagnosticTools functions work as expected", {
  testthat::expect_named(x3pPlt,expected = c("name1","name2"))

  testthat::expect_true(all(unlist(purrr::map(x3pPlt, ~ class(.) == c("gg","ggplot")))))


  testthat::expect_named(cmcPlt,
                         expected = c("originalMethodCMCs_reference_v_target",
                                      "originalMethodCMCs_target_v_reference",
                                      "highCMC_reference_v_target",
                                      "highCMC_target_v_reference"))

  testthat::expect_true(all(unlist(purrr::map(cmcPlt, ~ class(.) == c("gg","ggplot")))))

  #Returning each plot individually rather than faceted:
  testthat::expect_named(cmcPlt_list,
                         expected = c("originalMethodCMCs_reference_v_target",
                                      "originalMethodCMCs_target_v_reference",
                                      "highCMC_reference_v_target",
                                      "highCMC_target_v_reference"))

  #individual plots should be named appropriately
  testthat::expect_true(all(purrr::map2_lgl(cmcPlt_list,
                                            list(c("name1","name2"),
                                                 c("name2","name1"),
                                                 c("name1","name2"),
                                                 c("name2","name1")),
                                            ~ assertthat::are_equal(names(.x),.y))))

  #add more "expect failure" tests?
})

