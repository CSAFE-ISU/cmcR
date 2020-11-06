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
  x3ptools::x3p_sample()

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
  x3ptools::x3p_sample()

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
                                                           thetaThresh = 6),
                highCMCClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                    x = x,
                                                    y = y,
                                                    theta = theta,
                                                    corr = pairwiseCompCor,
                                                    xThresh = 20,
                                                    corrThresh = .5,
                                                    thetaThresh = 6,
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
                                                           thetaThresh = 6),
                highCMCClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                    x = x,
                                                    y = y,
                                                    theta = theta,
                                                    corr = pairwiseCompCor,
                                                    xThresh = 20,
                                                    corrThresh = .5,
                                                    thetaThresh = 6,
                                                    tau = 1))

originalMethod_cmcCounts <- dplyr::bind_rows(cellTibble %>%
                                               dplyr::mutate(direction = "direction_1to2"),
                                             cellTibble_rev %>%
                                               dplyr::mutate(direction = "direction_2to1")) %>%
  dplyr::filter(originalMethodClassif == "CMC") %>%
  dplyr::group_by(direction) %>%
  dplyr::tally() %>%
  dplyr::pull(n)

high_cmcCounts <- dplyr::bind_rows(cellTibble %>%
                                     dplyr::mutate(direction = "direction_1to2"),
                                   cellTibble_rev %>%
                                     dplyr::mutate(direction = "direction_2to1")) %>%
  dplyr::filter(highCMCClassif == "CMC") %>%
  dplyr::group_by(direction) %>%
  dplyr::tally() %>%
  dplyr::pull(n)

#some cells will likely not be considered congruent even if the overall
#cartridge case pair passes the High CMC criterion.
nonCMC_butPassing <- dplyr::bind_rows(cellTibble %>%
                                        dplyr::mutate(direction = "direction_1to2"),
                                      cellTibble_rev %>%
                                        dplyr::mutate(direction = "direction_2to1")) %>%
  dplyr::filter(highCMCClassif != "CMC")

cellTibble_failed <- cmcR::comparison_allTogether(x3p1,x3p2,
                                                  theta = -24,
                                                  numCells = 64,
                                                  maxMissingProp = .85) %>%
  dplyr::mutate(originalMethodClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                           x = x,
                                                           y = y,
                                                           theta = theta,
                                                           corr = pairwiseCompCor,
                                                           xThresh = 20,
                                                           corrThresh = 1,
                                                           thetaThresh = 6),
                highCMCClassif = cmcR::decision_CMC(cellIndex = cellIndex,
                                                    x = x,
                                                    y = y,
                                                    theta = theta,
                                                    corr = pairwiseCompCor,
                                                    xThresh = 20,
                                                    corrThresh = 1,
                                                    thetaThresh = 6,
                                                    tau = 1))

testthat::test_that("decision_ functions work as expected", {
  #cmc counts should be equal for this limited example considering only 1 theta
  #value
  testthat::expect_equal(originalMethod_cmcCounts,high_cmcCounts)

  testthat::expect_true(all(nonCMC_butPassing$highCMCClassif == "non-CMC (passed)"))

  testthat::expect_true(all(cellTibble_failed$originalMethodClassif == "non-CMC"))
  testthat::expect_true(all(cellTibble_failed$highCMCClassif == "non-CMC (failed)"))

})
