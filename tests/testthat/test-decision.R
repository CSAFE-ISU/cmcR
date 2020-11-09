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

#make correlation unrealistically large so that comparison should yield 0 CMCs
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

# make correlation unrealistically large for the original method, but not for
# the High CMC -- the decision_combineDirections function with compareThetas =
# FALSE should yield same CMCs as cellTibble above, but 0 original CMCs in the
# reference_v_target direction
cellTibble_failedOriginal <- cmcR::comparison_allTogether(x3p1,x3p2,
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
                                                    corrThresh = .5,
                                                    thetaThresh = 6,
                                                    tau = 1))

#combining the CMCs under the various correlation thresholds used should change
#depending on the value of compareThetas
combinedCMCs <- cmcR::decision_combineDirections(cellTibble,cellTibble_rev)

combinedCMCs_failed <- cmcR::decision_combineDirections(cellTibble_failed,
                                                        cellTibble_rev,
                                                        compareThetas = TRUE)

combinedCMCs_failedOriginal <- cmcR::decision_combineDirections(cellTibble_failedOriginal,
                                                                cellTibble_rev,
                                                                compareThetas = FALSE)


testthat::test_that("decision_ functions work as expected", {
  #cmc counts should be equal for this limited example considering only 1 theta
  #value
  testthat::expect_equal(originalMethod_cmcCounts,high_cmcCounts)

  testthat::expect_true(all(nonCMC_butPassing$highCMCClassif == "non-CMC (passed)"))

  testthat::expect_true(all(cellTibble_failed$originalMethodClassif == "non-CMC"))
  testthat::expect_true(all(cellTibble_failed$highCMCClassif == "non-CMC (failed)"))

  #The reference_v_target direction original CMC count and High CMC count should
  #be 0 if compareThetas = TRUE and the correlation is excessively high
  testthat::expect_equal(nrow(combinedCMCs_failed$originalMethodCMCs[[1]]),0)
  testthat::expect_equal(nrow(combinedCMCs_failed$highCMCs),0)

  #setting compareThetas = FALSE means the High CMCs should be affected even if
  #original CMCs in reference_v_target direction are 0
  testthat::expect_equal(nrow(combinedCMCs$highCMCs),nrow(combinedCMCs_failedOriginal$highCMCs))

  #The target_v_reference direction original CMCs should be unaffected by
  #raising the reference_v_target direction correlation
  testthat::expect_true(identical(combinedCMCs$originalMethodCMCs[[2]],
                                  combinedCMCs_failed$originalMethodCMCs[[2]]))
})
