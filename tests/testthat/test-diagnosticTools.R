`%>%` <- dplyr::`%>%`

x3p1 <- x3ptools::read_x3p(tmpfile1)
x3p2 <- x3ptools::read_x3p(tmpfile2)

if(!exists("skipPreprocess")){
  x3p1 <- x3p1 %>%
    cmcR::preProcess_crop(region = "exterior") %>%
    cmcR::preProcess_crop(region = "interior") %>%
    cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
    cmcR::preProcess_gaussFilter(wavelength = c(16,500),
                                 filtertype = "bp") %>%
    x3ptools::sample_x3p()

  x3p2 <- x3p2 %>%
    cmcR::preProcess_crop(region = "exterior") %>%
    cmcR::preProcess_crop(region = "interior") %>%
    cmcR::preProcess_removeTrend(statistic = "quantile",
                                 tau = .5,
                                 method = "fn") %>%
    cmcR::preProcess_gaussFilter(wavelength = c(16,500),
                                 filtertype = "bp") %>%
    x3ptools::sample_x3p()
}

cellTibble <- cmcR::comparison_allTogether(x3p1,x3p2,
                                           theta = -24,
                                           numCells = c(8,8),
                                           maxMissingProp = .85,
                                           returnX3Ps = TRUE) %>%
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
                                               numCells = c(8,8),
                                               maxMissingProp = .85,
                                               returnX3Ps = TRUE) %>%
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

cmcPlt <- cmcR::cmcPlot(reference = x3p1,
                        target = x3p2,
                        cmcClassifs = cellTibble %>%
                          dplyr::filter(highCMCClassif == "CMC"),
                        cmcCol = "highCMCClassif")

cmcPlt_list <- cmcR::cmcPlot(reference = x3p1,
                             target = x3p2,
                             cmcClassifs = cellTibble %>%
                               dplyr::filter(highCMCClassif == "CMC"),
                             cmcCol = "highCMCClassif",
                             type = "list")

testthat::test_that("diagnosticTools functions work as expected", {
  testthat::expect_named(x3pPlt,expected = c("name1","name2"))

  testthat::expect_true(all(unlist(purrr::map(x3pPlt, ~ class(.) == c("gg","ggplot")))))


  # testthat::expect_named(cmcPlt,
  #                        expected = c("originalMethodCMCs_reference_v_target",
  #                                     "originalMethodCMCs_target_v_reference",
  #                                     "highCMC_reference_v_target",
  #                                     "highCMC_target_v_reference"))

  testthat::expect_true(all(unlist(purrr::map(cmcPlt_list, ~ class(.) == c("gg","ggplot")))))

  #Returning each plot individually rather than faceted:
  testthat::expect_named(cmcPlt_list,
                         expected = c("reference","target","legend"))

  #individual plots should be named appropriately
  # testthat::expect_true(all(purrr::map2_lgl(cmcPlt_list,
  #                                           list(c("name1","name2"),
  #                                                c("name2","name1"),
  #                                                c("name1","name2"),
  #                                                c("name2","name1")),
  #                                           ~ assertthat::are_equal(names(.x),.y))))

  #add more "expect failure" tests?
})
