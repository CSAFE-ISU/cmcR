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

#Perform entire comparison procedure explicitly
cellTibble <- x3p1 %>%
  cmcR::comparison_cellDivision(numCells = 64) %>%
  dplyr::mutate(regionHeightValues = cmcR::comparison_getTargetRegions(cellHeightValues = cellHeightValues,
                                                          target = x3p2,
                                                          theta = 0)) %>%
  dplyr::mutate(cellPropMissing = cmcR::comparison_calcPropMissing(cellHeightValues),
         regionPropMissing = cmcR::comparison_calcPropMissing(regionHeightValues)) %>%
  dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
  dplyr::mutate(cellHeightValues = cmcR::comparison_standardizeHeights(cellHeightValues),
         regionHeightValues = cmcR::comparison_standardizeHeights(regionHeightValues)) %>%
  dplyr::mutate(cellHeightValues_replaced = cmcR::comparison_replaceMissing(cellHeightValues),
         regionHeightValues_replaced = cmcR::comparison_replaceMissing(regionHeightValues)) %>%
  dplyr::mutate(fft_ccf_df = cmcR::comparison_fft_ccf(cellHeightValues = cellHeightValues_replaced,
                                         regionHeightValues = regionHeightValues_replaced))  %>%
  dplyr::mutate(pairwiseCompCor = cmcR::comparison_cor(cellHeightValues,regionHeightValues,fft_ccf_df)) %>%
  tidyr::unnest(fft_ccf_df) %>%
  dplyr::mutate(theta = 0)

#Use comparison_allTogether and compare results

cellTibble2 <- cmcR::comparison_allTogether(x3p1,x3p2,
                                            theta = 0,numCells = 64,
                                            maxMissingProp = .85) %>%
  dplyr::select(-cellIndex)

testthat::test_that("comparison_ functions work as expected", {

  testthat::has_names(cellTibble,
                      c("cellIndex","cellHeightValues","regionHeightValues","cellPropMissing","regionPropMissing","cellHeightValues_replaced","regionHeightValues_replaced","fft_ccf","x","y","pairwiseCompCor"))
  testthat::expect_s3_class(cellTibble,c("tbl_df","tbl","data.frame"))

  testthat::expect_true(identical(dplyr::select(cellTibble,c("x","y","fft_ccf","pairwiseCompCor","theta")),
                                  cellTibble2))

  #Add more "expect failure" tests?
})
