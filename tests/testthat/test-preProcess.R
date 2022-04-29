`%>%` <- dplyr::`%>%`

testthat::test_that("preProcess_ functions work as expected", {

  testthat::skip_on_cran()
  testthat::skip_if_offline()

  if(identical(Sys.getenv("NOT_CRAN"), "true")){

    x3p1_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
    # x3p2_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"

    suppressWarnings({
      x3p1 <- x3ptools::read_x3p(file = x3p1_url)
      # x3p2 <- x3ptools::read_x3p(file = x3p2_url)
    })

    x3p1_raw <- x3p1

    x3p1_raw_dim <- dim(x3p1$surface.matrix)
    x3p1_raw_missing <- sum(is.na(x3p1$surface.matrix))

    #preProcess_ransacLevel reduces the variability in scan
    set.seed(4292022)

    x3p1_ransac <- x3p1_raw %>%
      x3ptools::sample_x3p(m = 4) %>%
      cmcR::preProcess_ransacLevel()

    preRansacVar <- var(c(x3p1_raw$surface.matrix),na.rm = TRUE)*1e12
    postRansacVar <- var(c(x3p1_ransac$surface.matrix),na.rm = TRUE)*1e12

    x3p1 <- x3p1 %>%
      cmcR::preProcess_crop(region = "exterior",
                            offset = -50)

    #cropping exterior should reduce dimension, but remove NAs on exterior of scan
    x3p1_extCrop_dim <- dim(x3p1$surface.matrix)
    x3p1_extCrop_missing <- sum(is.na(x3p1$surface.matrix))

    #cropping with robust estimate shouldn't affect the scan dimensions
    x3p1_cropped <- x3p1 %>%
      x3ptools::sample_x3p() %>%
      cmcR:::preProcess_cropWS(croppingThresh = .5,robust = TRUE)

    x3p1_center <- x3p1 %>%
      x3ptools::sample_x3p(2) %>%
      cmcR:::fpCenterCalc()

    #remove firing pin observations
    x3p1 <- x3p1 %>%
      cmcR::preProcess_crop(region = "interior")

    x3p1_intCrop_dim <- dim(x3p1$surface.matrix)
    x3p1_intCrop_missing <- sum(is.na(x3p1$surface.matrix))

    x3p1_preDeTrend_var <- var(x3p1$surface.matrix[!is.na(x3p1$surface.matrix)])*1e12

    x3p1_preDeTrend_missingIndices <- which(is.na(x3p1$surface.matrix))

    #remove conditional/median
    x3p1_meanDeTrend <- x3p1 %>%
      cmcR::preProcess_removeTrend(statistic = "mean")

    x3p1_medianDeTrend <- x3p1 %>%
      cmcR::preProcess_removeTrend(statistic = "quantile",
                                   tau = .5,
                                   method = "fn")

    x3p1_meanDeTrend_var <- var(x3p1_meanDeTrend$surface.matrix[!is.na(x3p1$surface.matrix)])*1e12
    x3p1_medianDeTrend_var <- var(x3p1_medianDeTrend$surface.matrix[!is.na(x3p1$surface.matrix)])*1e12

    #Pass Gaussian filter over scan
    x3p1_lp <- x3p1_medianDeTrend %>%
      cmcR::preProcess_gaussFilter(wavelength = c(16),
                                   filtertype = "lp")
    x3p1_hp <- x3p1_medianDeTrend %>%
      cmcR::preProcess_gaussFilter(wavelength = c(500),
                                   filtertype = "hp")
    x3p1_bp <- x3p1_medianDeTrend %>%
      cmcR::preProcess_gaussFilter(wavelength = c(16,500),
                                   filtertype = "bp")

    postFilterVar_lp <- var(x3p1_lp$surface.matrix[!is.na(x3p1_lp$surface.matrix)])*1e12
    postFilterVar_hp <- var(x3p1_hp$surface.matrix[!is.na(x3p1_hp$surface.matrix)])*1e12
    postFilterVar_bp <- var(x3p1_bp$surface.matrix[!is.na(x3p1_bp$surface.matrix)])*1e12

    x3p1 <- x3p1 %>%
      x3ptools::sample_x3p()

    x3p1_eroded <- x3p1 %>%
      cmcR::preProcess_erode(region = "interior") %>%
      cmcR::preProcess_erode(region = "exterior") %>%
      cmcR:::preProcess_cropWS(croppingProp = .5)

    message("Skipping downloading of remote scan on CRAN.")
  }

  # if(x3p1$cmcR.info$skipPreprocess == 1){
  #   testthat::skip()
  # }

  testthat::expect_true(all(x3p1_extCrop_dim <= x3p1_raw_dim))
  testthat::expect_true(x3p1_extCrop_missing <= x3p1_raw_missing)

    # if(x3p1$cmcR.info$skipPreprocess == 1){
    #   testthat::skip()
    # }

    # ransac-leveled scan should reduce in extremes and variance of height
    # values
    testthat::expect_true(postRansacVar < preRansacVar &
                            abs(min(c(x3p1_raw$surface.matrix),na.rm = TRUE)) < abs(min(c(x3p1_ransac$surface.matrix),na.rm = TRUE)) &
                            abs(max(c(x3p1_raw$surface.matrix),na.rm = TRUE)) > abs(max(c(x3p1_ransac$surface.matrix),na.rm = TRUE)))

    testthat::expect_true(all(x3p1_extCrop_dim <= x3p1_raw_dim))
    testthat::expect_true(x3p1_extCrop_missing <= x3p1_raw_missing)

    #cropping interior should not change dimension, but should introduct more NAs
    testthat::expect_equal(x3p1_intCrop_dim, x3p1_extCrop_dim)
    testthat::expect_true(x3p1_intCrop_missing >=  x3p1_extCrop_missing)

    #de-trending shouldn't affect which indices contain NAs or observed values
    testthat::expect_true(all(which(is.na(x3p1_meanDeTrend$surface.matrix)) == x3p1_preDeTrend_missingIndices))
    testthat::expect_true(all(which(is.na(x3p1_medianDeTrend$surface.matrix)) == x3p1_preDeTrend_missingIndices))

    #de-trending should reduce (large scale) variability in height values
    testthat::expect_true(x3p1_meanDeTrend_var <= x3p1_preDeTrend_var)
    testthat::expect_true(x3p1_medianDeTrend_var <= x3p1_preDeTrend_var)

    # robust preProcess_cropWS shouldn't affect dimensions after cropping
    testthat::expect_true(all(dim(x3p1$surface.matrix) == dim(x3p1_cropped$surface.matrix)))

    # fpCenterCalc should work as expected
    testthat::expect_true(all(x3p1_center == c(215,259)))

    #Now check older preProcess functions:

    x3p1_downSampled <- x3p1_raw %>%
      x3ptools::sample_x3p()

    x3p1_downSampled_dim <- dim(x3p1_downSampled$surface.matrix)
    x3p1_downSampled_missing <- sum(is.na(x3p1_downSampled$surface.matrix))

    x3p1_downSampled <- x3p1_downSampled %>%
      cmcR::preProcess_ransacLevel()

    x3p1_ransacLeveled_dim <- dim(x3p1_downSampled$surface.matrix)
    x3p1_ransacLeveled_missing <- sum(is.na(x3p1_downSampled$surface.matrix))

    x3p1_ransacLeveled_var <- var(x3p1_downSampled$surface.matrix[!is.na(x3p1_downSampled$surface.matrix)])

    x3p1_fpCircleRemoved <- x3p1_downSampled %>%
      cmcR::preProcess_removeFPCircle()

    x3p1_fpCircleRemoved_dim <- dim(x3p1_fpCircleRemoved$surface.matrix)
    x3p1_fpCircleRemoved_missing <- sum(is.na(x3p1_fpCircleRemoved$surface.matrix))
    x3p1_fpCircleRemoved_var <- var(x3p1_fpCircleRemoved$surface.matrix[!is.na(x3p1_fpCircleRemoved$surface.matrix)])

    # if(x3p1$cmcR.info$skipPreprocess == 1){
    #   testthat::skip()
    # }

    #applying RANSAC method shouldn't affect dimension, but should introduce more
    #NAs
    testthat::expect_true(all(x3p1_ransacLeveled_dim == x3p1_downSampled_dim))
    testthat::expect_true(x3p1_ransacLeveled_missing >= x3p1_downSampled_missing)

    #applying RANSAC method shouldn't affect dimension, but should introduce more
    #NAs
    testthat::expect_true(all(x3p1_ransacLeveled_dim == x3p1_downSampled_dim))
    testthat::expect_true(x3p1_ransacLeveled_missing >= x3p1_downSampled_missing)

    #removing firing pin circle should introduce more NAs and reduce variability
    testthat::expect_true(all(x3p1_fpCircleRemoved_dim == x3p1_downSampled_dim))
    testthat::expect_true(x3p1_fpCircleRemoved_missing >= x3p1_ransacLeveled_missing)
    testthat::expect_true(x3p1_fpCircleRemoved_var <= x3p1_ransacLeveled_var)

    #preProcess_erode + preProcess_cropWS reduces the dimension of the surface matrix as expected
    testthat::expect_true(all(dim(x3p1_eroded$surface.matrix) == c(423,429)))

    #Add more "expect failure" tests?
})
