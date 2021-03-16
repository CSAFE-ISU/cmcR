`%>%` <- dplyr::`%>%`

x3p1 <- x3ptools::read_x3p(tmpfile1)

if(!exists("skipPreprocess")){

  x3p1_raw <- x3p1

  x3p1_raw_dim <- dim(x3p1$surface.matrix)
  x3p1_raw_missing <- sum(is.na(x3p1$surface.matrix))

  x3p1 <- x3p1 %>%
    cmcR::preProcess_crop(region = "exterior",
                          radiusOffset = -30)

  #cropping exterior should reduce dimension, but remove NAs on exterior of scan
  x3p1_extCrop_dim <- dim(x3p1$surface.matrix)
  x3p1_extCrop_missing <- sum(is.na(x3p1$surface.matrix))

  #remove firing pin observations
  x3p1 <- x3p1 %>%
    cmcR::preProcess_crop(region = "interior",
                          radiusOffset = 200)

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

  testthat::test_that("preProcess_ functions work as expected", {
    if(x3p1$cmcR.info$skipPreprocess == 1){
      testthat::skip()
    }

    testthat::expect_true(all(x3p1_extCrop_dim <= x3p1_raw_dim))
    testthat::expect_true(x3p1_extCrop_missing <= x3p1_raw_missing)

    #croppint interior should not change dimension, but should introduct more NAs
    testthat::expect_equal(x3p1_intCrop_dim, x3p1_extCrop_dim)
    testthat::expect_true(x3p1_intCrop_missing >=  x3p1_extCrop_missing)

    #de-trending shouldn't affect which indices contain NAs or observed values
    testthat::expect_true(all(which(is.na(x3p1_meanDeTrend$surface.matrix)) == x3p1_preDeTrend_missingIndices))
    testthat::expect_true(all(which(is.na(x3p1_medianDeTrend$surface.matrix)) == x3p1_preDeTrend_missingIndices))

    #de-trending should reduce (large scale) variability in height values
    testthat::expect_true(x3p1_meanDeTrend_var <= x3p1_preDeTrend_var)
    testthat::expect_true(x3p1_medianDeTrend_var <= x3p1_preDeTrend_var)

    #any Gaussian filter should attenuate affect of certain frequencies (and thus
    #reduce variability)
    testthat::expect_true(postFilterVar_lp < x3p1_medianDeTrend_var)
    testthat::expect_true(postFilterVar_hp < x3p1_medianDeTrend_var)
    testthat::expect_true(postFilterVar_bp < x3p1_medianDeTrend_var)

    #final downsampled scan should have exactly half the dimension of original
    #(since original has even dimension)
    testthat::expect_equal(dim(x3p1$surface.matrix), x3p1_intCrop_dim/2)

    #Add more "expect failure" tests?
  })


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

  testthat::test_that("Legacy preProcess_ functions work as expected", {
    if(x3p1$cmcR.info$skipPreprocess == 1){
      testthat::skip()
    }

    #applying RANSAC method shouldn't affect dimension, but should introduce more
    #NAs
    testthat::expect_true(all(x3p1_ransacLeveled_dim == x3p1_downSampled_dim))
    testthat::expect_true(x3p1_ransacLeveled_missing >= x3p1_downSampled_missing)

    #removing firing pin circle should introduce more NAs and reduce variability
    testthat::expect_true(all(x3p1_fpCircleRemoved_dim == x3p1_downSampled_dim))
    testthat::expect_true(x3p1_fpCircleRemoved_missing >= x3p1_ransacLeveled_missing)
    testthat::expect_true(x3p1_fpCircleRemoved_var <= x3p1_ransacLeveled_var)

    #Add more "expect failure" tests?
  })


}
