x3p1_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"

tmp <- selectBFImpression_sample_x3p(x3p1_url,croppingThresh = -1)

test_that("selectBFImpression_ functions work as expected", {
  testthat::expect_equal(tmp$x3p$header.info$incrementY,6.25e-06)
  testthat::expect_equal(dim(tmp$x3p$surface.matrix),c(616,612))
  testthat::expect_equal(tmp$params$croppingThresh,-1)

  testthat::expect_error(selectBFImpression_sample_x3p(file.path("does","not","exist.x3p")))
})
