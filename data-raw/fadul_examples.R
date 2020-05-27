# code to prepare `fadul_examples` dataset goes here
set.seed(4132020)

fadul1.1 <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d") %>%
  x3ptools::sample_x3p(m = 8)

fadul1.1$surface.matrix <- fadul1.1$surface.matrix %>%
  cmcR::preProcess_ransac() %>%
  cmcR::preProcess_levelBF(useResiduals = TRUE)  %>%
  cmcR::preProcess_cropWS(croppingThresh =  1)  %>%
  cmcR::preProcess_removeFPCircle(aggregation_function = function(x,na.rm){min(x,na.rm = na.rm) - 15}) %>%
  cmcR::preProcess_gaussFilter(res = fadul1.1$header.info$incrementY)

fadul1.2 <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8") %>%
  x3ptools::sample_x3p(m = 8)

fadul1.2$surface.matrix <- fadul1.2$surface.matrix %>%
  cmcR::preProcess_ransac() %>%
  cmcR::preProcess_levelBF(useResiduals = TRUE)  %>%
  cmcR::preProcess_cropWS(croppingThresh =  1)  %>%
  cmcR::preProcess_removeFPCircle(aggregation_function = function(x,na.rm){min(x,na.rm = na.rm) - 15}) %>%
  cmcR::preProcess_gaussFilter(res = fadul1.2$header.info$incrementY)

usethis::use_data(fadul1.1,fadul1.2, overwrite = TRUE)
