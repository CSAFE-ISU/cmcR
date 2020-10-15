# code to prepare `fadul_examples` dataset goes here

fadul1.1_raw <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d")

fadul1.1_processed <- fadul1.1_raw %>%
  cmcR::cropBFExterior(radiusOffset = -30,
                       agg_function = median) %>%
  cmcR::filterBFInterior(radiusOffset = 200) %>%
  cmcR::levelByConditionalStatistic(statistic = "quantile",
                                    tau = .5,
                                    method = "fn") %>%
  x3ptools::sample_x3p(m = 2)

fadul1.1_processed$surface.matrix <- cmcR::preProcess_gaussFilter(surfaceMat = fadul1.1_processed$surface.matrix,
                                                                  res = fadul1.1_processed$header.info$incrementY)

fadul1.2_raw <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8")

fadul1.2_processed <- fadul1.2_raw %>%
  cmcR::cropBFExterior(radiusOffset = -30,
                       agg_function = median) %>%
  cmcR::filterBFInterior(radiusOffset = 200) %>%
  cmcR::levelByConditionalStatistic(statistic = "quantile",
                                    tau = .5,
                                    method = "fn") %>%
  x3ptools::sample_x3p(m = 2)

fadul1.2_processed$surface.matrix <- cmcR::preProcess_gaussFilter(surfaceMat = fadul1.2_processed$surface.matrix,
                                                                  res = fadul1.2_processed$header.info$incrementY)

fadul2.1_raw <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/8ae0b86d-210a-41fd-ad75-8212f9522f96")

fadul2.1_processed <- fadul2.1_raw %>%
  cmcR::cropBFExterior(radiusOffset = -30,
                       agg_function = median) %>%
  cmcR::filterBFInterior(radiusOffset = 200) %>%
  cmcR::levelByConditionalStatistic(statistic = "quantile",
                                    tau = .5,
                                    method = "fn") %>%
  x3ptools::sample_x3p(m = 2)

fadul2.1_processed$surface.matrix <- cmcR::preProcess_gaussFilter(surfaceMat = fadul2.1_processed$surface.matrix,
                                                                  res = fadul2.1_processed$header.info$incrementY)

usethis::use_data(fadul1.1_processed,fadul1.2_processed, fadul2.1_processed,overwrite = TRUE)
