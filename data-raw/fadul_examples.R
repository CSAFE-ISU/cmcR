# code to prepare `fadul_examples` dataset goes here

download.file("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d",
              destfile = "data/fadul1-1.x3p",
              mode = "wb")

download.file("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8",
              destfile = "data/fadul1-2.x3p",
              mode = "wb")

fadul1.1_raw <- x3ptools::read_x3p("data/fadul1-1.x3p")

fadul1.2_raw <- x3ptools::read_x3p("data/fadul1-2.x3p")

fadul1.1_processed <- fadul1.1_raw %>%
  cmcR::preProcess_cropExterior(radiusOffset = -30,
                                agg_function = median) %>%
  cmcR::preProcess_filterInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                               tau = .5,
                               method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

fadul1.2_processed <- fadul1.2_raw %>%
  cmcR::preProcess_cropExterior(radiusOffset = -30,
                                agg_function = median) %>%
  cmcR::preProcess_filterInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend(statistic = "quantile",
                               tau = .5,
                               method = "fn") %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

usethis::use_data(fadul1.1_processed,
                  fadul1.2_processed,
                  overwrite = TRUE)

try(file.remove("data/fadul1-1.x3p"))
try(file.remove("data/fadul1-2.x3p"))