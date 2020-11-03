`%>%` <- dplyr::`%>%`

x3p1_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
x3p2_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"

x3p1 <- x3ptools::x3p_read(file = x3p1_url) %>%
  x3ptools::sample_x3p(m = 8)

x3p1 <- x3p1 %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_filterInterior(region = "interior",
                                  radiusOffset = 200) %>%
  cmcR::preProcess_removeTrend()

x3p2 <- x3ptools::x3p_read(file = x3p2_url) %>%
  x3ptools::sample_x3p(m = 8)

x3p2 <- x3p2 %>%
  cmcR::preProcess_crop(region = "exterior",
                        radiusOffset = -30) %>%
  cmcR::preProcess_filterInterior(region = "interior",
                                  radiusOffset = 200)
  cmcR::preProcess_removeTrend()

tmpfile1 <- tempfile(fileext = ".x3p")
tmpfile2 <- tempfile(fileext = ".x3p")

x3ptools::write_x3p(x3p1,file = tmpfile1)
x3ptools::write_x3p(x3p2,file = tmpfile2)
