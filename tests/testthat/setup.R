`%>%` <- dplyr::`%>%`

set.seed(4132020)

x3p1_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
x3p2_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"

x3p1 <- x3ptools::read_x3p(file = x3p1_url) %>%
  x3ptools::sample_x3p(m = 8)

x3p1$surface.matrix <- x3p1$surface.matrix %>%
  cmcR::preProcess_ransac() %>%
  cmcR::preProcess_levelBF(useResiduals = TRUE)  %>%
  cmcR::preProcess_cropWS(croppingThresh =  1)  %>%
  cmcR::preProcess_removeFPCircle(aggregation_function = function(x,na.rm){min(x,na.rm = na.rm) - 15})

x3p1$header.info$sizeY <- ncol(x3p1$surface.matrix)
x3p1$header.info$sizeX <- nrow(x3p1$surface.matrix)

x3p2 <- x3ptools::read_x3p(file = x3p2_url) %>%
  x3ptools::sample_x3p(m = 8)

x3p2$surface.matrix <- x3p2$surface.matrix %>%
  cmcR::preProcess_ransac() %>%
  cmcR::preProcess_levelBF(useResiduals = TRUE)  %>%
  cmcR::preProcess_cropWS(croppingThresh =  1)  %>%
  cmcR::preProcess_removeFPCircle(aggregation_function = function(x,na.rm){min(x,na.rm = na.rm) - 15})

x3p2$header.info$sizeY <- ncol(x3p2$surface.matrix)
x3p2$header.info$sizeX <- nrow(x3p2$surface.matrix)

tmpfile1 <- tempfile(fileext = ".x3p")
tmpfile2 <- tempfile(fileext = ".x3p")

x3ptools::write_x3p(x3p1,file = tmpfile1)
x3ptools::write_x3p(x3p2,file = tmpfile2)
