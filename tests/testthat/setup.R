`%>%` <- dplyr::`%>%`

x3p1_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
x3p2_url <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8"

tryCatch({

  suppressWarnings({
    x3p1 <- x3ptools::read_x3p(file = x3p1_url)
    x3p2 <- x3ptools::read_x3p(file = x3p2_url)
  })

  tmpfile1 <<- tempfile(fileext = ".x3p")
  tmpfile2 <<- tempfile(fileext = ".x3p")

  x3ptools::write_x3p(x3p1,file = tmpfile1)
  x3ptools::write_x3p(x3p2,file = tmpfile2)

},error = function(e){

  load("../../data/fadul1.1_processed.rda")
  load("../../data/fadul1.2_processed.rda")

  x3p1 <- fadul1.1_processed
  x3p2 <- fadul1.2_processed

  skipPreprocess <<- 1

  x3p1$cmcR.info$skipPreprocess <- 1
  x3p2$cmcR.info$skipPreprocess <- 1

  tmpfile1 <<- tempfile(fileext = ".x3p")
  tmpfile2 <<- tempfile(fileext = ".x3p")

  x3ptools::write_x3p(x3p1,file = tmpfile1)
  x3ptools::write_x3p(x3p2,file = tmpfile2)
})