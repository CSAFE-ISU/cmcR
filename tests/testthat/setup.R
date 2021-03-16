`%>%` <- dplyr::`%>%`

utils::data("fadul1.1_processed")
utils::data("fadul1.2_processed")

x3p1 <- fadul1.1_processed
x3p2 <- fadul1.2_processed

tmpfile1 <<- tempfile(fileext = ".x3p")
tmpfile2 <<- tempfile(fileext = ".x3p")

x3ptools::write_x3p(x3p1,file = tmpfile1)
x3ptools::write_x3p(x3p2,file = tmpfile2)