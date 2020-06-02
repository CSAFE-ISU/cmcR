#' @name gaussianKernel
#'
#' @seealso https://rdrr.io/bioc/EBImage/src/R/morphology.R
#'
#' @keywords internal

gaussianKernel <- function(size,sig){
  x = seq(-(size-1)/2, (size-1)/2, length=size)
  x = matrix(x, size, size)
  z = exp(- (x^2 + t(x)^2) / (2*sig^2))
  z = z / sum(z)

  return(z)
}

#' @name gaussianFilter
#'
#' @seealso https://www.mathworks.com/matlabcentral/fileexchange/61003-filt2-2d-geospatial-data-filter?focused=7181587&tab=example
#' @keywords internal

utils::globalVariables(".")

gaussianFilter <- function(surfaceMat,
                           res,
                           wavelength,
                           filtertype = "bp",
                           filterSize = NULL){

  if(filtertype %in% c("lp","hp")){
    sig <- (wavelength/res)/(2*pi)

  if(is.null(filterSize)){
    filterSize <- {2*ceiling(2.6*sig) + 1} %>%
      magrittr::add(. %% 2)
  }

    kern <- gaussianKernel(size = filterSize,
                           sig = sig)

    surfaceMatP <- surfaceMat #may need to pad rows/cols with 0s to make even dims

    if((filterSize - dim(surfaceMat)[1]) %% 2 == 1){
      surfaceMatP <- surfaceMatP %>%
        rbind(rep(0,times = ncol(.)))
    }
    if((filterSize - dim(surfaceMat)[2]) %% 2 == 1){
      surfaceMatP <- surfaceMatP %>%
        cbind(rep(0,times = nrow(.)))
    }

    dimPadded <- dim(surfaceMatP) + filterSize

    imPadded <- surfaceMatP %>%
      imager::as.cimg() %>%
      imager::pad(nPix = dimPadded[1] - dim(surfaceMatP)[1],
                  axes = "y",
                  pos = 0) %>%
      imager::pad(nPix = dimPadded[2] - dim(surfaceMatP)[2],
                  axes = "x",
                  pos = 0) %>%
      as.matrix()

    kernPadded <- kern %>%
      imager::as.cimg() %>%
      imager::pad(nPix = {dimPadded[1] - filterSize} %>%
                    magrittr::add(. %% 2),
                  axes = "x",
                  pos = 0) %>%
      imager::pad(nPix = {dimPadded[2] - filterSize} %>%
                    magrittr::add(. %% 2),
                  axes = "y",
                  pos = 0) %>%
      as.matrix()
  }

  if(filtertype == "lp"){
    imFiltered <- imPadded %>%
      fft() %>%
      magrittr::multiply_by(fft(kernPadded)) %>%
      fft(inverse = TRUE) %>%
      magrittr::divide_by(prod(dimPadded)) %>%
      fftshift()

    imFiltered <- Re(imFiltered) %>%
      imager::as.cimg() %>%
      imager::crop.borders(nx = {dimPadded[1] - dim(surfaceMatP)[1]} %>%
                             magrittr::add(-1*(. %% 2)) %>%
                             magrittr::divide_by(2),
                           ny = {dimPadded[2] - dim(surfaceMatP)[2]} %>%
                             magrittr::add(-1*(. %% 2)) %>%
                             magrittr::divide_by(2)) %>%
      as.matrix()
  }
  else if(filtertype == "hp"){
    imFilteredLP <- imPadded %>%
      fft() %>%
      magrittr::multiply_by(fft(kernPadded)) %>%
      fft(inverse = TRUE) %>%
      magrittr::divide_by(prod(dimPadded)) %>%
      fftshift()

    imFiltered <- imPadded - imFilteredLP

    imFiltered <- Re(imFiltered) %>%
      imager::as.cimg() %>%
      imager::crop.borders(nx = {dimPadded[1] - dim(surfaceMatP)[1]} %>%
                             magrittr::add(-1*(. %% 2)) %>%
                             magrittr::divide_by(2),
                           ny = {dimPadded[2] - dim(surfaceMatP)[2]} %>%
                             magrittr::add(-1*(. %% 2)) %>%
                             magrittr::divide_by(2)) %>%
      as.matrix()
  }
  else if(filtertype == "bp"){
    imFiltered <- gaussianFilter(surfaceMat = surfaceMat %>%
                                   gaussianFilter(res = res,
                                                  wavelength = max(wavelength),
                                                  filtertype = "hp"),
                                 res = res,
                                 wavelength = min(wavelength),
                                 filtertype = "lp")
  }

  #remove extra row/col that was artificially added to make dim of surfaceMat even, if needed
  if((dim(imFiltered)[1] - dim(surfaceMat)[1]) %% 2 == 1){
    imFiltered <- imFiltered[1:(nrow(imFiltered) - 1),1:ncol(imFiltered)]
  }
  if((dim(imFiltered)[2] - dim(surfaceMat)[2]) %% 2 == 1){
    imFiltered <- imFiltered[1:nrow(imFiltered),1:(ncol(imFiltered) - 1)]
  }

  return(imFiltered)
}