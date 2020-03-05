#' @name ccfMap
#'
#' @keywords internal

ccfMap <- function(mat1,mat2){
  #Test that mat1,mat2 contain no NAs
  if(any(is.na(as.vector(mat1))) | any(is.na(as.vector(mat2)))){
    mat1[is.na(mat1)] <- 0
    mat2[is.na(mat2)] <- 0
  }

  ccfMat <- filterViaFFT(mat1,mat2) / (sqrt(sum(mat1^2)) * sqrt(sum(mat2^2)))

  return(Re(ccfMat))

  #this function should only be used on matrices that have been pre-processed - at the very least, had NAs replaced by another value prior to use.

  #Call "mat1" whichever the the smaller mat is and "mat2" the larger
  # matSmall <- list(mat1,mat2)[[all(dim(mat1) > dim(mat2)) + 1]] #will assign smaller of two matrices to mat1
  # matBig <- list(mat1,mat2)[[all(dim(mat1) < dim(mat2)) + 1]] #will assign larger of two matrices to mat2
  #
  # # size of full filter
  # dim(matSmall) <- dim(matSmall)
  # dim(matBig) <- dim(matBig)
  #
  # dimPadded <- dim(matSmall) + dim(matBig) - 1
  #
  # # pad images with 0 so that we do not have circular issues with FFT
  # padMatSmall <- matrix(0, nrow = dimPadded[1], ncol = dimPadded[2])
  # padMatBig <- matrix(0, nrow = dimPadded[1], ncol = dimPadded[2])
  #
  # padMatSmall[1:dim(matSmall)[1], 1:dim(matSmall)[2]] <- matSmall
  # padMatBig[1:dim(matBig)[1], 1:dim(matBig)[2]] <- matBig
  #
  # # Filter in frequency domain
  # ccfMat <- fftw::FFT(fftw::FFT(padMatSmall)*Conj(fftw::FFT(padMatBig)), inverse = TRUE)/(prod(dimPadded))
  #
  # ccfMat <- matrix(ccfMat,nrow = dimPadded[1],ncol = dimPadded[2])
  #
  # ccfMat <- circshift(fftshift(ccfMat), round2((dim(matBig) - dim(matSmall))/2, 0))
  #
  # # halfDimSmall <- cmcR:::round2((dim(matBig) - dim(matSmall))/2, 0)
  #
  # ccfMat_validRegion <- ccfMat
  #
  # # ccfMat_validRegion <- ccfMat[halfDimSmall[1]:(halfDimSmall[1] + dim(matBig)[1] - 1),
  # #                              halfDimSmall[2]:(halfDimSmall[2] + dim(matBig)[2] - 1)]
  #
  # if (all.equal(c(Im(ccfMat_validRegion)), rep(0, prod(dim(ccfMat_validRegion)))) == FALSE) {
  #   stop("Non-zero imaginary part")
  # }

  # ccfMat <- ccfMat_validRegion / (sqrt(sum(mat1^2)) * sqrt(sum(mat2^2)))
}

#' @name ccfMapPlot
#'
#' @export

#TODO: add theta and cellID to the plot table

ccfMapPlot <- function(mat1,
                       mat2,
                       params,
                       returnGrob = FALSE){

  ccfMat <- cmcR:::ccfMap(mat1,mat2)

  ccfDF <- ccfMat %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(dx = x - max(x)/2,
                  dy = y - max(y)/2) %>%
    dplyr::rename(CCF = value)

  ccfMaxInfo <- ccfDF %>%
    dplyr::filter(CCF == max(CCF)) %>%
    dplyr::mutate(CCF = round(CCF,3))

  mat1 <- (mat1 - mean(mat1,na.rm = TRUE))/(sd(mat1,na.rm = TRUE))
  mat2 <- (mat2 - mean(mat2,na.rm = TRUE))/(sd(mat2,na.rm = TRUE))

  mat1Plot <- mat1 %>%
    t() %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    # mutate(x = x + floor(max(nrow(mat1),nrow(mat2))/4),
    #        y = y + floor(max(ncol(mat1),ncol(mat2))/4)) %>%
    ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient2(low = "grey0",
                                  mid = "grey50",
                                  high = "grey100") +
    ggplot2::coord_fixed(xlim = c(0,max(nrow(mat1),nrow(mat2))),
                         ylim = c(0,max(ncol(mat1),ncol(mat2)))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())

  mat2Plot <- mat2 %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    imager::as.cimg() %>%
    as.data.frame() %>%
    ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient2(low = "grey0",
                                  mid = "grey50",
                                  high = "grey100",
                                  midpoint = 0) +
    ggplot2::coord_fixed(xlim = c(0,max(nrow(mat1),nrow(mat2))),
                         ylim = c(0,max(ncol(mat1),ncol(mat2)))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::geom_tile(ggplot2::aes(
      x = max(x)/2 - ccfMaxInfo$dx,
      y = max(y)/2 - ccfMaxInfo$dy,
      width = ncol(mat1),
      height = nrow(mat1)),
      alpha = 0,
      colour = "orange")

  ccfPlot <- ccfDF %>%
    dplyr::mutate(dx = rev(dx),
                  dy = rev(dy)) %>%
    ggplot2::ggplot(ggplot2::aes(x = dx,y = dy,fill = CCF)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(low = "purple",
                                  mid = "white",
                                  high = "orange",
                                  midpoint = 0) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  ccfMaxSummary <- ccfMaxInfo %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::mutate(dx = -dx,
                  dy = -dy) %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    gridExtra::tableGrob(rows = c("CCFmax","dx","dy"),
                         cols = "CCFmax.Summary")

  layoutMat <- matrix(c(1,2,3,4),ncol = 2,byrow = TRUE)

  gridPlot <- gridExtra::arrangeGrob(mat1Plot,
                                     mat2Plot,
                                     ccfPlot,
                                     ccfMaxSummary,
                                     layout_matrix = layoutMat)

  plot(gridPlot)

  if(returnGrob){
    return(gridPlot)
  }
}