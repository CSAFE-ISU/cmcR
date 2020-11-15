# These functions were copied from the cartridges3D R package currently (as of
# 6/13/2020) available on GitHub: https://github.com/xhtai/cartridges3D.

#These functions were copied because they are internal to both the cartridges3D
#and cmcR packages and referring to another package's internal functions is not
#SOP in R package development (devtools::check warnings when I had it
#differently)


# Calculates cross-correlation between two matrices using FFTs
#
# @name filterViaFFT
#
#
# @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
# @keywords internal
#
#' @importFrom stats fft

filterViaFFT <- function(A, B) {
  # size of full filter
  m <- dim(A)
  n <- dim(B)
  x <- m + n - 1

  # pad images with 0 so that we do not have circular issues with FFT
  padA <- matrix(0, nrow = x[1], ncol = x[2])
  padB <- matrix(0, nrow = x[1], ncol = x[2])
  padA[1:m[1], 1:m[2]] <- A
  padB[1:n[1], 1:n[2]] <- B

  # Filter in frequency domain
  C <- fft(fft(padA)*Conj(fft(padB)), inverse = TRUE)/(prod(x))

  C <- matrix(C,nrow = x[1],ncol = x[2])

  C <- circshift(fftshift(C), round2((n - m)/2, 0))

  half_m <- round2(m/2, 0)
  C <- C[half_m[1]:(half_m[1] + n[1] - 1), half_m[2]:(half_m[2] + n[2] - 1)]
  if (all.equal(c(Im(C)), rep(0, prod(dim(C)))) == FALSE) {
    stop("Non-zero imaginary part")
  }
  return(Re(C))
}

# Copies behavior of round() function in MATLAB
#
# @name round2
# @seealso http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
# @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
# @keywords internal

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

# Shifts CCF matrix output so that Nyquist frequency is in middle of matrix
#
# @name fftshift
# @seealso http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
# @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
# @keywords internal

fftshift <- function(input_matrix) {

  input_matrix <- as.matrix(input_matrix)

  rows <- dim(input_matrix)[1]
  cols <- dim(input_matrix)[2]

  swap_up_down <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
  }

  swap_left_right <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
  }

  input_matrix <- swap_up_down(input_matrix)
  return(swap_left_right(input_matrix))

}

# Performs a circular shift based on periodic boundary conditions in a vector
#
# @name circshift
# @seealso http://stackoverflow.com/questions/18791212/equivalent-to-numpy-roll-in-r
# @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
# @keywords internal
#
#' @importFrom utils head tail

circshift <- function(x, vec) {
  dimx <- dim(x)
  # row first
  if (vec[1] != 0) {
    # out <- rbind(x[(dimx[1] - vec[1] + 1):dimx[1], ], x[1:(dimx[1] - vec[1]), ])
    tmp <- c(t(x))
    x <- matrix(c( tail(tmp, vec[1]*dimx[2]) , head(tmp, -vec[1]*dimx[2]) ), byrow = TRUE, nrow = dimx[1])
  }
  # col
  if (vec[2] != 0) {
    tmp <- c(x)
    x <- matrix(c( tail(tmp, vec[2]*dimx[1]) , head(tmp, -vec[2]*dimx[1]) ), nrow = dimx[1])
  }
  return(x)
}

# pairwise_ccf <- function(a, b, offsetx = 0, offsety = 0) {
#   if ("x3p" %in% class(a)) a <- a$surface.matrix
#   if ("x3p" %in% class(b)) b <- b$surface.matrix
#
#   aadj <- b_pad <- matrix(NA,
#                           nrow = nrow(a) + nrow(b),
#                           ncol = ncol(a) + ncol(b))
#
#   hdim <- (dim(a) + dim(b))/4
#
#   b_pad[(hdim[1] + 1):(hdim[1] + dim(b)[1]),
#         (hdim[2] + 1):(hdim[2] + dim(b)[2])] <- b
#
#   idx <- (offsetx + 1):(offsetx + ncol(a))
#   idy <- (offsety + 1):(offsety + nrow(a))
#   aadj[idy,idx] <- a
#
#   if(all(is.na(aadj)) | all(is.na(b_pad))){
#     return(NA)
#   }
#   else{
#     suppressWarnings(
#       suppressMessages(
#         list("ccf" = cor(as.numeric(aadj), as.numeric(b_pad), use = "pairwise.complete.obs"),
#              "nonMissing" = sum(!is.na(as.numeric(aadj)) & !is.na(b_pad)))
#
#       )
#     )
#   }
# }

# Computes the location of the maximum CCF value in a CCF map between two
# matrices
#
# @name ccfComparison
#
# @param ccfMethod implements 3 different methods to calculate the CCF -- all which yield differing CCF values
# @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
# @keywords internal

ccfComparison <- function(mat1, mat2, ccfMethod = "fft") {
  stopifnot(ccfMethod %in% c("bruteForceReweighted","imager","fft"))

  # if(ccfMethod == "bruteForceReweighted"){
  #   resp <- expand.grid(offsetx = 1:(max(ncol(mat1),ncol(mat2))),
  #                       offsety = 1:(max(nrow(mat1),nrow(mat2)))) %>%
  #     purrr::pmap_dfr(~ {
  #       pwiseCCF <- pairwise_ccf(mat1,mat2,offsetx = ..1,offsety = ..2)
  #
  #       data.frame(offsetx = ..1,
  #                  offsety = ..2,
  #                  ccf = pwiseCCF$ccf,
  #                  nonMissing = pwiseCCF$nonMissing)
  #     }) %>%
  #     dplyr::mutate(ccfReweighted = nonMissing*ccf/max(nonMissing),
  #                   offsetx = offsetx - max(offsetx)/2 - ncol(mat1)/2,
  #                   offsety = offsety - max(offsety)/2 - nrow(mat1)/2) %>%
  #     dplyr::filter(ccfReweighted == max(ccfReweighted,na.rm = TRUE))
  #
  #   return(data.frame("ccf" = resp$ccfReweighted,"x" = resp$offsetx,"y" = resp$offsety))
  # }
  #
  # else if(ccfMethod == "imager"){
  #   resp <- imager::correlate(im = imager::as.cimg(mat2),
  #                             filter = imager::as.cimg(mat1),
  #                             normalise = TRUE) %>%
  #     as.matrix()
  # }

  # else if(ccfMethod == "fft"){
  resp <- filterViaFFT(mat1, mat2) / (sqrt(sum(mat1^2)) * sqrt(sum(mat2^2)))
  # }

  corr <- max(resp)
  tmp <- which(resp == corr, arr.ind = TRUE)[1, ]
  d_offset <- floor(dim(mat2)/2)

  dx <- tmp[["col"]] - d_offset[2] - 1
  dy <- tmp[["row"]] - d_offset[1] - 1

  ret <- data.frame("fft_ccf" = corr, "x" = dx, "y" = dy)
  return(ret)
}

# @name calcRawCorr
#
# @description Given the dy,dx values at which CCF_max occurs, this function
#   extract from the larger "region" matrix a matrix of the same dimension as
#   "cell" with appropriate shifting/scaling relative to the center based on
#   the dx,dy arguments. It is possible that the "cell" shaped matrix, after
#   shifting by dx,dy, lies outside of the bounds of the region matrix (this is
#   because both "cell" and "region" are padded with 0s when calculating the
#   CCF). If this occurs, the "cell"-sized matrix that would be extracted from
#   "region" is padded with NAs to appropriately reflect the conditions under
#   which the CCF_max value was determined. Even with this padding, however,
#   there may be a slight mis-match (like one row/col) between the dimensions
#   of "cell" and of the matrix extracted from "region," typically because
#   "region" has an odd dimension while "cell" has an even dimension or vice
#   versa. Thus, the "center" of one of the matrices might not be well-defined
#   as a single index. Thus, additional padding/cropping needs to be performed
#   in the "cell"-sized matrix extracted from "region." However, we don't know
#   whether these rows/cols should pre/post padded/cropped, so all possible
#   combinations of padding/cropping are considered. To determine which of
#   these these combinations should be ultimately chosen as the final
#   "cell"-sized matrix, the tieBreaker argument can be used to determine, for
#   example, which "cell"-sized matrix has the highest correlation with "cell"
#   (tieBreaker = which.max).
#
#
# @keywords internal
#
#' @importFrom stats cor

calcRawCorr <- function(cell,
                        region,
                        dx = 0,
                        dy = 0,
                        m1 = 0,
                        m2 = 0,
                        sd1 = 1,
                        sd2 = 1,
                        tieBreaker = which.max,
                        use = "pairwise.complete.obs"){
  cell <- (cell - m1)/sd1
  region <- (region - m2)/sd2

  #pad the region as it was when the CCF was calculated using the FFT method:
  regionPadded <- matrix(NA,
                         nrow = nrow(cell) + nrow(region) - 1,
                         ncol = ncol(cell) + ncol(region) - 1)
  regionPadded[1:nrow(region),1:ncol(region)] <- region

  regionCenter <- floor(dim(regionPadded)/2)

  alignedCenter <- regionCenter - c(dy + floor(nrow(cell)/2) + 1,dx + floor(ncol(cell)/2) + 1)

  alignedRows <- c((alignedCenter[1] - floor(dim(cell)[1]/2)),(alignedCenter[1] + floor(dim(cell)[1]/2)))

  alignedCols <- c((alignedCenter[2] - floor(dim(cell)[2]/2)),(alignedCenter[2] + floor(dim(cell)[2]/2)))

  regionCroppedInitial <- regionPadded[max(1,alignedRows[1]):min(nrow(region),alignedRows[2]),
                                       max(1,alignedCols[1]):min(ncol(region),alignedCols[2])]

  #build-in some contingencies in case the regionCropped isn't same size as
  #the cell. If the dimensions don't agree, then we need to consider copies
  #of regionCroppedInitial that have had rows/cols (whichever
  #applicable) cropped from the top/bottom or left/right
  regionCroppedList <- list("initial" = regionCroppedInitial)

  if(any(dim(regionCroppedInitial) != dim(cell))){
    #these make sure that the indices that we've cropped the region by aren't
    #less than 1 or larger than the dimension of the region
    if(alignedRows[1] < 1){
      rowsToPad <- -1*alignedRows[1] + 1

      regionCroppedInitial <- rbind(matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedInitial)),
                                    regionCroppedInitial)

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedRows[2] > nrow(region)){
      rowsToPad <- alignedRows[2] - nrow(region)

      regionCroppedInitial <- rbind(regionCroppedInitial,
                                    matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedInitial)))

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedCols[1] < 1){
      colsToPad <- -1*alignedCols[1] + 1

      regionCroppedInitial <- cbind(matrix(NA,nrow = nrow(regionCroppedInitial),ncol = colsToPad),
                                    regionCroppedInitial)

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedCols[2] > ncol(region)){
      colsToPad <- alignedCols[2] - ncol(region)

      regionCroppedInitial <- cbind(regionCroppedInitial,
                                    matrix(NA,nrow = nrow(regionCroppedInitial),ncol = colsToPad))

      regionCroppedList$initial <- regionCroppedInitial
    }

    #two copies if rows are off
    if(nrow(regionCroppedInitial) > nrow(cell) & ncol(regionCroppedInitial) == ncol(cell)){
      rowsToCrop <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedRowPre <- regionCroppedInitial[-rowsToCrop,]
      regionCroppedRowPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),]

      regionCroppedList$rowPre <- regionCroppedRowPre
      regionCroppedList$rowPost <- regionCroppedRowPost
    }

    #two copies if rows are off
    if(nrow(regionCroppedInitial) < nrow(cell) & ncol(regionCroppedInitial) == ncol(cell)){
      rowsToPad <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedRowPre <- rbind(matrix(NA,
                                          nrow = rowsToPad,
                                          ncol = ncol(regionCroppedInitial)),
                                   regionCroppedInitial)

      regionCroppedRowPost <- rbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedInitial)))

      regionCroppedList$rowPre <- regionCroppedRowPre
      regionCroppedList$rowPost <- regionCroppedRowPost
    }

    #2 copies if cols are off
    if(ncol(regionCroppedInitial) > ncol(cell) & nrow(regionCroppedInitial) == nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedInitial) - ncol(cell))

      regionCroppedColPre <- regionCroppedInitial[,-colsToCrop]
      regionCroppedColPost <- regionCroppedInitial[,-(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedList$colPre <- regionCroppedColPre
      regionCroppedList$colPost <- regionCroppedColPost
    }

    #2 copies if cols are off
    if(ncol(regionCroppedInitial) < ncol(cell) & nrow(regionCroppedInitial) == nrow(cell)){
      colsToPad <- abs(ncol(regionCroppedInitial) - ncol(cell))

      regionCroppedColPre <- cbind(matrix(NA,
                                          nrow = nrow(regionCroppedInitial),
                                          ncol = colsToPad),
                                   regionCroppedInitial)

      regionCroppedColPost <- cbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad))

      regionCroppedList$colPre <- regionCroppedColPre
      regionCroppedList$colPost <- regionCroppedColPost
    }

    #4 different copies if both dimensions are too large
    if(ncol(regionCroppedInitial) > ncol(cell) & nrow(regionCroppedInitial) > nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedInitial) - ncol(cell))

      rowsToCrop <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedBothPre <- regionCroppedInitial[-rowsToCrop,
                                                   -colsToCrop]

      regionCroppedRowPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),
                                                   -colsToCrop]

      regionCroppedColPost <- regionCroppedInitial[-rowsToCrop,
                                                   -(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedBothPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),
                                                    -(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedList$bothPre <- regionCroppedBothPre
      regionCroppedList$rowPost <- regionCroppedRowPost
      regionCroppedList$colPost <- regionCroppedColPost
      regionCroppedList$bothPost <- regionCroppedBothPost
    }

    #4 different copies if both dimensions are too small
    if(ncol(regionCroppedInitial) < ncol(cell) & nrow(regionCroppedInitial) < nrow(cell)){
      colsToPad <- abs(ncol(regionCroppedInitial) - ncol(cell))

      rowsToPad <- abs(nrow(regionCroppedInitial) - nrow(cell))

      #both rows and cols pre-padded
      regionCroppedBothPre <- cbind(matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad),
                                    regionCroppedInitial)
      regionCroppedBothPre <- rbind(matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedBothPre)),
                                    regionCroppedBothPre)

      #rows post-padded and cols pre-padded
      regionCroppedRowPost <- cbind(matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad),
                                    regionCroppedInitial)
      regionCroppedRowPost <- rbind(regionCroppedRowPost,
                                    matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedRowPost)))

      #rows pre-padded and cols post-padded
      regionCroppedColPost <- cbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad))
      regionCroppedColPost <- rbind(matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedColPost)),
                                    regionCroppedColPost)

      #rows and cols both post-padded
      regionCroppedBothPost <- cbind(regionCroppedInitial,
                                     matrix(NA,
                                            nrow = nrow(regionCroppedInitial),
                                            ncol = colsToPad))
      regionCroppedBothPost <- rbind(regionCroppedBothPost,
                                     matrix(NA,
                                            nrow = rowsToPad,
                                            ncol = ncol(regionCroppedBothPost)))

      regionCroppedList$bothPre <- regionCroppedBothPre
      regionCroppedList$rowPost <- regionCroppedRowPost
      regionCroppedList$colPost <- regionCroppedColPost
      regionCroppedList$bothPost <- regionCroppedBothPost
    }
  }

  #return NA if cor fails.
  corSafe <- purrr::safely(cor,otherwise = NA,quiet = TRUE)

  corrVals <- purrr::map_dbl(regionCroppedList,function(croppedRegion){

    #Very infrequently (like once per 100,000 cell/region comparisons) the
    #standard deviation is zero -- for example, in cases where a the two
    #matrices only overlap by one non-NA value. cor() will return NA in such
    #cases which we handle later on in the processing procedures. We'll suppress
    #the warnings.
    suppressWarnings(
      corVal <- corSafe(as.vector(cell),
                        as.vector(croppedRegion),
                        use = use)
    )
    as.numeric(corVal[1])
  })

  maxCorr <- as.numeric(corrVals[tieBreaker(corrVals)])
  if(length(maxCorr) == 0){
    return(NA)
  }
  else(return(maxCorr))
}
