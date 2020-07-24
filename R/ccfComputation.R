# These functions were copied from the cartridges3D R package currently (as of
# 6/13/2020) available on GitHub: https://github.com/xhtai/cartridges3D.

#These functions were copied because they are internal to both the cartridges3D
#and cmcR packages and referring to another package's internal functions is not
#SOP in R package development (devtools::check warnings when I had it
#differently)


#' Calculates cross-correlation between two matrices using FFTs
#'
#' @name filterViaFFT
#'
#'
#' @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
#' @keywords internal
#'
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

#' Copies behavior of round() function in MATLAB
#'
#' @name round2
#' @seealso http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
#' @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
#' @keywords internal

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

#' Shifts CCF matrix output so that Nyquist frequency is in middle of matrix
#'
#' @name fftshift
#' @seealso http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
#' @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
#' @keywords internal

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

#' Performs a circular shift based on periodic boundary conditions in a vector
#'
#' @name circshift
#' @seealso http://stackoverflow.com/questions/18791212/equivalent-to-numpy-roll-in-r
#' @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
#' @keywords internal
#'
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

#' calculates the "brute force" pairwise-complete correlation between matrix a
#' and matrix b translated by (offsetx,offsety)
#'
#' @name pairwise_ccf
#'
#' @keywords internal

pairwise_ccf <- function(a, b, offsetx = 0, offsety = 0) {
  if ("x3p" %in% class(a)) a <- a$surface.matrix
  if ("x3p" %in% class(b)) b <- b$surface.matrix

  aadj <- b_pad <- matrix(NA,
                          nrow = nrow(a) + nrow(b),
                          ncol = ncol(a) + ncol(b))

  hdim <- (dim(a) + dim(b))/4

  b_pad[(hdim[1] + 1):(hdim[1] + dim(b)[1]),
        (hdim[2] + 1):(hdim[2] + dim(b)[2])] <- b

  idx <- (offsetx + 1):(offsetx + ncol(a))
  idy <- (offsety + 1):(offsety + nrow(a))
  aadj[idy,idx] <- a

  if(all(is.na(aadj)) | all(is.na(b_pad))){
    return(NA)
  }
  else{
    suppressWarnings(
      suppressMessages(
        list("ccf" = cor(as.numeric(aadj), as.numeric(b_pad), use = "pairwise.complete.obs"),
             "nonMissing" = sum(!is.na(as.numeric(aadj)) & !is.na(b_pad)))

      )
    )
  }
}

#' Computes the location of the maximum CCF value in a CCF map between two
#' matrices
#'
#' @name ccfComparison
#'
#' @param ccfMethod implements 3 different methods to calculate the CCF -- all which yield differing CCF values
#' @seealso cartridges3D package \url{https://github.com/xhtai/cartridges3D}
#' @keywords internal

ccfComparison <- function(im1, im2, ccfMethod = "fftThenPairwise") {
  if(ccfMethod == "bruteForceReweighted"){
    resp <- expand.grid(offsetx = 1:(max(ncol(im1),ncol(im2))),
                        offsety = 1:(max(nrow(im1),nrow(im2)))) %>%
      purrr::pmap_dfr(~ {
        pwiseCCF <- pairwise_ccf(im1,im2,offsetx = ..1,offsety = ..2)

        data.frame(offsetx = ..1,
                   offsety = ..2,
                   ccf = pwiseCCF$ccf,
                   nonMissing = pwiseCCF$nonMissing)
      }) %>%
      dplyr::mutate(ccfReNorm = nonMissing*ccf/max(nonMissing),
                    offsetx = offsetx - max(offsetx)/2 - ncol(im1)/2,
                    offsety = offsety - max(offsety)/2 - nrow(im1)/2) %>%
      dplyr::filter(ccfReNorm == max(ccfReNorm,na.rm = TRUE))

    return(list("ccf" = resp$ccfReNorm,"dx" = resp$offsetx,"dy" = resp$offsety))
  }

  else if(ccfMethod == "imager"){
    resp <- imager::correlate(im = imager::as.cimg(im2),
                              filter = imager::as.cimg(im1),
                              normalise = TRUE) %>%
      as.matrix()
  }

  else if(ccfMethod == "fftThenPairwise"){
    resp <- filterViaFFT(im1, im2) / (sqrt(sum(im1^2)) * sqrt(sum(im2^2)))
  }

  corr <- max(resp)
  tmp <- which(resp == corr, arr.ind = TRUE)[1, ]
  d_offset <- floor(dim(im2)/2)

  dx <- tmp[["col"]] - d_offset[2] - 1
  dy <- tmp[["row"]] - d_offset[1] - 1

  ret <- list("ccf" = corr, "dx" = dx, "dy" = dy)
  return(ret)
}
