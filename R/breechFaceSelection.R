#' Finds plane of breechface marks using RANSAC
#'
#' Given input depths (in microns), find best-fitting plane using RANSAC. This
#' should be the plane that the breechface marks are on. This a modified version
#' of the findPlaneRansac function available in the cartridges3D package on
#' GitHub.
#'
#' @param surfaceMat matrix of input depths in microns.
#' @param inlierThreshold threshold to declare an observed value close to the
#'   fitted plane an "inlier". A smaller value will yield a more stable
#'   estimate.
#' @param finalSelectionThreshold once the RANSAC plane is fitted based on the
#'   inlierThreshold, this argument dictates which observations are selected as
#'   the final breech face estimate.
#' @param iters number of candidate planes to fit (higher value yields more
#'   stable breech face estimate)
#'
#' @return List object containing fitted plane (as an lm object) and selected
#'   breechface marks (matrix of same size as original matrix containing all
#'   inliers close to fitted plane).
#'
#' @note The function will throw an error if the final plane estimate is
#'   rank-deficient (which is relatively unlikely, but theoretically possible).
#'   Re-run the function (possibly setting a different seed) if this occurs.
#'
#' @examples
#' \dontrun{
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#'
#' @seealso https://github.com/xhtai/cartridges3D
#'
#' @keywords internal
#'
#' @importFrom stats lm predict

findPlaneRansac <- function(surfaceMat,
                            inlierTreshold = 1e-5, # 1 micron
                            finalSelectionThreshold = 2*1e-5, # 2 micron
                            iters = 150,...) {
  assertthat::not_empty(surfaceMat)
  testthat::expect_true(is.matrix(surfaceMat))
  assertthat::is.number(inlierTreshold)
  testthat::expect_gt(inlierTreshold,0)
  assertthat::is.number(finalSelectionThreshold)
  testthat::expect_gt(finalSelectionThreshold,0)
  assertthat::is.number(iters)
  testthat::expect_gt(iters,0)

  inlierCount <- 0

  # sample from this
  observedPixelLocations <- data.frame(which(!is.na(surfaceMat),
                                             arr.ind = TRUE)) %>%
    dplyr::mutate(depth = surfaceMat[!is.na(surfaceMat)])

  assertthat::not_empty(observedPixelLocations)

  for (iter in 1:iters) {
    rowsToSample <- sample(nrow(observedPixelLocations),3)

    candidatePlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[rowsToSample, ])

    suppressWarnings(
      preds <- predict(candidatePlane, observedPixelLocations)
    )

    errors <- abs(preds - observedPixelLocations$depth)
    inlierBool <- errors < inlierTreshold

    #if candidate plane is closer to more observed values, make this the new
    #fitted plane
    if (sum(inlierBool) > inlierCount) {
      finalPlaneErrors <- errors
      inlierCount <- sum(inlierBool)
      inliers <- inlierBool
    }
  }

  # final coefs only computed using inliers. fit the plane based on what we've
  # identified to be inliers. We will not allow this final estimate to be
  # rank-deficient
  finalRansacPlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[inliers, ],
                         singular.ok = FALSE)

  #Once the plane is fitted based on the inliers identified, we want to take a
  #potentially larger band of observations around the fitted plane than just the
  #inlier threshold:
  finalInliers <- finalPlaneErrors < finalSelectionThreshold

  inlierLocations <- cbind(observedPixelLocations$row[finalInliers],
                           observedPixelLocations$col[finalInliers])

  #Now populate a new matrix to contain the estimated breech face
  estimatedBreechFace <- matrix(NA, nrow = nrow(surfaceMat), ncol = ncol(surfaceMat))

  estimatedBreechFace[inlierLocations] <- observedPixelLocations$depth[finalInliers]

  testthat::expect_s3_class(finalRansacPlane,class = "lm")
  testthat::expect_true(is.matrix(estimatedBreechFace))

  return(list(ransacPlane = finalRansacPlane,
              estimatedBreechFace = estimatedBreechFace))
}

#' Helper for selectBFImpression. This ia modified version of the levelBF3D
#' function available in the cartridges3D package on GitHub.
#' @name levelBFImpression
#'
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#'
#' @seealso https://github.com/xhtai/cartridges3D
#'
#' @keywords internal
#'
#' @importFrom stats predict

levelBFImpression <- function(ransacFit,
                              useResiduals = TRUE,...){

  if(useResiduals){ #if the residuals from the RANSAC method are desired...
    esimatedBFdf <- data.frame(which(!is.na(ransacFit$estimatedBreechFace),
                                     arr.ind = TRUE)) %>%
      dplyr::mutate(depth = ransacFit$estimatedBreechFace[!is.na(ransacFit$estimatedBreechFace)])

    preds <- predict(ransacFit$ransacPlane,
                     newdata = esimatedBFdf)

    fittedPlane <- ransacFit$estimatedBreechFace
    fittedPlane[!is.na(fittedPlane)] <- preds

    # then take residuals
    resids <- ransacFit$estimatedBreechFace - fittedPlane

    return(resids)
  }
  else{
    #otherwise, just return estimated breech face centered vertically at 0
    bfEstim <- ransacFit$estimatedBreechFace -
      mean(as.vector(ransacFit$estimatedBreechFace),
           na.rm = TRUE)

    return(bfEstim)
  }
}

#' Crop out rows/columns outside of the breech face impression in a cartridge
#' case scan. Helper for selectBFImpression
#' @name cropScanWhitespace
#'
#' @param surfaceMat matrix of input depths in microns.
#' @param croppingThresh number of pixels that need to be observed in rows/cols
#'   before cropping from the exterior stops
#'
#' @keywords internal

cropScanWhitespace <- function(surfaceMat,
                               croppingThresh = 1,...){
  #Look at the middle 20% of columns and count the number of non-NA pixels in
  #each
  colSum <- surfaceMat[(nrow(surfaceMat)/2 - .1*nrow(surfaceMat)):
                         (nrow(surfaceMat)/2 + .1*nrow(surfaceMat)),] %>%
    is.na() %>%
    magrittr::not() %>%
    colSums()

  #Look at the middle 20% of rows and count the number of non-NA pixels in each
  rowSum <- surfaceMat[,(ncol(surfaceMat)/2 - .1*ncol(surfaceMat)):
                         (ncol(surfaceMat)/2 + .1*ncol(surfaceMat))] %>%
    is.na() %>%
    magrittr::not() %>%
    rowSums()

  #Crop out any rows/columns containing only NA pixels
  surfaceMatCropped <- surfaceMat[min(which(rowSum > croppingThresh)):
                                    max(which(rowSum > croppingThresh)),
                                  min(which(colSum > croppingThresh)):
                                    max(which(colSum > croppingThresh))]

  return(surfaceMatCropped)
}

#' @name removeFPImpressionCircle
#'
#' @keywords internal

utils::globalVariables(c("x","y","value","."))

removeFPImpressionCircle <- function(bfImpression,fpImpressionCircle){
  breechFace_firingPinFiltered <- bfImpression %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(value = ifelse(test = (x - fpImpressionCircle$x)^2 + (y - fpImpressionCircle$y)^2 >= (fpImpressionCircle$r)^2,
                                 yes = value,
                                 no = NA)) %>%
    imager::as.cimg(dim = c(max(.$x),max(.$y),1,1)) %>%
    as.matrix()

  return(breechFace_firingPinFiltered)
}

#' Select breech face impression from a cartridge case scan
#'
#' @name selectBFImpression
#'
#' @description Given a string representing a path to a cartridge case scan .x3p
#'   file, this function will read and process the scan including (1) use the
#'   RANSAC iterative plane fitting method to detect the height value at which
#'   the breech face impressions are in the scan, (2) crop out extra whitespace,
#'   NA values from the exterior of the resulting breech face impression surface
#'   matrix, and (3) detect and filter the inner firing pin impression circle
#'   using a circular Hough transform. Optionally, a Gaussian filter can be
#'   applied. See the cmcR::preProcess_gaussFilter function for more
#'   information.
#'
#'   Note that the RANSAC method is applied twice to the scan: once to get a
#'   rough estimate of the breech face impression height values (controlled with
#'   the ransacInlierThresh and ransacIters arguments) and again using
#'   .1*ransacInlierThresh and 2*ransacIters to get a more precise estimate.
#'
#' @param x3p_path path to an .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacFinalSelectThresh once the RANSAC plane is fitted based on the
#'   inlierThreshold, this argument dictates which observations are selected as
#'   the final breech face estimate.
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizeBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param gaussFilterRes sampling resolution of the surface matrix (given by
#'   incrementX or incrementY of the header.info element of an x3p object)
#' @param gaussFilterWavelength cut-off wavelength to be attenuated
#' @param gaussFilterType specifies whether a low pass, "lp", high pass, "hp",
#'   or bandpass, "bp" filter is to be used. Note that setting filterype = "bp"
#'   means that wavelength should be a vector of two numbers. In this case, the
#'   max of these two number will be used for the high pass filter and the min
#'   for the low pass filter.
#'
#' @return x3p object containing the processed breech face
#'
#' @examples
#' \dontrun{
#'  processedBF1 <- cmcR::selectBFImpression(x3p_path = "path/to/file.x3p")
#' }
#'
#' @seealso cartridges3D package (LINK)
#' @seealso x3ptools (LINK)
#'
#' @export
#'
#' @importFrom stats sd

utils::globalVariables("preProcess")

selectBFImpression <- function(x3p_path,
                               ransacInlierThresh = 1e-6,
                               ransacFinalSelectThresh = 2*1e-5,
                               ransacIters = 300,
                               useResiduals = TRUE,
                               croppingThresh = 1,
                               standardizeBF = FALSE,
                               gaussFilterRes = NULL,
                               gaussFilterWavelength = c(16,250),
                               gaussFilterType = "bp"){

  x3p <- x3ptools::read_x3p(x3p_path)

  #First, we want to find the approximate height value(s) of the breech face
  #impression within the cartridge case scan. We can find this using the RANSAC
  #method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = 10*ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ceiling(ransacIters/2)) %>%
    levelBFImpression(useResiduals = useResiduals) %>%
    #do it again, but with more restrictive thresholds:
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ransacIters) %>%
    #either returns residuals between fitted RANSAC plane and observed cartridge
    #case scan values or just returns the raw values of the estimated bf
    #impression
    levelBFImpression(useResiduals = useResiduals) %>%
    #also crop out whitespace on exterior of cartridge case scan
    cropScanWhitespace(croppingThresh = croppingThresh)

  stopifnot(is.matrix(bfImpression_ransacSelected),
            all(dim(bfImpression_ransacSelected) > c(0,0)),
            sum(!is.na(bfImpression_ransacSelected)) > 0)

  #Some additional, unwanted pixels remain in the middle of the cartridge scan
  #even after the RANSAC method has selected the breech face impression height
  #values. We can remove these unwanted pixels by identifying the equation of
  #the firing pin impression circle and filtering out any pixels within that
  #circle. The following returns the estimated center and radius of the firing
  #pin impression circle:
  fpImpressionCircle <-  preProcess_detectFPCircle(surfaceMat = bfImpression_ransacSelected,
                                                   aggregation_function = mean,
                                                   smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                                   meshSize = 1,
                                                   houghScoreQuant = .9)

  stopifnot(is.data.frame(fpImpressionCircle), nrow(fpImpressionCircle) == 1)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  stopifnot(is.matrix(bfImpressionFinal))

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

  if(length(gaussFilterWavelength) == 2 & is.null(gaussFilterType)){
    gaussFilterType <- "bp"
  }

  if(missing(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    gaussFilterRes <- x3p$header.info$incrementY
  }

  if(!is.null(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    bfImpressionFinal <- bfImpressionFinal %>%
      preProcess_gaussFilter(res = gaussFilterRes,
                             wavelength = gaussFilterWavelength,
                             filtertype = gaussFilterType)
  }

  stopifnot(is.matrix(bfImpressionFinal),
            all(dim(bfImpressionFinal) > c(0,0)),
            sum(!is.na(bfImpressionFinal)) > 0)

  x3p$surface.matrix <- bfImpressionFinal

  return(list("params" = list(ransacInlierThresh = ransacInlierThresh,
                              ransacIters = ransacIters,
                              useResiduals = useResiduals,
                              croppingThresh = croppingThresh,
                              standardizeBF = standardizeBF,
                              preProcess = preProcess),
              "x3p" = x3p))
}

#' Select breech face impression from a cartridge case scan after using
#' x3ptools::sample_x3p to downsample the scan.
#'
#' @name selectBFImpression_sample_x3p
#' @param x3p_path path to an .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacFinalSelectThresh once the RANSAC plane is fitted based on the
#'   inlierThreshold, this argument dictates which observations are selected as
#'   the final breech face estimate.
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizeBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param m integer value - every mth value is included in the sample
#' @param mY integer value - every mth value is included in the sample in the x
#'   direction and every mYth value is included in y direction
#' @param offset integer value between 0 and m-1 to specify offset of the sample
#' @param offsetY integer value between 0 and mY-1 to specify different offsets
#'   for x and y direction.
#' @param gaussFilterRes sampling resolution of the surface matrix
#'   (given by incrementX or incrementY of the header.info element of an x3p
#'   object). If not given, then this is determined automatically within the
#'   function
#' @param gaussFilterWavelength cut-off wavelength(s) to be attenuated
#' @param gaussFilterType specifies whether a low pass, "lp", high
#'   pass, "hp", or bandpass, "bp" filter is to be used. Note that setting
#'   filterype = "bp" means that wavelength should be a vector of two numbers.
#'   In this case, the max of these two number will be used for the high pass
#'   filter and the min for the low pass filter.
#'
#' @note Given a matrix, x3ptools populates an x3p object's surface matrix
#'   starting in the top left corner moving right by reading the the matrix from
#'   the bottom left corner and moving up. Effectively rotating the matrix by 90
#'   degrees clockwise.
#'
#' @return x3p object containing the processed, down-sampled breech face
#'
#' @examples
#' \dontrun{
#'  #downsample x3p to a quarter the original size (every other row/col selected)
#'  processedBF1 <- cmcR::selectBFImpression_sample_x3p(x3p_path = "path/to/file.x3p",m = 2)
#' }
#' @export
#'
#' @seealso cartridges3D package
#' @seealso x3ptools package
#'
#' @importFrom stats sd

utils::globalVariables(".")

selectBFImpression_sample_x3p <- function(x3p_path,
                                          ransacInlierThresh = 1e-6,
                                          ransacFinalSelectThresh = 2*1e-5,
                                          ransacIters = 300,
                                          useResiduals = TRUE,
                                          croppingThresh = 1,
                                          standardizeBF = FALSE,
                                          m = 2,
                                          mY = m,
                                          offset = 0,
                                          offsetY = offset,
                                          gaussFilterRes = NULL,
                                          gaussFilterWavelength = c(16,250),
                                          gaussFilterType = "bp"){

  x3p <- x3p_path %>%
    x3ptools::read_x3p() %>%
    x3ptools::sample_x3p(x3p = .,
                         m = m,
                         mY = mY,
                         offset = offset,
                         offsetY = offsetY)

  #First, we want to find the approximate height value(s) of the breech face
  #impression within the cartridge case scan. We can find this using the RANSAC
  #method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = 10*ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ceiling(ransacIters/2)) %>%
    #either returns residuals between fitted RANSAC plane and observed cartridge
    #case scan values or just returns the raw values of the estimated bf
    #impression
    levelBFImpression(useResiduals = useResiduals) %>%
    #do it again, but with more restrictive thresholds:
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ransacIters) %>%
    levelBFImpression(useResiduals = useResiduals) %>%
    #also crop out whitespace on exterior of cartridge case scan
    cropScanWhitespace(croppingThresh = croppingThresh)

  stopifnot(is.matrix(bfImpression_ransacSelected),
            all(dim(bfImpression_ransacSelected) > c(0,0)),
            sum(!is.na(bfImpression_ransacSelected)) > 0)

  #Some additional, unwanted pixels remain in the middle of the cartridge scan
  #even after the RANSAC method has selected the breech face impression height
  #values. We can remove these unwanted pixels by identifying the equation of
  #the firing pin impression circle and filtering out any pixels within that
  #circle. The following returns the estimated center and radius of the firing
  #pin impression circle:
  fpImpressionCircle <- preProcess_detectFPCircle(surfaceMat = bfImpression_ransacSelected,
                                                  aggregation_function = mean,
                                                  smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                                  meshSize = 1,
                                                  houghScoreQuant = .9)

  stopifnot(is.data.frame(fpImpressionCircle),
            nrow(fpImpressionCircle) == 1)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

  if(length(gaussFilterWavelength) == 2 & is.null(gaussFilterType)){
    gaussFilterType <- "bp"
  }

  if(missing(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    gaussFilterRes <- x3p$header.info$incrementY
  }

  if(!is.null(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    bfImpressionFinal <- bfImpressionFinal %>%
      preProcess_gaussFilter(res = gaussFilterRes,
                             wavelength = gaussFilterWavelength,
                             filtertype = gaussFilterType)
  }

  stopifnot(is.matrix(bfImpressionFinal),
            all(dim(bfImpressionFinal) > c(0,0)),
            sum(!is.na(bfImpressionFinal)) > 0)

  #Important note: x3ptools and imager read an image from a surface matrix
  #starting from the bottom left corner and moving upwards. Thus, a matrix ends
  #up being treated like its transpose within these packages. For example, the
  #"Y" axis of an x3p object actually refers to the columns in the surface
  #matrix
  x3p$header.info$sizeY <- ncol(bfImpressionFinal)
  x3p$header.info$sizeX <- nrow(bfImpressionFinal)

  x3p$surface.matrix <- bfImpressionFinal

  return(list("params" = list(ransacInlierThresh = ransacInlierThresh,
                              ransacIters = ransacIters,
                              useResiduals = useResiduals,
                              croppingThresh = croppingThresh,
                              standardizeBF = standardizeBF,
                              m = m,
                              mY = mY,
                              offset = offset,
                              offsetY = offsetY),
              "x3p" = x3p))
}

#' Select breech face impression from a cartridge case scan after using
#' imager::resize to resize the scan.
#'
#' @name selectBFImpression_resize
#'
#' @param x3p_path path to an .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacFinalSelectThresh threshold to declare an observed value close
#'   to the fitted plane an "inlier" for the RANSAC method
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizeBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param size_x Number of columns (new size along the X-axis). Note: the
#'   "X"-axis according to imager is actually the rows of the matrix
#' @param size_y Number of rows (new size along the Y-axis). Note: the "Y"-axis
#'   according to imager is actually the number of cols of the matrix
#' @param interpolation_type Method of interpolation: -1 = no interpolation: raw
#'   memory resizing. 0 = no interpolation: additional space is filled according
#'   to boundary_conditions. 1 = nearest-neighbor interpolation. 2 = moving
#'   average interpolation. 3 = linear interpolation. 4 = grid interpolation. 5
#'   = cubic interpolation. 6 = lanczos interpolation.
#' @param boundary_conditions Border condition type.
#' @param gaussFilterRes **Optional** sampling resolution of the surface matrix
#'   (given by incrementX or incrementY of the header.info element of an x3p
#'   object)
#' @param gaussFilterWavelength **Optional** cut-off wavelength to be attenuated
#' @param gaussFilterType **Optional** specifies whether a low pass, "lp", high
#'   pass, "hp", or bandpass, "bp" filter is to be used. Note that setting
#'   filterype = "bp" means that wavelength should be a vector of two numbers.
#'   In this case, the max of these two number will be used for the high pass
#'   filter and the min for the low pass filter.
#'
#' @note imager treats a matrix as its transpose (i.e., x and y axes are
#'   swapped). As such the size_x argument corresponds to changing the number of
#'   rows in the original matrix and size_y the number of columns.
#'
#' @return x3p object containing the processed, resized breech face
#'
#' @examples
#' \dontrun{
#'  #use moving average interpolation to resize matrix to be 560x560
#'  processedBF1 <- cmcR::selectBFImpression_resize(x3p_path="path/to/file.x3p",
#'                                                  size_x = 560,
#'                                                  size_y = 560,
#'                                                  interpolation_type = 2)
#' }
#' @export
#'
#' @seealso cartridges3D package
#' @seealso imager package
#'
#' @importFrom stats sd

selectBFImpression_resize <- function(x3p_path,
                                      ransacInlierThresh = 1e-6,
                                      ransacFinalSelectThresh = 2*1e-5,
                                      ransacIters = 300,
                                      useResiduals = TRUE,
                                      croppingThresh = 1,
                                      standardizeBF = FALSE,
                                      size_x,
                                      size_y = size_x,
                                      interpolation_type = 1,
                                      boundary_conditions = 0,
                                      gaussFilterRes = NULL,
                                      gaussFilterWavelength = c(16,250),
                                      gaussFilterType = "bp"){

  x3p <- x3p_path %>%
    x3ptools::read_x3p()

  x3p$surface.matrix <- x3p$surface.matrix %>%
    imager::as.cimg() %>%
    imager::resize(size_x = size_x,
                   size_y = size_y,
                   interpolation_type = interpolation_type,
                   boundary_conditions = boundary_conditions) %>%
    as.matrix()

  #First, we want to find the approximate height value(s) of the breech face
  #impression within the cartridge case scan. We can find this using the RANSAC
  #method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = 10*ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ceiling(ransacIters/2)) %>%
    levelBFImpression(useResiduals = useResiduals) %>%
    #do it again, but with more restrictive thresholds:
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    finalSelectionThreshold = ransacFinalSelectThresh,
                    iters = ransacIters) %>%
    #either returns residuals between fitted RANSAC plane and observed cartridge
    #case scan values or just returns the raw values of the estimated bf
    #impression
    levelBFImpression(useResiduals = useResiduals) %>%
    #also crop out whitespace on exterior of cartridge case scan
    cropScanWhitespace(croppingThresh = croppingThresh)

  stopifnot(is.matrix(bfImpression_ransacSelected),
            all(dim(bfImpression_ransacSelected) > c(0,0)),
            sum(!is.na(bfImpression_ransacSelected)) > 0)

  #Some additional, unwanted pixels remain in the middle of the cartridge scan
  #even after the RANSAC method has selected the breech face impression height
  #values. We can remove these unwanted pixels by identifying the equation of
  #the firing pin impression circle and filtering out any pixels within that
  #circle. The following returns the estimated center and radius of the firing
  #pin impression circle:
  fpImpressionCircle <- preProcess_detectFPCircle(surfaceMat = bfImpression_ransacSelected,
                                                  aggregation_function = mean,
                                                  smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                                  meshSize = 1,
                                                  houghScoreQuant = .9)

  stopifnot(is.data.frame(fpImpressionCircle),
            nrow(fpImpressionCircle) == 1)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

  if(length(gaussFilterWavelength) == 2 & is.null(gaussFilterType)){
    gaussFilterType <- "bp"
  }

  if(missing(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    gaussFilterRes <- (x3p$header.info$incrementY/nrow(x3p$surface.matrix))*nrow(bfImpressionFinal)
  }

  if(!is.null(gaussFilterRes) & !is.null(gaussFilterWavelength) & !is.null(gaussFilterType)){
    bfImpressionFinal <- bfImpressionFinal %>%
      preProcess_gaussFilter(res = gaussFilterRes,
                             wavelength = gaussFilterWavelength,
                             filtertype = gaussFilterType)
  }

  #Important note: x3ptools and imager read an image from a surface matrix
  #starting from the bottom left corner and moving upwards. Thus, a matrix ends
  #up being treated like its transpose within these packages. For example, the
  #"Y" axis of an x3p object actually refers to the columns in the surface
  #matrix
  x3p$header.info <- list(sizeY = ncol(bfImpressionFinal),
                          sizeX = nrow(bfImpressionFinal),
                          incrementY = (x3p$incrementY/x3p$sizeY)*size_y,
                          incrementX = (x3p$incrementX/x3p$sizeX)*size_x)

  stopifnot(is.matrix(bfImpressionFinal),
            all(dim(bfImpressionFinal) > c(0,0)),
            sum(!is.na(bfImpressionFinal)) > 0)

  x3p$surface.matrix <- bfImpressionFinal

  return(list("params" = list(ransacInlierThresh = ransacInlierThresh,
                              ransacIters = ransacIters,
                              useResiduals = useResiduals,
                              croppingThresh = croppingThresh,
                              standardizeBF = standardizeBF,
                              size_x = size_x,
                              size_y = size_y,
                              interpolation_type = interpolation_type,
                              boundary_conditions = boundary_conditions),
              "x3p" = x3p))
}
