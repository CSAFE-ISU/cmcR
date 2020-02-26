#' Finds plane of breechface marks using RANSAC
#'
#' Given input depths (in microns), find best-fitting plane using RANSAC. This
#' should be the plane that the breechface marks are on. Adapted from
#' cartridges3D::findPlaneRansac function
#'
#' @param surfaceMat matrix of input depths in microns.
#' @param inlierThreshold threshold to declare an observed value close to the
#'   fitted plane an "inlier".
#' @param iters number of candidate planes to fit (higher value yields more
#'   stable breech face estimate)
#'
#' @return List object containing fitted plane (as an lm object) and selected
#'   breechface marks (matrix of same size as original matrix containing all
#'   inliers close to fitted plane).
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#'
#' @seealso cartridges3D package (LINK)

findPlaneRansac <- function(surfaceMat,
                            inlierTreshold = 2*(10^(-5)),
                            iters = 150,...) {
  inlierCount <- 0

  # sample from this
  observedPixelLocations <- data.frame(which(!is.na(surfaceMat),
                                             arr.ind = TRUE)) %>%
    dplyr::mutate(depth = surfaceMat[!is.na(surfaceMat)])
  # observedPixelLocations$value <- surfaceMat[!is.na(surfaceMat)]

  for (iter in 1:iters) {
    rowsToSample <- sample(nrow(observedPixelLocations),
                           3)

    candidatePlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[rowsToSample, ])

    preds <- predict(candidatePlane, observedPixelLocations)

    errors <- abs(preds - observedPixelLocations$depth)
    inlierBool <- errors < inlierTreshold

    if (sum(inlierBool) > inlierCount) {
      inlierCount <- sum(inlierBool)
      inliers <- inlierBool
    }
  }

  # final coefs only computed using inliers
  finalRansacPlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[inliers, ])

  estimatedBreechFace <- matrix(NA, nrow = nrow(surfaceMat), ncol = ncol(surfaceMat))

  inlierLocations <- cbind(observedPixelLocations$row[inliers],
                           observedPixelLocations$col[inliers])

  estimatedBreechFace[inlierLocations] <- observedPixelLocations$depth[inliers]

  return(list(ransacPlane = finalRansacPlane,
              estimatedBreechFace = estimatedBreechFace))
}

#' Helper for selectBFImpression
#' @name levelBFImpression
#'
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value

levelBFImpression <- function(ransacFit,
                              useResiduals = FALSE,...){

  if(useResiduals){ #if the residuals from the RANSAC method are desired...
    fittedPlane <- ransacFit$estimatedBreechFace

    fittedPlane[!is.na(fittedPlane)] <- ransacFit$ransacPlane$fitted.values

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
#' @param surfaceMat
#' @param croppingThresh

cropScanWhitespace <- function(surfaceMat,
                               croppingThresh = 1,...){
  #Look at the middle 20% of columns and count the number of non-NA pixels in each
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
removeFPImpressionCircle <- function(bfImpression,fpImpressionCircle){
  breechFace_firingPinFiltered <- bfImpression %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(value = ifelse(test = (x - fpImpressionCircle$x)^2 + (y - fpImpressionCircle$y)^2 >= (fpImpressionCircle$r)^2,
                                 yes = value,
                                 no = NA)) %>%
    imager::as.cimg(dim = c(max(.$x),max(.$y),1,1)) %>%
    as.matrix()
}

#' @name selectBFImpression
#' @param x3p_path path to a .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizedBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param preProcess if FALSE, then no pre-processing is performed. Equivalent
#'   to calling x3ptools::read_x3p(x3p_path)
#' @return x3p object containing the processed breech face
#' @examples
#' \dontrun{
#'  processedBF1 <- cmcR::selectBFImpression(x3p_path = "path/to/file.x3p")
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)
#' @seealso x3ptools (LINK)

selectBFImpression <- function(x3p_path,
                               ransacInlierThresh = 2*(10^(-5)),
                               ransacIters = 150,
                               useResiduals = FALSE,
                               croppingThresh = 1,
                               standardizeBF = FALSE,
                               preProcess = TRUE,...){

  x3p <- x3ptools::read_x3p(x3p_path)

  if(!preProcess){
    return(x3p)
  }

  #First, we want to find the approximate height value(s) of the breech face impression within the cartridge case scan. We can find this using the RANSAC method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    iters = ransacIters) %>%
    levelBFImpression(useResiduals = useResiduals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThresh = croppingThresh) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan
  #even after the RANSAC method has selected the breech face impression height
  #values. We can remove these unwanted pixels by identifying the equation of
  #the firing pin impression circle and filtering out any pixels within that
  #circle. The following returns the estimated center and radius of the firing
  #pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

  x3p$surface.matrix <- bfImpressionFinal

  return(list("params" = list(ransacInlierThresh = ransacInlierThresh,
                              ransacIters = ransacIters,
                              useResiduals = useResiduals,
                              croppingThresh = croppingThresh,
                              standardizeBF = standardizeBF,
                              preProcess = preProcess),
              "x3p" = x3p))
}

#' Implements x3ptools::sample_x3p to downsample an x3p object's surface matrix.
#'
#' @name selectBFImpression_sample_x3p
#' @param x3p_path path to a .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizedBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param m integer value - every mth value is included in the sample
#' @param mY integer value - every mth value is included in the sample in the x
#'   direction and every mYth value is included in y direction
#' @param offset integer value between 0 and m-1 to specify offset of the sample
#' @param offsetY integer value between 0 and mY-1 to specify different offsets
#'   for x and y direction
#'
#' @note Given a matrix, x3ptools populates an x3p object's surface matrix
#'   starting in the top left corner moving right by reading the the matrix from
#'   the bottom left corner and moving up. Effectively rotating the matrix by 90
#'   degrees clockwise.
#' @return
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)
#' @seelalso x3ptools package (LINK)

selectBFImpression_sample_x3p <- function(x3p_path,
                                          ransacInlierThresh = 2*(10^(-5)),
                                          ransacIters = 150,
                                          useResiduals = FALSE,
                                          croppingThresh = 1,
                                          standardizeBF = FALSE,
                                          m = 2,
                                          mY = m,
                                          offset = 0,
                                          offsetY = offset){

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
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    iters = ransacIters) %>%
    levelBFImpression(useResiduals = useResiduals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThresh = croppingThresh) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan
  #even after the RANSAC method has selected the breech face impression height
  #values. We can remove these unwanted pixels by identifying the equation of
  #the firing pin impression circle and filtering out any pixels within that
  #circle. The following returns the estimated center and radius of the firing
  #pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

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

#' Implements imager::resize to resize an x3p object's surface matrix to an
#' arbitrary size.
#' @name selectBFImpression_resize
#' @param x3p_path path to a .x3p file
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier" for the RANSAC method
#' @param ransacIters number of candidate planes to fit for the RANSAC method
#'   (higher value yields more stable breech face estimate)
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#' @param croppingThresh treshold for cropping exterior whitespace from the
#'   breech face impression
#' @param standardizedBF subtract mean and divide by standard deviation of
#'   processed surface matrix
#' @param size_x,
#' @param size_y,
#' @param interpolation_type = 1,
#' @param boundary_conditions = 0
#'
#' @note imager treats a matrix as its transpose (i.e., x and y axes are
#'   swapped). As such the size_x argument corresponds to changing the number of
#'   rows in the original matrix and size_y the number of columns.
#' @return
#' @examples
#' \dontrun{
#'
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)
#' @seelalso imager package (LINK)

selectBFImpression_resize <- function(x3p_path,
                                      ransacInlierThresh = 2*(10^(-5)),
                                      ransacIters = 150,
                                      useResiduals = FALSE,
                                      croppingThresh = 1,
                                      standardizeBF = FALSE,
                                      size_x,
                                      size_y = size_x,
                                      interpolation_type = 1,
                                      boundary_conditions = 0){

  x3p <- x3p_path %>%
    x3ptools::read_x3p()

  x3p$surface.matrix <- x3p$surface.matrix %>%
    imager::as.cimg() %>%
    imager::resize(size_x = size_x,
                   size_y = size_y,
                   interpolation_type = interpolation_type,
                   boundary_conditions = boundary_conditions) %>%
    as.matrix()

  #First, we want to find the approximate height value(s) of the breech face impression within the cartridge case scan. We can find this using the RANSAC method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    iters = ransacIters) %>%
    levelBFImpression(useResiduals = useResiduals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThresh = croppingThresh) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan even after the RANSAC method has selected the breech face impression height values. We can remove these unwanted pixels by identifying the equation of the firing pin impression circle and filtering out any pixels within that circle. The following returns the estimated center and radius of the firing pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  if(standardizeBF){
    bfImpressionFinal <- (bfImpressionFinal - mean(bfImpressionFinal,na.rm = TRUE))/sd(bfImpressionFinal,na.rm = TRUE)
  }

  #Important note: x3ptools and imager read an image from a surface matrix starting from the bottom left corner and moving upwards. Thus, a matrix ends up being treated like its transpose within these packages. For example, the "Y" axis of an x3p object actually refers to the columns in the surface matrix
  x3p$header.info <- list(sizeY = ncol(bfImpressionFinal),
                          sizeX = nrow(bfImpressionFinal),
                          incrementY = (x3p$incrementY/x3p$sizeY)*size_y,
                          incrementX = (x3p$incrementX/x3p$sizeX)*size_x)

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
