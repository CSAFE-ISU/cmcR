#' Finds plane of breechface marks using RANSAC
#'
#' Given input depths (in microns), find best-fitting plane using RANSAC. This
#' should be the plane that the breechface marks are on. Adapted from cartridges3D::findPlaneRansac function
#'
#' @param surfaceMat matrix of input depths in microns.
#' @return List object containing fitted plane and selected breechface marks.
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)

findPlaneRansac <- function(surfaceMat,
                            inlierTreshold = 10^(-5),
                            iters = 75,...) {
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

#' @name levelBFImpression
levelBFImpression <- function(ransacFit,
                              use_residuals = TRUE,...){

  if(use_residuals){ #if the residuals from the RANSAC method are desired...
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
#' case scan
#' @name cropScanWhitespace
#'
#' @param surfaceMat
#' @param croppingThreshold

cropScanWhitespace <- function(surfaceMat,
                               croppingThreshold = 10,...){
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
  surfaceMatCropped <- surfaceMat[min(which(rowSum > croppingThreshold)):
                                    max(which(rowSum > croppingThreshold)),
                                  min(which(colSum > croppingThreshold)):
                                    max(which(colSum > croppingThreshold))]

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
#' @param
#' @return
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)

selectBFImpression <- function(x3p_path,
                               ransacInlierThresh = 10^(-5),
                               ransacIters = 75,
                               use_residuals = TRUE,
                               croppingThreshold = 10,...){

  x3p <- x3ptools::read_x3p(x3p_path)

  #First, we want to find the approximate height value(s) of the breech face impression within the cartridge case scan. We can find this using the RANSAC method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    iters = ransacIters) %>%
    levelBFImpression(use_residuals = use_residuals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThreshold = croppingThreshold) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan even after the RANSAC method has selected the breech face impression height values. We can remove these unwanted pixels by identifying the equation of the firing pin impression circle and filtering out any pixels within that circle. The following returns the estimated center and radius of the firing pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  x3p$surface.matrix <- bfImpressionFinal

  return(x3p)
}

#' @name selectBFImpression_sample_x3p
#' @param
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
                               ransacInlierThresh = 10^(-5),
                               ransacIters = 75,
                               use_residuals = TRUE,
                               croppingThreshold = 10,
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

  #First, we want to find the approximate height value(s) of the breech face impression within the cartridge case scan. We can find this using the RANSAC method
  bfImpression_ransacSelected <- x3p$surface.matrix %>%
    findPlaneRansac(inlierTreshold = ransacInlierThresh,
                    iters = ransacIters) %>%
    levelBFImpression(use_residuals = use_residuals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThreshold = croppingThreshold) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan even after the RANSAC method has selected the breech face impression height values. We can remove these unwanted pixels by identifying the equation of the firing pin impression circle and filtering out any pixels within that circle. The following returns the estimated center and radius of the firing pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  x3p$surface.matrix <- bfImpressionFinal

  return(x3p)
}

#' @name selectBFImpression_resize
#' @param
#' @return
#' @examples
#' \dontrun{
#'     set.seed(0)
#'     testImage <- findPlaneRansac(surfaceMat)
#' }
#' @export
#'
#' @seealso cartridges3D package (LINK)
#' @seelalso imager package (LINK)

selectBFImpression_resize <- function(x3p_path,
                                      size_x,
                                      size_y,
                                      interpolation_type = 1,
                                      boundary_conditions = 0,
                                      ransacInlierThresh = 10^(-5),
                                      ransacIters = 75,
                                      use_residuals = TRUE,
                                      croppingThreshold = 10){

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
    levelBFImpression(use_residuals = use_residuals) %>% #either returns residuals between fitted RANSAC plane and observed cartridge case scan values or just returns the raw values of the estimated bf impression
    cropScanWhitespace(croppingThreshold = croppingThreshold) #also crop out whitespace on exterior of cartridge case scan

  #Some additional, unwanted pixels remain in the middle of the cartridge scan even after the RANSAC method has selected the breech face impression height values. We can remove these unwanted pixels by identifying the equation of the firing pin impression circle and filtering out any pixels within that circle. The following returns the estimated center and radius of the firing pin impression circle:
  fpImpressionCircle <- fpCircleHoughDetection(bfImpression_ransacSelected,
                                               aggregation_function = mean,
                                               smootherSize = 2*round((.1*nrow(bfImpression_ransacSelected)/2)) + 1,
                                               meshSize = 1,
                                               houghScoreQuant = .9)

  #This then filters out any pixels within the firing pin
  bfImpressionFinal <- removeFPImpressionCircle(bfImpression = bfImpression_ransacSelected,
                                                fpImpressionCircle = fpImpressionCircle)

  x3p$surface.matrix <- bfImpressionFinal

  return(x3p)
}