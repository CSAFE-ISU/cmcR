#' Levels a breech face impression matrix basedo on a RANSAC-fitted plane
#'
#' @name preProcess_levelBF
#'
#' @description Given the output of preProcess_ransacLevel, extracts values (either
#'   raw or residual) from the surface matrix to which the RANSAC plane was fit.
#'   Adapted from the cartridges3D::levelBF3D function. This ia modified version
#'   of the levelBF3D function available in the cartridges3D package on GitHub.
#'
#' @param ransacFit output from the cmcR::preProcess_ransacLevel function.
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#'
#' @return a surface matrix of either "raw" breech face values that are inliers
#'   to the RANSAC-fitted plane or residuals between the fitted plane and
#'   observed values.
#'
#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#'   x3ptools::sample_x3p(m = 2)
#'
#' raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#'  cmcR::preProcess_ransacLevel() %>%
#'  cmcR::preProcess_levelBF(useResiduals = TRUE)
#' }
#'
#' @seealso https://github.com/xhtai/cartridges3D
#' @keywords internal
#'
#' @importFrom stats predict

preProcess_levelBF <- function(ransacFit,
                               useResiduals = TRUE){

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

#' Finds plane of breechface marks using the RANSAC method
#'
#' @description Given input depths (in microns), find best-fitting plane using
#'   RANSAC. This should be the plane that the breechface marks are on. Adapted
#'   from cartridges3D::findPlaneRansac function. This a modified version of the
#'   findPlaneRansac function available in the cartridges3D package on GitHub.
#'
#' @name preProcess_ransacLevel
#'
#' @param x3p an x3p object containing a surface matrix
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier". A smaller value will yield a more stable
#'   estimate.
#' @param ransacFinalSelectThresh once the RANSAC plane is fitted based on the
#'   ransacInlierThresh, this argument dictates which observations are selected as
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

#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#'   x3ptools::sample_x3p(m = 2)
#'
#' fittedPlane <- raw_x3p$surface.matrix %>%
#'   preProcess_ransacLevel(ransacInlierThresh = 10^(-5),
#'                     ransacFinalSelectThresh = 2*10^(-5),
#'                     iters = 150)
#' }
#'
#' @seealso
#'   https://github.com/xhtai/cartridges3D
#' @export
#'
#' @importFrom stats lm predict

preProcess_ransacLevel <- function(x3p,
                                   ransacInlierThresh = 1e-6,
                                   ransacFinalSelectThresh = 2e-5,
                                   iters = 300,
                                   returnResiduals = TRUE) {

  surfaceMat <- x3p$surface.matrix

  inlierCount <- 0

  # sample from this
  observedPixelLocations <- data.frame(which(!is.na(surfaceMat),
                                             arr.ind = TRUE)) %>%
    dplyr::mutate(depth = surfaceMat[!is.na(surfaceMat)])

  for (iter in 1:iters) {
    rowsToSample <- sample(nrow(observedPixelLocations),
                           3)

    candidatePlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[rowsToSample, ])

    #it's not important if a handful of the many iterations yields singular
    #matrices, only if the final approximation does -- this suppresses
    #intermediate warnings
    suppressWarnings(
      preds <- predict(candidatePlane, observedPixelLocations)
    )

    errors <- abs(preds - observedPixelLocations$depth)
    inlierBool <- errors < ransacInlierThresh

    if (sum(inlierBool) > inlierCount) { #if candidate plane is closer to more observed values, make this the new fitted plane
      finalPlaneErrors <- errors
      inlierCount <- sum(inlierBool)
      inliers <- inlierBool
    }
  }

  # final coefs only computed using inliers
  finalRansacPlane <- lm(depth ~ row + col,
                         data = observedPixelLocations[inliers, ]) #fit the plane based on what we've identified to be inliers

  #Once the plane is fitted based on the inliers identified, we want to take a potentially larger band of observations around the fitted plane than just the inlier threshold:
  finalInliers <- finalPlaneErrors < ransacFinalSelectThresh

  inlierLocations <- cbind(observedPixelLocations$row[finalInliers],
                           observedPixelLocations$col[finalInliers])

  #Now populate a new matrix to contain the estimated breech face
  estimatedBreechFace <- matrix(NA, nrow = nrow(surfaceMat), ncol = ncol(surfaceMat))

  estimatedBreechFace[inlierLocations] <- observedPixelLocations$depth[finalInliers]

  #Level the surface either by considering residuals or returning the surface matrix vertically-shifted to mean 0
  ransacFit <- list("ransacPlane" = finalRansacPlane,
                   "estimatedBreechFace" = estimatedBreechFace)

  estimatedBreechFace <- preProcess_levelBF(ransacFit = ransacFit,
                                            useResiduals = returnResiduals)

  x3p$surface.matrix <- estimatedBreechFace

  return(x3p)
}

#' Crop out rows/columns outside of the breech face impression in a cartridge
#' case scan.
#' @name preProcess_cropWS
#'
#' @param x3p an x3p object containing a surface matrix
#' @param croppingThresh minimum number of non-NA pixels that need to be in a
#'   row/column for it to not be cropped out of the surface matrix
#'
#' @return a surface matrix with outer rows/columns removed depending on
#'   croppingThresh
#'
#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#'   x3ptools::sample_x3p(m = 2)
#'
#' raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#'   cmcR::preProcess_ransacLevel() %>%
#'   cmcR::preProcess_levelBF() %>%
#'   cmcR::preProcess_cropWS(croppingThresh = 2)
#' }
#' @export

preProcess_cropWS <- function(x3p,
                              croppingThresh = 1){

  surfaceMat <- x3p$surface.matrix

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
  surfaceMatCropped <- surfaceMat[min(which(rowSum >= croppingThresh)):
                                    max(which(rowSum >= croppingThresh)),
                                  min(which(colSum >= croppingThresh)):
                                    max(which(colSum >= croppingThresh))]

  x3p$surface.matrix <- surfaceMatCropped

  #need to update metainformation now that rows/cols have been removed
  x3p$header.info$sizeX <- nrow(surfaceMatCropped)
  x3p$header.info$sizeY <- ncol(surfaceMatCropped)

  return(x3p)
}

#' Detect the radius and center of a firing pin impression circle in a breech
#' face impression scan using a circular Hough transform.
#'
#' @name preProcess_detectFPCircle
#'
#' @param x3p an x3p object containing a surface matrix
#' @param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#'   used to determine where the non-NA pixel count-per-row series attains a
#'   mode.
#' @param aggregationFunction function to aggregate all 12 radius estimates
#'   determined under the initial radius estimation procedure those calculated
#'   using fpRadiusGridSearch
#' @param gridGranularity granularity of radius grid used to determine the best
#'   fitting circle to the surface matrix via the Hough transform method
#' @param houghScoreQuant quantile cut-off to be used when determining a final
#'   radius estimate using the score values returned by the imager::hough_circle
#'   function
#'
#' @description This function estimates the radius of a firing pin impression
#'   within a breech face impression scan. It does so by detecting local maxima
#'   in the non-NA pixel count by row/col within the scan. To make the algorithm
#'   more robust, multiple radii estimates are considered by rotating the image
#'   15, 30, 45, 60, and 75 degrees and again detecting local maxima in the
#'   non-NA pixel count by row/col. Based on the argument passed to
#'   aggregationFunction, these radii estimates are reduced to a single, rough
#'   radius estimate (e.g., minimum was determined to be an effective
#'   aggregation function in preliminary tests). A grid of radii values centered
#'   on this estimate are then tested to determine which a final estimate. The
#'   grid mesh size is determined by the argument gridGranularity. A hough transform is
#'   applied to the breech face impression scan for each radius value in the
#'   grid. A final estimate is determined by finding the longest consecutive
#'   sequence of radii values with high associated hough scores. How we
#'   determine "high" hough scores is determined by the houghScoreQuant
#'   argument. Once the longest sequence of high hough score radii values is
#'   found, the average of these radii values is used as the final radius
#'   estimate.
#'
#' @keywords internal
#'
#' @importFrom stats quantile

utils::globalVariables(c(".","value","x","y","r"))

preProcess_detectFPCircle <- function(surfaceMat,
                                      aggregationFunction = mean,
                                      smootherSize = 2*round((.1*nrow(surfaceMat)/2)) + 1,
                                      gridSize = 40,
                                      gridGranularity = 1,
                                      houghScoreQuant = .9){

  firingPinRadiusEstimate <- fpRadiusGridSearch(surfaceMat = surfaceMat,smootherSize = smootherSize,
                                                aggregationFunction = aggregationFunction) %>%
    .$radiusEstim

  firingPinRadiusGrid <- seq(from = firingPinRadiusEstimate - floor(gridSize/2),
                             to = firingPinRadiusEstimate + floor(gridSize/2),
                             by = gridGranularity)

  surfaceMat_cannyEdges <- surfaceMat %>%
    is.na() %>%
    magrittr::not() %>%
    imager::as.cimg() %>%
    imager::cannyEdges()

  surfaceMat_houghCircleLocations <- purrr::map_dfr(firingPinRadiusGrid,
                                                    function(rad){
                                                      surfaceMat_houghCircle <-
                                                        imager::hough_circle(surfaceMat_cannyEdges,
                                                                             radius = rad) %>%
                                                        nms(50) %>%
                                                        as.data.frame() %>%
                                                        dplyr::top_n(1,wt = value) %>%
                                                        dplyr::group_by(value) %>%
                                                        dplyr::summarise(x=mean(x),y=mean(y),r = rad)
                                                    })

  q3_htValue <- surfaceMat_houghCircleLocations$value %>%
    quantile(houghScoreQuant)

  highValue_radii <- surfaceMat_houghCircleLocations$r[which(surfaceMat_houghCircleLocations$value >= q3_htValue)]
  breaks <- c(0,which(diff(highValue_radii) > gridGranularity),length(highValue_radii))

  consecutiveRadii <- sapply(seq(length(breaks) - 1),
                             function(i) highValue_radii[(breaks[i] + 1):breaks[i+1]])

  consecSeqLengths <- consecutiveRadii %>%
    purrr::map_int(length)

  finalRadiusEstimate <- consecutiveRadii[which(consecSeqLengths == max(consecSeqLengths))] %>%
    unlist() %>%
    mean()

  finalRadiusEstimate <- floor(finalRadiusEstimate)

  houghCircleLoc <- surfaceMat_houghCircleLocations %>%
    dplyr::filter(r == finalRadiusEstimate)

  return(houghCircleLoc)
}

#' Given a surface matrix, estimates and filters any pixels within the estimated
#' firing pin impression circle
#'
#' @name preProcess_removeFPCircle
#'
#' @param x3p an x3p object containing a surface matrix
#' @param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#' @param aggregationFunction function to select initial radius estimate from
#'   those calculated using fpRadiusGridSearch
#' @param gridSize size of grid, centered on the initial radius estimate, to be
#'   used to determine the best fitting circle to the surface matrix via the
#'   Hough transform method
#' @param gridGranularity granularity of radius grid used to determine the best
#'   fitting circle to the surface matrix via the Hough transform method
#' @param houghScoreQuant quantile cut-off to be used when determining a final
#'   radius estimate using the score values returned by the imager::hough_circle
#'
#' @note imager treats a matrix as its transpose (i.e., x and y axes are
#'   swapped). As such, relative to the original surface matrix, the x and y
#'   columns in the data frame fpImpressionCircle actually correspond to the row
#'   and column indices at which the center of the firing pin impression circle
#'   is estiamted to be.
#'
#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#' x3ptools::sample_x3p(m = 2)
#'
#' raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#'   cmcR::preProcess_ransacLevel() %>%
#'   cmcR::preProcess_levelBF() %>%
#'   cmcR::preProcess_cropWS() %>%
#'   cmcR::preProcess_removeFPCircle(aggregationFunction = mean,
#'                                   smootherSize = 2*round((.1*nrow(surfaceMat)/2)) + 1,
#'                                   gridGranularity = 1,
#'                                   houghScoreQuant = .9)
#' }
#'
#' @export

preProcess_removeFPCircle <- function(x3p,
                                      aggregationFunction = mean,
                                      smootherSize = 2*round((.1*nrow(surfaceMat)/2)) + 1,
                                      gridSize = 40,
                                      gridGranularity = 1,
                                      houghScoreQuant = .9){

  surfaceMat <- x3p$surface.matrix

  fpImpressionCircle <- preProcess_detectFPCircle(surfaceMat = surfaceMat,
                                                  aggregationFunction = aggregationFunction,
                                                  smootherSize = smootherSize,
                                                  gridSize = gridSize,
                                                  gridGranularity = gridGranularity,
                                                  houghScoreQuant = houghScoreQuant)

  breechFace_firingPinFiltered <- surfaceMat %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(value = ifelse(test = (x - fpImpressionCircle$x)^2 + (y - fpImpressionCircle$y)^2 >= (fpImpressionCircle$r)^2,
                                 yes = value,
                                 no = NA)) %>%
    imager::as.cimg(dim = c(max(.$x),max(.$y),1,1)) %>%
    as.matrix()

  x3p$surface.matrix <- breechFace_firingPinFiltered

  return(x3p)
}

#' Performs a low, high, or bandpass Gaussian filter on a surface matrix with a
#' particular cut-off wavelength.
#' @name preProcess_gaussFilter
#'
#' @param x3p an x3p object containing a surface matrix
#' @param wavelength cut-off wavelength
#' @param filtertype specifies whether a low pass, "lp", high pass, "hp", or bandpass,
#'   "bp" filter is to be used. Note that setting filterype = "bp" means that
#'   wavelength should be a vector of two numbers. In this case, the max of
#'   these two number will be used for the high pass filter and the min for the
#'   low pass filter.
#'
#' @examples
#' \dontrun{
#' raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#' x3ptools::sample_x3p(m = 2)
#'
#' raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#'   cmcR::preProcess_ransacLevel() %>%
#'   cmcR::preProcess_levelBF() %>%
#'   cmcR::preProcess_cropWS() %>%
#'   cmcR::preProcess_removeFPCircle() %>%
#'   cmcR::preProcess_gaussFilter(wavelength = c(16,500),
#'                                filtertype = "bp")
#' }
#'
#' @seealso
#' https://www.mathworks.com/matlabcentral/fileexchange/61003-filt2-2d-geospatial-data-filter?focused=7181587&tab=example
#'
#' @export

preProcess_gaussFilter <- function(x3p,
                                   wavelength = c(16,500),
                                   filtertype = "bp"){

  surfaceMat <- x3p$surface.matrix
  res <- x3p$header.info$incrementY

  if(res < .0001){ #rescale surface matrix for intermediate calculations (FFT seems to struggle with excessively small values)
    res <- res*1e6
  }

  surfaceMatMissing <- is.na(surfaceMat)

  surfaceMatFake <- surfaceMat - mean(as.vector(surfaceMat),na.rm=TRUE)
  surfaceMatFake[is.na(surfaceMatFake)] <- 0

  surfaceMatFake <- surfaceMatFake*1e6

  surfaceMatFiltered <- gaussianFilter(surfaceMat = surfaceMatFake,
                                       res = res,
                                       wavelength = wavelength,
                                       filtertype = filtertype)

  surfaceMatFiltered[surfaceMatMissing] <- NA

  surfaceMatFiltered <- surfaceMatFiltered/(1e6)

  x3p$surface.matrix <- surfaceMatFiltered

  return(x3p)
}
