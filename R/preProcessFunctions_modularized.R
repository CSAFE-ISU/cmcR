#' Finds plane of breechface marks using RANSAC
#'
#' Given input depths (in microns), find best-fitting plane using RANSAC. This
#' should be the plane that the breechface marks are on. Adapted from
#' cartridges3D::findPlaneRansac function
#'
#' @name preProcess_ransac
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
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
#' @examples
#' \dontrun{
#'     testImage <- preProcess_ransac(surfaceMat)
#' }
#'
#' @seealso cartridges3D package (LINK)
#' @export

preProcess_ransac <- function(surfaceMat,
                                   inlierTreshold = (10^(-5)), # 1 micron
                                   finalSelectionThreshold = 2*(10^(-5)), # 2 micron
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
  finalInliers <- finalPlaneErrors < finalSelectionThreshold

  inlierLocations <- cbind(observedPixelLocations$row[finalInliers],
                           observedPixelLocations$col[finalInliers])

  #Now populate a new matrix to contain the estimated breech face
  estimatedBreechFace <- matrix(NA, nrow = nrow(surfaceMat), ncol = ncol(surfaceMat))

  estimatedBreechFace[inlierLocations] <- observedPixelLocations$depth[finalInliers]

  return(list(ransacPlane = finalRansacPlane,
              estimatedBreechFace = estimatedBreechFace))
}

#' Given the output of preProcess_ransac, extracts values (either raw or residual) from the surface matrix to which the RANSAC plane was fit.
#'
#' @name preProcess_levelBF
#'
#' @param useResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#'
#' @export

preProcess_levelBF <- function(ransacFit,
                               useResiduals = FALSE,...){

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
#' case scan.
#' @name preProcess_cropWS
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param croppingThresh minimum number of non-NA pixels that need to be in a
#'   row/column for it to not be cropped out of the surface matrix
#'
#' @export

preProcess_cropWS <- function(surfaceMat,
                              croppingThresh = 2,...){
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

  return(surfaceMatCropped)
}

#' Detect the radius and center of a firing pin impression circle in a breech face
#' impression scan using a circular Hough transform.
#'
#' @name preProcess_detectFPCircle
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#' @param aggregation_function function to select initial radius estimate from
#'   those calculated using fpRadiusGridSearch
#' @param meshSize size of radius mesh to be used for further refinement of the
#'   radius estimate obtained from fpRadiusGridSearch
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
#'   aggregation_function, these radii estimates are reduced to a single, rough
#'   radius estimate (e.g., minimum was determined to be an effective
#'   aggregation function in preliminary tests). A grid of radii values centered
#'   on this estimate are then tested to determine which a final estimate. The
#'   grid mesh size is determined by the argument meshSize. A hough transform is
#'   applied to the breech face impression scan for each radius value in the
#'   grid. A final estimate is determined by finding the longest consecutive
#'   sequence of radii values with high associated hough scores. How we
#'   determine "high" hough scores is determined by the houghScoreQuant
#'   argument. Once the longest sequence of high hough score radii values is
#'   found, the average of these radii values is used as the final radius
#'   estimate.
#' @export

preProcess_detectFPCircle <- function(surfaceMat,
                                      aggregation_function = mean,
                                      smootherSize = 2*round((.1*nrow(surfaceMat)/2)) + 1,
                                      meshSize = 1,
                                      houghScoreQuant = .9){

  firingPinRadiusEstimate <- fpRadiusGridSearch(surfaceMat,smootherSize = smootherSize,
                                                aggregation_function = aggregation_function) %>%
    .$radiusEstim

  firingPinRadiusGrid <- seq(from = firingPinRadiusEstimate - 20,
                             to = firingPinRadiusEstimate + 20,
                             by = meshSize)

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
  breaks <- c(0,which(diff(highValue_radii) > meshSize),length(highValue_radii))

  consecutiveRadii <- sapply(seq(length(breaks) - 1),
                             function(i) highValue_radii[(breaks[i] + 1):breaks[i+1]])

  consecSeqLengths <- consecutiveRadii %>%
    purrr::map_int(length)

  finalRadiusEstimate <- consecutiveRadii[which(consecSeqLengths == max(consecSeqLengths))] %>%
    unlist() %>%
    mean() %>%
    floor()

  houghCircleLoc <- surfaceMat_houghCircleLocations %>%
    dplyr::filter(r == finalRadiusEstimate)

  return(houghCircleLoc)
}

#' Given a surface matrix and the output of preProcess_detectFPCircle, filters
#' any pixels within the estimated firing pin impression circle
#'
#' @name preProcess_removeFPCircle
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param fpImpressionCircle data frame containing 3 columns: "x" and "y"
#'   containing the estimated center of the firing pin impression circle and "r"
#'   containing the estimated radius
#'
#' @note imager treats a matrix as its transpose (i.e., x and y axes are
#'   swapped). As such, relative to the original surface matrix, the x and y
#'   columns in the data frame fpImpressionCircle actually correspond to the row
#'   and column indices at which the center of the firing pin impression circle
#'   is estiamted to be.
#'
#' @export

preProcess_removeFPCircle <- function(surfaceMat,fpImpressionCircle){
  breechFace_firingPinFiltered <- surfaceMat %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(value = ifelse(test = (x - fpImpressionCircle$x)^2 + (y - fpImpressionCircle$y)^2 >= (fpImpressionCircle$r)^2,
                                 yes = value,
                                 no = NA)) %>%
    imager::as.cimg(dim = c(max(.$x),max(.$y),1,1)) %>%
    as.matrix()

  return(breechFace_firingPinFiltered)
}

#' Performs a low, high, or bandpass Gaussian filter on a surface matrix with a
#' particular cut-off wavelength.
#' @name preProcess_gaussFilter
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param res sampling resolution of the surface matrix
#' @param wavelength cut-off wavelength
#' @filtertype specifies whether a low pass, "lp", high pass, "hp", or bandpass,
#'   "bp" filter is to be used. Note that setting filterype = "bp" means that
#'   wavelength should be a vector of two numbers. In this case, the max of
#'   these two number will be used for the high pass filter and the min for the
#'   low pass filter.
#'
#' @seealso https://www.mathworks.com/matlabcentral/fileexchange/61003-filt2-2d-geospatial-data-filter?focused=7181587&tab=example
#' @export

preProcess_gaussFilter <- function(surfaceMat,
                                   res = selectedBF_x3p$header.info$incrementY,
                                   wavelength,
                                   filtertype = "bp"){

  if(res < .00001){ #if resolution measured in meters:
    res <- res*(10^(6)) #rescale to microns
  }

  surfaceMatMissing <- is.na(selectedBF_x3p$surface.matrix)

  surfaceMatFake <- selectedBF_x3p$surface.matrix - mean(as.vector(selectedBF_x3p$surface.matrix),na.rm=TRUE)
  surfaceMatFake[is.na(surfaceMatFake)] <- 0

  surfaceMatFake <- surfaceMatFake*(10^6) #scale to microns (avoids small number numerical issues?)

  surfaceMatFiltered <- cmcR:::gaussianFilter(surfaceMat = surfaceMatFake,
                                              res = res,
                                              wavelength = wavelength,
                                              filtertype = filtertype)

  surfaceMatFiltered[surfaceMatMissing] <- NA

  surfaceMatFiltered <- surfaceMatFiltered/(10^6)

  selectBFImpression_output$x3p$surface.matrix <- surfaceMatFiltered
  return(selectBFImpression_output)
}
