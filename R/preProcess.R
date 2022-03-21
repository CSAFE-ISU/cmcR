# Levels a breech face impression matrix basedo on a RANSAC-fitted plane
#
# @name preProcess_levelBF
#
# @description Given the output of preProcess_ransacLevel, extracts values
#   (either raw or residual) from the surface matrix to which the RANSAC plane
#   was fit. Adapted from the cartridges3D::levelBF3D function. This ia
#   modified version of the levelBF3D function available in the cartridges3D
#   package on GitHub.
#
# @param ransacFit output from the cmcR::preProcess_ransacLevel function.
# @param useResiduals dictates whether the difference between the estimated
#   breech face and fitted plane are returned (residuals) or if the estimates
#   breech face is simply shifted down by its mean value
#
# @return a surface matrix of either "raw" breech face values that are inliers
#   to the RANSAC-fitted plane or residuals between the fitted plane and
#   observed values.
#
# @examples
# \dontrun{
# raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#   x3ptools::sample_x3p(m = 2)
#
# raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#  cmcR::preProcess_ransacLevel() %>%
#  cmcR::preProcess_levelBF(useResiduals = TRUE)
# }
#
# @seealso https://github.com/xhtai/cartridges3D
# @keywords internal
#
# @importFrom stats predict

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
#' @note Given input depths (in microns), find best-fitting plane using
#'   RANSAC. This should be the plane that the breechface marks are on. Adapted
#'   from cartridges3D::findPlaneRansac function. This a modified version of the
#'   findPlaneRansac function available in the cartridges3D package on GitHub.
#'
#' @param x3p an x3p object containing a surface matrix
#' @param ransacInlierThresh threshold to declare an observed value close to the
#'   fitted plane an "inlier". A smaller value will yield a more stable
#'   estimate.
#' @param ransacFinalSelectThresh once the RANSAC plane is fitted based on the
#'   ransacInlierThresh, this argument dictates which observations are selected
#'   as the final breech face estimate.
#' @param iters number of candidate planes to fit (higher value yields more
#'   stable breech face estimate)
#' @param returnResiduals dictates whether the difference between the estimated
#'   breech face and fitted plane are returned (residuals) or if the estimates
#'   breech face is simply shifted down by its mean value
#'@return an x3p object containing the leveled surface matrix.
#' @note The preProcess_ransacLevel function will throw an error if the final
#'   plane estimate is rank-deficient (which is relatively unlikely, but
#'   theoretically possible). Re-run the function (possibly setting a different
#'   seed) if this occurs.

#' @examples
#' \dontrun{
#' nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#' fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
#'
#' fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))
#'
#' fadul1.1_ransacLeveled <- fadul1.1 %>%
#'                      preProcess_crop(region = "exterior",
#'                                      radiusOffset = -30) %>%
#'                      preProcess_crop(region = "interior",
#'                                      radiusOffset = 200) %>%
#'                      preProcess_removeTrend(statistic = "quantile",
#'                                             tau = .5,
#'                                             method = "fn")
#'
#' x3pListPlot(list("Original" = fadul1.1,
#'                  "RANSAC Leveled" = fadul1.1_ransacLeveled),type = "list")
#'}
#'
#' @seealso
#'   https://github.com/xhtai/cartridges3D
#' @export
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

# Crop out rows/columns outside of the breech face impression in a cartridge
# case scan.
# @name preProcess_cropWS
#
# @param x3p an x3p object containing a surface matrix
# @param croppingThresh minimum number of non-NA pixels that need to be in a
#   row/column for it to not be cropped out of the surface matrix
#
# @return a surface matrix with outer rows/columns removed depending on
#   croppingThresh
#
# @examples
# \dontrun{
# raw_x3p <- x3ptools::read_x3p("path/to/file.x3p") %>%
#   x3ptools::sample_x3p(m = 2)
#
# raw_x3p$surface.matrix <- raw_x3p$surface.matrix %>%
#   cmcR::preProcess_ransacLevel() %>%
#   cmcR::preProcess_levelBF() %>%
#   cmcR::preProcess_cropWS(croppingThresh = 2)
# }
#
# @keywords internal

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

  minCroppingValues <- c("row" = min(which(rowSum >= croppingThresh)),
                         "col" = min(which(colSum >= croppingThresh)))

  x3p$cmcR.info$exteriorCroppingValues <- minCroppingValues

  return(x3p)
}

# Detect the radius and center of a firing pin impression circle in a breech
# face impression scan using a circular Hough transform.
#
# @name preProcess_detectFPCircle
#
# @param x3p an x3p object containing a surface matrix
# @param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#   used to determine where the non-NA pixel count-per-row series attains a
#   mode.
# @param aggregationFunction function to aggregate all 12 radius estimates
#   determined under the initial radius estimation procedure those calculated
#   using fpRadiusGridSearch
# @param gridGranularity granularity of radius grid used to determine the best
#   fitting circle to the surface matrix via the Hough transform method
# @param houghScoreQuant quantile cut-off to be used when determining a final
#   radius estimate using the score values returned by the imager::hough_circle
#   function
#
# @description This function estimates the radius of a firing pin impression
#   within a breech face impression scan. It does so by detecting local maxima
#   in the non-NA pixel count by row/col within the scan. To make the algorithm
#   more robust, multiple radii estimates are considered by rotating the image
#   15, 30, 45, 60, and 75 degrees and again detecting local maxima in the
#   non-NA pixel count by row/col. Based on the argument passed to
#   aggregationFunction, these radii estimates are reduced to a single, rough
#   radius estimate (e.g., minimum was determined to be an effective
#   aggregation function in preliminary tests). A grid of radii values centered
#   on this estimate are then tested to determine which a final estimate. The
#   grid mesh size is determined by the argument gridGranularity. A hough
#   transform is applied to the breech face impression scan for each radius
#   value in the grid. A final estimate is determined by finding the longest
#   consecutive sequence of radii values with high associated hough scores. How
#   we determine "high" hough scores is determined by the houghScoreQuant
#   argument. Once the longest sequence of high hough score radii values is
#   found, the average of these radii values is used as the final radius
#   estimate.
#
# @keywords internal
# @importFrom stats quantile
# @importFrom rlang .data

preProcess_detectFPCircle <- function(surfaceMat,
                                      aggregationFunction = mean,
                                      smootherSize = 2*round((.1*nrow(surfaceMat)/2)) + 1,
                                      gridSize = 40,
                                      gridGranularity = 1,
                                      houghScoreQuant = .9){

  firingPinRadiusEstimate <- fpRadiusGridSearch(surfaceMat = surfaceMat,smootherSize = smootherSize,
                                                aggregationFunction = aggregationFunction) %>%
    {.$radiusEstim}

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
                                                        dplyr::top_n(1,wt = .data$value) %>%
                                                        dplyr::group_by(.data$value) %>%
                                                        dplyr::summarise(x=mean(.data$x),y=mean(.data$y),r = rad)

                                                      return(surfaceMat_houghCircle)
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
    dplyr::filter(.data$r == finalRadiusEstimate)

  return(houghCircleLoc)
}

#'Given a surface matrix, estimates and filters any pixels within the estimated
#'firing pin impression circle
#'
#'@name preProcess_removeFPCircle
#'
#'@param x3p an x3p object containing a surface matrix
#'@param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#'@param aggregationFunction function to select initial radius estimate from
#'  those calculated using fpRadiusGridSearch
#'@param gridSize size of grid, centered on the initial radius estimate, to be
#'  used to determine the best fitting circle to the surface matrix via the
#'  Hough transform method
#'@param gridGranularity granularity of radius grid used to determine the best
#'  fitting circle to the surface matrix via the Hough transform method
#'@param houghScoreQuant quantile cut-off to be used when determining a final
#'  radius estimate using the score values returned by the imager::hough_circle
#'
#'@note imager treats a matrix as its transpose (i.e., x and y axes are
#'  swapped). As such, relative to the original surface matrix, the x and y
#'  columns in the data frame fpImpressionCircle actually correspond to the row
#'  and column indices at which the center of the firing pin impression circle
#'  is estiamted to be.
#'@return An x3p object containing a surface matrix with the estimated firing
#'  pin circle pixels replaced with NAs.
#' @examples
#' \dontrun{
#' nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#' fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
#'
#' fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))
#'
#' fadul1.1_labelCropped <- fadul1.1 %>%
#'                      preProcess_crop(region = "exterior",
#'                                      radiusOffset = -30) %>%
#'                      preProcess_crop(region = "interior",
#'                                      radiusOffset = 200) %>%
#'                      preProcess_removeTrend(statistic = "quantile",
#'                                             tau = .5,
#'                                             method = "fn")
#'
#' fadul1.1_houghCropped <- fadul1.1 %>%
#'                           x3ptools::x3p_sample() %>%
#'                           preProcess_ransacLevel() %>%
#'                           preProcess_crop(region = "exterior",
#'                                           radiusOffset = -30) %>%
#'                           preProcess_removeFPCircle()
#'
#' x3pListPlot(list("Original" = fadul1.1,
#'                  "Cropped by Labeling" = fadul1.1_labelCropped,
#'                  "Cropped by Hough" = fadul1.1_houghCropped),type = "list")
#'}
#'@importFrom rlang .data
#'@export

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
    dplyr::mutate(value = ifelse(test = (.data$x - fpImpressionCircle$x)^2 + (.data$y - fpImpressionCircle$y)^2 >= (fpImpressionCircle$r)^2,
                                 yes = .data$value,
                                 no = NA)) %>%
    imager::as.cimg(dim = c(max(.$x),max(.$y),1,1)) %>%
    as.matrix()

  x3p$surface.matrix <- breechFace_firingPinFiltered

  return(x3p)
}

#'Performs a low, high, or bandpass Gaussian filter on a surface matrix with a
#'particular cut-off wavelength.
#'@name preProcess_gaussFilter
#'
#'@param x3p an x3p object containing a surface matrix
#'@param wavelength cut-off wavelength
#'@param filtertype specifies whether a low pass, "lp", high pass, "hp", or
#'  bandpass, "bp" filter is to be used. Note that setting filterype = "bp"
#'  means that wavelength should be a vector of two numbers. In this case, the
#'  max of these two number will be used for the high pass filter and the min
#'  for the low pass filter.
#'@return An x3p object containing the Gaussian-filtered surface matrix.
#' @examples
#' data(fadul1.1_processed)
#'
#' #Applying the function to fadul1.1_processed (note that this scan has already
#' #  been Gaussian filtered)
#' cmcR::preProcess_gaussFilter(fadul1.1_processed)
#'
#' #As a part of the recommended preprocessing pipeline (take > 5 sec to run):
#' \dontrun{
#' nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#' fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
#'
#' fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))
#' fadul1.1_extCropped <- preProcess_crop(x3p = fadul1.1,
#'                                        region = "exterior",
#'                                        radiusOffset = -30)
#'
#' fadul1.1_intCroped <- preProcess_crop(x3p = fadul1.1_extCropped,
#'                                       region = "interior",
#'                                       radiusOffset = 200)
#'
#' fadul1.1_leveled <- preProcess_removeTrend(x3p = fadul1.1_intCroped,
#'                                            statistic = "quantile",
#'                                            tau = .5,
#'                                            method = "fn")
#' fadul1.1_filtered <- preProcess_gaussFilter(x3p = fadul1.1_leveled,
#'                                             wavelength = c(16,500),
#'                                             filtertype = "bp")
#'
#' x3pListPlot(list("Original" = fadul1.1,
#'                  "Ext. & Int. Cropped" = fadul1.1_intCroped,
#'                  "Cropped and Leveled" = fadul1.1_leveled,
#'                  "Filtered" = fadul1.1_filtered),type = "list")
#'}
#'
#'@seealso
#'https://www.mathworks.com/matlabcentral/fileexchange/61003-filt2-2d-geospatial-data-filter?focused=7181587&tab=example
#'@export

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

# @name estimateBFRadius
#
# @keywords internal
# @importFrom rlang .data

estimateBFRadius <- function(mat,
                             scheme = 3,
                             high_connectivity = FALSE,
                             tolerance = 0,
                             angle = 0,
                             interpolation = 0,
                             boundary = 0,
                             roughEstimate = FALSE,
                             agg_function = median){

  matFake <- (mat*1e5) + 1 #scale and shift all non-NA pixels up 1 (meter)
  matFakeRotated <- matFake %>%
    imager::as.cimg() %>%
    imager::imrotate(angle = angle,
                     interpolation = interpolation,
                     cx = floor(nrow(.)/2), #imager treats the rows as the "x" axis of an image
                     cy = floor(ncol(.)/2),
                     boundary = boundary) %>% #pad boundary with 0s (dirichlet condition)
    as.matrix()

  matFakeRotated[matFakeRotated == 0] <- NA
  #shift all of the legitimate pixels back down by 1:
  mat <- (matFakeRotated - 1)/(1e5)

  mat_binarized <- mat
  mat_binarized[!is.na(mat_binarized)] <- 1
  mat_binarized[is.na(mat_binarized)] <- 0

  mat_labeled <- mat_binarized %>%
    imager::as.cimg() %>%
    imager::imgradient(axes = "xy",
                       scheme = scheme) %>%
    imager::enorm() %>%
    imager::add() %>%
    imager::label(high_connectivity = high_connectivity,
                  tolerance = tolerance) %>%
    as.matrix()

  exteriorLabel <- data.frame(firstRow = mat_labeled[1,]) %>%
    dplyr::group_by(.data$firstRow) %>%
    dplyr::tally() %>%
    dplyr::top_n(n = 1,wt= .data$n) %>%
    dplyr::pull(.data$firstRow)

  mat_segmented <- mat_labeled
  mat_segmented[mat_segmented != exteriorLabel] <- -1
  mat_segmented[mat_segmented == exteriorLabel] <- -2

  mat_radiusInitialEstim <- round(sqrt(sum(mat_segmented == -1)/pi))

  mat_bfRegionCenter <- round(colMeans(which(mat_segmented == -1,arr.ind = TRUE)))

  mat_radiusEstimate <- mat_radiusInitialEstim

  if(!roughEstimate){
    mat_segmentedEdges <- mat_segmented %>%
      imager::as.cimg() %>%
      imager::imgradient(axes = "xy",
                         scheme = scheme) %>%
      imager::enorm() %>%
      imager::add() %>%
      as.matrix()

    mat_segmentedEdgesMidRow <- mat_segmentedEdges[mat_bfRegionCenter[1],]

    if(sum(mat_segmentedEdgesMidRow > 0) < 2){
      return(list("centerEstimate" = mat_bfRegionCenter,
                  "radiusEstimate" = mat_radiusEstimate))
    }

    mat_radiusEstimateCandidate <- data.frame(y = mat_segmentedEdgesMidRow) %>%
      dplyr::mutate(x = 1:nrow(.)) %>%
      dplyr::filter(.data$y != 0) %>%
      dplyr::mutate(x_lag = c(.data$x[2:(nrow(.))],NA)) %>%
      dplyr::mutate(x_diff = abs(.data$x - .data$x_lag)) %>%
      dplyr::top_n(n = 1,wt = .data$x_diff) %>%
      dplyr::summarise(x_lag = agg_function(.data$x_lag),
                       x = agg_function(.data$x)) %>%
      dplyr::summarise(radEstimate = round((.data$x_lag - .data$x)/2)) %>%
      dplyr::pull(.data$radEstimate) %>%
      agg_function(na.rm = TRUE)

    if(mat_radiusEstimateCandidate >= min(dim(mat)/4) & mat_radiusEstimateCandidate <= min(dim(mat)/2)){
      mat_radiusEstimate <- mat_radiusEstimateCandidate
    }
  }

  return(list("centerEstimate" = mat_bfRegionCenter,
              "radiusEstimate" = mat_radiusEstimate))
}

# Crop the exterior of a breech face impression surface matrix
#
# @name preProcess_cropExterior
#
# @param x3p an x3p object containing the surface matrix of a cartridge case
#   scan
# @param scheme argument for imager::imgradient
# @param high_connectivity argument for imager::label
# @param tolerance argument for imager::label
# @param radiusOffset number of pixels to add to estimated breech face radius.
#   This is commonly a negative value (e.g., -30) to trim the cartridge case
#   primer roll-off from the returned, cropped surface matrix.
# @param croppingThresh argument for cmcR::preProcess_cropWS
# @param agg_function the breech face radius estimation procedure returns
#   numerous radius estimates. This argument dictates the function used to
#   aggregate these into a final estimate.
#
# @return An x3p object containing the surface matrix of a breech face
#   impression scan where the rows/columns on the exterior of the breech face
#   impression have been cropped-out.
#
# @keywords internal
# @description The radius estimation procedure tends to over-estimate the
#   desired radius values. As such, a lot of the breech face impression
#   "roll-off" is included in the final scan. Excessive roll-off can bias the
#   calculation of the CCF. As such, we can manually shrink the radius estimate
#   so that little to no roll-off is included in the final processed scan.

preProcess_cropExterior <- function(x3p,
                                    scheme = 3,
                                    high_connectivity = FALSE,
                                    tolerance = 0,
                                    radiusOffset = 0,
                                    croppingThresh = 1,
                                    roughEstimate = FALSE,
                                    agg_function = median,
                                    ...){
  mat <- x3p$surface.matrix

  mat_estimates <- estimateBFRadius(mat = mat,
                                    scheme = scheme,
                                    high_connectivity = high_connectivity,
                                    tolerance = tolerance,
                                    angle = 0,
                                    roughEstimate = roughEstimate,
                                    agg_function = agg_function)

  if(all(!is.na(mat_estimates))){
    mat_radiusEstimate <- mat_estimates$radiusEstimate %>%
      magrittr::add(radiusOffset)

    mat_centerEstimateRow <- mat_estimates$centerEstimate[1]
    mat_centerEstimateCol <- mat_estimates$centerEstimate[2]
  }

  #the edges of some cartridge case scans aren't prominent, so the radius
  #estimate obtained above might not be accurate. The estimateBFRadius function
  #has built-in logic to determine if the radius estimates are obviously
  #incorrect (e.g., if the radius estimate is 0, then an NA is returned). We can
  #obtain a more precise estimate of the radius by rotating the matrix and
  #estimating the radius per rotation. Since this is computationally more
  #expensive, we only want to do this if necessary (i.e., if the initial radius
  #estimate came back NA).

  if(all(is.na(mat_estimates)) | !roughEstimate){
    mat_estimates <- purrr::map(seq(-180,180,by = 20),
                                ~ estimateBFRadius(mat = mat,
                                                   scheme = scheme,
                                                   high_connectivity = high_connectivity,
                                                   tolerance = tolerance,
                                                   angle = .,
                                                   roughEstimate = roughEstimate,
                                                   agg_function = agg_function)) %>%
      purrr::discard(.p = ~ all(is.na(.))) %>%
      purrr::compact()

    mat_radiusEstimate <- mat_estimates %>%
      purrr::map_dbl(~ .$radiusEstimate) %>%
      agg_function(na.rm = TRUE) %>%
      magrittr::add(radiusOffset)

    mat_centerEstimateRow <- mat_estimates %>%
      purrr::map_dbl(~ .$centerEstimate[1]) %>%
      agg_function(na.rm = TRUE)

    mat_centerEstimateCol <- mat_estimates %>%
      purrr::map_dbl(~ .$centerEstimate[2]) %>%
      agg_function(na.rm = TRUE)
  }

  optionalInput <- list(...)
  if(length(optionalInput[[1]]) > 0 & !all(unlist(optionalInput[[1]]) == -1)){
    mat_centerEstimateRow <- unlist(optionalInput)[[1]]
    mat_centerEstimateCol <- unlist(optionalInput)[[2]]
  }

  exteriorIndices <- expand.grid(row = 1:nrow(mat),col = 1:ncol(mat)) %>%
    dplyr::filter((row - mat_centerEstimateRow)^2 + (col - mat_centerEstimateCol)^2 > mat_radiusEstimate^2) %>%
    as.matrix()

  mat_interior <- mat
  mat_interior[exteriorIndices] <- NA

  x3p$surface.matrix <- mat_interior

  x3p <- preProcess_cropWS(x3p,croppingThresh = croppingThresh)

  return(x3p)
}

# Filter-out the firing pin impression region of a breech face impression scan
#
# @name preProcess_filterInterior
#
# @param x3p an x3p object containing the surface matrix of a cartridge case
#   scan
# @param scheme argument for imager::imgradient
# @param high_connectivity argument for imager::label
# @param tolerance argument for imager::label
# @param radiusOffset number of pixels to add to estimated firing pin hole
#   radius. This is commonly a positive value (e.g., 200) to not only remove
#   observations within the firing pin impression hole, but also the plateaued
#   region surrounding the hole that does not come into contact with the breech
#   face.
#
# @return An x3p object containing the surface matrix of a breech face
#   impression scan where the observations on the interior of the firing pin
#   impression hole have been filtered out.
# @keywords internal
# @note The radius estimation procedure effectively estimates the radius
#   of the firing pin hole. Unfortunately, it is often desired that more than
#   just observations in firing pin hole are removed. In particular, the
#   plateaued region surrounding the firing pin impression hole does not come
#   into contact with the breech face of a firearm and is thus unwanted in the
#   final, processed scan. The radiusOffset argument must be tuned (around 200
#   seems to work well for the Fadul cartridge cases) to remove these unwanted
#   observations.

preProcess_filterInterior <- function(x3p,
                                      scheme = 3,
                                      high_connectivity = FALSE,
                                      tolerance = 0,
                                      radiusOffset = 0,
                                      ...){
  mat <- x3p$surface.matrix

  mat_bfRegion <- mat

  mat_bfRegionBinarized <- mat_bfRegion
  mat_bfRegionBinarized[!is.na(mat_bfRegionBinarized)] <- 1
  mat_bfRegionBinarized[is.na(mat_bfRegionBinarized)] <- 0

  mat_bfRegionLabeled <- mat_bfRegionBinarized %>%
    imager::as.cimg() %>%
    imager::imgradient(scheme = scheme) %>%
    imager::enorm() %>%
    imager::add() %>%
    imager::label(high_connectivity = high_connectivity,
                  tolerance = tolerance) %>%
    as.matrix()

  mat_bfRegioncenterLabel <- mat_bfRegionLabeled[round(nrow(mat_bfRegionLabeled)/2),round(ncol(mat_bfRegionLabeled)/2)]

  mat_bfRegionfpHoleIndices <- which(mat_bfRegionLabeled == mat_bfRegioncenterLabel,arr.ind = TRUE)

  mat_bfRegionfpHoleCenter <- round(colMeans(mat_bfRegionfpHoleIndices))

  mat_bfRegionfpHoleRadiusEstim <- {mat_bfRegionLabeled == mat_bfRegioncenterLabel} %>%
    sum() %>%
    magrittr::divide_by(pi) %>%
    sqrt() %>%
    round() %>%
    magrittr::add(radiusOffset)

  optionalInput <- list(...)
  if(length(optionalInput[[1]]) > 0 & !all(unlist(optionalInput[[1]]) == -1) & !is.null(x3p$cmcR.info$exteriorCroppingValues)){
    mat_bfRegionfpHoleCenter[1] <- unlist(optionalInput)[[1]] - x3p$cmcR.info$exteriorCroppingValues[1]
    mat_bfRegionfpHoleCenter[2] <- unlist(optionalInput)[[2]] - x3p$cmcR.info$exteriorCroppingValues[2]
  }

  interiorIndices <- expand.grid(row = 1:nrow(mat_bfRegion),col = 1:ncol(mat_bfRegion)) %>%
    dplyr::filter((row - mat_bfRegionfpHoleCenter[1])^2 + (col - mat_bfRegionfpHoleCenter[2])^2 <= mat_bfRegionfpHoleRadiusEstim^2) %>%
    as.matrix()

  mat_bfRegionFiltered <- mat_bfRegion
  mat_bfRegionFiltered[interiorIndices] <- NA

  x3p_clone <- x3p
  x3p_clone$surface.matrix <- mat_bfRegionFiltered
  x3p_clone$header.info$sizeY <- ncol(mat_bfRegionFiltered)
  x3p_clone$header.info$sizeX <- nrow(mat_bfRegionFiltered)


  return(x3p_clone)
}

# #' Remove observations from the exterior of interior of a breech face scan
# #'
# #' @name preProcess_crop
# #'
# #' @param x3p an x3p object containing the surface matrix of a cartridge case
# #'  scan
# #' @param region dictates whether the observations on the "exterior" or
# #'  "interior" of the scan are removed
# #' @param radiusOffset the procedures used to detect the interior/exterior of the
# #'  scan tend to under/overestimate the radius of the region of interest,
# #'  respectively. The number passed to this argument is added to the radius
# #'  estimate as an offset.
# #' @param croppingThresh minimum number of non-NA pixels that need to be in a
# #'  row/column for it to not be cropped out of the breech face scan exterior
# #' @param agg_function the breech face radius estimation procedure returns a
# #'  number of radius estimates. This argument dictates the function used to
# #'  aggregate these into a final estimate.
# #' @param scheme argument for imager::imgradient
# #' @param high_connectivity argument for imager::label
# #' @param tolerance argument for imager::label
# #' @param radiusOffset number of pixels to add to estimated breech face radius.
# #'  This is commonly a negative value (e.g., -30 for region = "exterior") to
# #'  trim the cartridge case primer roll-off from the returned, cropped surface
# #'  matrix or a positive value (e.g., 200 for region = "interior") to remove
# #'  observations around the firing pin impression hole.
# #'
# #' @return An x3p object containing the surface matrix of a breech face
# #'  impression scan where the observations on the exterior/interior of the
# #'  breech face scan surface.
# #'
# #' @note The radius estimation procedure tends to over-estimate the desired
# #'  radius values. As such, a lot of the breech face impression "roll-off" is
# #'  included in the final scan. Excessive roll-off can bias the calculation of
# #'  the CCF. As such, we can manually shrink the radius estimate (-30 or -30
# #'  seems to work well for the Fadul cartridge cases) so that little to no
# #'  roll-off is included in the final processed scan.
# #'
# #' @note The radius estimation procedure is effective at estimating the radius of
# #'  the firing pin hole. Unfortunately, it is often desired that more than just
# #'  observations in firing pin hole are removed. In particular, the plateaued
# #'  region surrounding the firing pin impression hole does not come into contact
# #'  with the breech face of a firearm and is thus unwanted in the final,
# #'  processed scan. The radiusOffset argument must be tuned (around 200 seems to
# #'  work well for the Fadul cartridge cases) to remove these unwanted
# #'  observations.
# #' @examples
# #'
# #' #Process fadul1.1 "from scratch" (takes > 5 seconds to run)
# #' \dontrun{
# #' nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
# #' fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
# #'
# #' fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))
# #'
# #' fadul1.1_extCropped <- preProcess_crop(x3p = fadul1.1,
# #'                                        radiusOffset = -30,
# #'                                        region = "exterior")
# #'
# #' fadul1.1_extIntCropped <- preProcess_crop(x3p = fadul1.1_extCropped,
# #'                                           radiusOffset = 200,
# #'                                           region = "interior")
# #'
# #' x3pListPlot(list("Original" = fadul1.1,
# #'                  "Exterior Cropped" = fadul1.1_extCropped,
# #'                  "Exterior & Interior Cropped" = fadul1.1_extIntCropped ))
# #' }
# #' @export

# preProcess_crop <- function(x3p,
#                             region = "exterior",
#                             radiusOffset = 0,
#                             croppingThresh = 1,
#                             agg_function = median,
#                             scheme = 3,
#                             high_connectivity = FALSE,
#                             tolerance = 0,
#                             roughEstimateExterior = FALSE,
#                             ...){
#   #test that region is "exterior" or "interior"
#
#   optionalInput <- list(...)
#
#   if(region == "exterior"){
#     x3p <- preProcess_cropExterior(x3p = x3p,
#                                    radiusOffset = radiusOffset,
#                                    high_connectivity = high_connectivity,
#                                    tolerance = tolerance,
#                                    croppingThresh = croppingThresh,
#                                    roughEstimate = roughEstimateExterior,
#                                    agg_function = agg_function,
#                                    scheme = scheme,
#                                    optionalInput)
#
#     return(x3p)
#   }
#   if(region == "interior"){
#     x3p <- preProcess_filterInterior(x3p,
#                                      radiusOffset = radiusOffset,
#                                      high_connectivity = high_connectivity,
#                                      tolerance = tolerance,
#                                      scheme = scheme,
#                                      optionalInput)
#
#     return(x3p)
#   }
# }

#' Level a breech face impression surface matrix by a conditional statistic
#'
#' @name preProcess_removeTrend
#' @param x3p an x3p object containing the surface matrix of a cartridge case
#'   scan
#' @param statistic either "mean" or "quantile"
#' @param ... arguments to be set in the quantreg::rq function if statistic =
#'   "quantile" is set. In this case, tau = .5 and method = "fn" are recommended
#' @return an x3p object containing the leveled cartridge case scan surface
#'   matrix.
#' @examples
#'
#' #Process fadul1.1 "from scratch" (takes > 5 seconds to run)
#' \dontrun{
#' nbtrd_link <- "https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/"
#' fadul1.1_link <- "DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d"
#'
#' fadul1.1 <- x3ptools::read_x3p(paste0(nbtrd_link,fadul1.1_link))
#' fadul1.1_extCropped <- preProcess_crop(x3p = fadul1.1,
#'                                        region = "exterior",
#'                                        radiusOffset = -30)
#'
#' fadul1.1_intCroped <- preProcess_crop(x3p = fadul1.1_extCropped,
#'                                       region = "interior",
#'                                       radiusOffset = 200)
#'
#' fadul1.1_leveled <- preProcess_removeTrend(x3p = fadul1.1_intCroped,
#'                                            statistic = "quantile",
#'                                            tau = .5,
#'                                            method = "fn")
#'
#' x3pListPlot(list("Original" = fadul1.1,
#'                  "Ext. Cropped" = fadul1.1_extCropped,
#'                  "Ext. & Int. Cropped" = fadul1.1_intCroped,
#'                  "Cropped and Leveled" = fadul1.1_leveled))
#'}
#' @export

preProcess_removeTrend <- function(x3p,
                                   statistic = "mean",
                                   ...){
  stopifnot(statistic %in% c("quantile","mean"))

  if(statistic == "quantile"){
    rqArgs <- list(...)
    stopifnot("tau" %in% names(rqArgs) & "method" %in% names(rqArgs))

    x3p_fit <- quantreg::rq(data = expand.grid(y = 1:nrow(x3p$surface.matrix),
                                               x = 1:ncol(x3p$surface.matrix)) %>%
                              dplyr::mutate(value = as.numeric(x3p$surface.matrix)),
                            formula = value ~ x + y,
                            tau = rqArgs$tau,
                            method = rqArgs$method)
  }

  else if(statistic == "mean"){
    surfaceMat <- x3p$surface.matrix

    observedPixelLocations <- data.frame(which(!is.na(surfaceMat),
                                               arr.ind = TRUE)) %>%
      dplyr::mutate(depth = surfaceMat[!is.na(surfaceMat)])

    x3p_fit <- lm(depth ~ row + col,
                  data = observedPixelLocations)
  }

  x3p_condStatRemoved <- x3p
  x3p_condStatRemoved$surface.matrix[!is.na(x3p_condStatRemoved$surface.matrix)] <- x3p_condStatRemoved$surface.matrix[!is.na(x3p_condStatRemoved$surface.matrix)] - x3p_fit$fitted.values

  return(x3p_condStatRemoved)
}

