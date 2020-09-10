#' @name estimateBFRadius
#' @keywords internal

estimateBFRadius <- function(mat,
                             scheme = 3,
                             high_connectivity = FALSE,
                             tolerance = 0,
                             angle = 0,
                             interpolation = 0,
                             boundary = 0,
                             agg.function = median){
  # print(angle)
  # if(angle == -100){
  # browser()
  # }

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
    dplyr::group_by(firstRow) %>%
    dplyr::tally() %>%
    dplyr::top_n(n = 1,wt= n) %>%
    dplyr::pull(firstRow)

  mat_segmented <- mat_labeled
  mat_segmented[mat_segmented != exteriorLabel] <- -1
  mat_segmented[mat_segmented == exteriorLabel] <- -2

  mat_segmentedEdges <- mat_segmented %>%
    imager::as.cimg() %>%
    imager::imgradient(axes = "xy",
                       scheme = scheme) %>%
    imager::enorm() %>%
    imager::add() %>%
    as.matrix()

  mat_segmentedEdgesMidRow <- mat_segmentedEdges[floor(nrow(mat_segmentedEdges)/2),]

  if(sum(mat_segmentedEdgesMidRow > 0) < 2){
    return(NA)
  }

  mat_radiusEstimate <- data.frame(y = mat_segmentedEdgesMidRow) %>%
    dplyr::mutate(x = 1:nrow(.)) %>%
    dplyr::filter(y != 0) %>%
    dplyr::mutate(x_lag = c(x[2:(nrow(.))],NA)) %>%
    dplyr::mutate(x_diff = abs(x - x_lag)) %>%
    dplyr::top_n(n = 1,wt = x_diff) %>%
    dplyr::summarise(x_lag = agg.function(x_lag),
                     x = agg.function(x)) %>%
    dplyr::summarise(radEstimate = round((x_lag - x)/2)) %>%
    dplyr::pull(radEstimate) %>%
    agg.function(na.rm = TRUE)

  if(all(2*mat_radiusEstimate < max(nrow(mat)/2,ncol(mat)/2))){
    return(NA)
  }

  return(mat_radiusEstimate)
}

#' @name cropBFExterior
#' @export
#' @description The radius estimation procedure tends to over-estimate the
#'   desired radius values. As such, a lot of the breech face impression
#'   "roll-off" is included in the final scan. Excessive roll-off can bias the
#'   calculation of the CCF. As such, we can manually shrink the radius estimate
#'   so that little to no roll-off is included in the final processed scan.

cropBFExterior <- function(x3p,
                           scheme = 3,
                           high_connectivity = FALSE,
                           tolerance = 0,
                           radiusOffset = 0,
                           croppingThresh = 1,
                           agg.function = median){
  mat <- x3p$surface.matrix

  mat_radiusEstimate <- estimateBFRadius(mat = mat,
                                         scheme = scheme,
                                         high_connectivity = high_connectivity,
                                         tolerance = tolerance,
                                         angle = 0,
                                         agg.function = agg.function) %>%
    magrittr::add(radiusOffset)

  #the edges of some cartridge case scans aren't prominent, so the radius
  #estimate obtained above might not be accurate. The estimateBFRadius function
  #has built-in logic to determine if the radius estimates are obviously
  #incorrect (e.g., if the radius estimate is 0, then an NA is returned). We can
  #obtain a more precise estimate of the radius by rotating the matrix and
  #estimating the radius per rotation. Since this is computationally more
  #expensive, we only want to do this if necessary (i.e., if the initial radius
  #estimate came back NA).
  if(is.na(mat_radiusEstimate)){
    mat_radiusEstimate <- purrr::map_dbl(seq(-180,180,by = 20),
                                         ~ estimateBFRadius(mat = mat,
                                                            scheme = scheme,
                                                            high_connectivity = high_connectivity,
                                                            tolerance = tolerance,
                                                            angle = .,
                                                            agg.function = agg.function)) %>%
      agg.function(na.rm = TRUE) %>%
      magrittr::add(radiusOffset)
  }

  exteriorIndices <- expand.grid(row = 1:nrow(mat),col = 1:ncol(mat)) %>%
    dplyr::filter((row - nrow(mat)/2)^2 + (col - ncol(mat)/2)^2 > mat_radiusEstimate^2) %>%
    as.matrix()

  mat_interior <- mat
  mat_interior[exteriorIndices] <- NA

  mat_cropped <- mat_interior %>%
    cmcR::preProcess_cropWS(croppingThresh = croppingThresh)

  x3p_clone <- x3p
  x3p_clone$surface.matrix <- mat_cropped
  x3p_clone$header.info$sizeY <- ncol(mat_cropped)
  x3p_clone$header.info$sizeX <- nrow(mat_cropped)


  return(x3p_clone)
}

#' @name filterBFInterior
#' @export

filterBFInterior <- function(x3p,
                             scheme = 3,
                             high_connectivity = FALSE,
                             tolerance = 0,
                             radiusOffset = 0){
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

#' @name levelByConditionalStatistic
#' @description For statistic = "quantile," tau = .5 and method = "fn" recommended by default
#' @export

levelByConditionalStatistic <- function(x3p,
                                        statistic = "mean",
                                        ...){
  stopifnot(statistic %in% c("quantile","mean"))

  if(statistic == "quantile"){
    rqArgs <- list(...)
    stopifnot("tau" %in% names(rqArgs) & "method" %in% names(rqArgs))

    x3p_fit <- quantreg::rq(data = expand.grid(y = 1:nrow(x3p$surface.matrix),
                                               x = 1:ncol(x3p$surface.matrix)) %>%
                              mutate(value = as.numeric(x3p$surface.matrix)),
                            formula = value ~ x + y,
                            tau = rqArgs$tau,
                            method = rqArgs$method)
  }

  else if(statistic == "mean"){
    x3p_fit <- lm(data = expand.grid(y = 1:nrow(x3p$surface.matrix),
                                     x = 1:ncol(x3p$surface.matrix)) %>%
                    mutate(value = as.numeric(x3p$surface.matrix)),
                  formula = value ~ x + y)
  }

  x3p_condStatRemoved <- x3p
  x3p_condStatRemoved$surface.matrix[!is.na(x3p_condStatRemoved$surface.matrix)] <- x3p_condStatRemoved$surface.matrix[!is.na(x3p_condStatRemoved$surface.matrix)] - x3p_fit$fitted.values

  return(x3p_condStatRemoved)
}
