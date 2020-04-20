#' @name fpRadiusEstimation
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param sum_over a string, either "row" or "col" indicating whether the
#'   non-NAs pixels should be summed over rows or over cols of the surface
#'   matrix
#' @param nonNA_sum_smootherSize size of average smoother (to be passed to
#'   zoo::roll_mean). Should be an odd.
#'
#' @description An x3p file contains a "surface matrix" whose elements represent
#'   individual pixels of a breech face impression microscopy scan. The
#'   dimension of the surface matrix indicates the width/depth (x,y, if viewed
#'   from a top-down perspective) of the scan while the actual values of each
#'   element of the surface matrix represent the height (z) of the breech face
#'   at that observed (x,y) location. So a particular "voxel" is characterized
#'   by its row/column index in the matrix along with the value of the matrix at
#'   that particular row/column combination. An element of the surface matrix
#'   may be NA if the scan did not detect the cartridge case at the particular
#'   (x,y) location. For example, most of the values on the perimeter of the
#'   surface matrix are typically NA because the scan itself is taken such that
#'   the breech face impression is approximately centered in the scanning area
#'   (i.e, there's a fair amount of breathing room between the actual breech
#'   face impression and the edges of the scan). All of this introduction is to
#'   explain that for a given row/col in the surface matrix we can add up the
#'   number of non-NA pixels and use this as an indicator for where the breech
#'   face impression sits in the image (e.g, the first row in which a non-NA
#'   pixel occurs is likely the beginning of the top of the breech face
#'   impression while the last row in which a non-NA pixel occurs is likely the
#'   bottom of the breech face impression). While it's fairly easy to detect the
#'   beginning/end of the breech face impression within the entire surface
#'   matrix, our goal is to determine the location of the firing pin impression
#'   circle in the surface matrix. The firing pin impression circle is typically
#'   a large circle punched through the interior of the breech face impression
#'   scan. The goal of the function below is to detect where the firing pin
#'   circle begins/ends in the surface matrix by considering the non-NA pixel
#'   counts per row (or col). The idea being that the two largest modes of this
#'   row/col-wise pixel count sequence should be the beginning/end of the firing
#'   pin impression circle. Identifying the distance between the two largest
#'   modes is essentially all that the function below accomplishes. It is often
#'   the case that a breech face impression scan has a few spotty missing values
#'   within the interior of the breech face impression (electron microscpy isn't
#'   perfect). As such, the row/col-wise non-NA pixel counts are not often a
#'   smooth sequence. Thus, the zoo::roll_mean function is used to smooth-out
#'   the non-NA pixel counts so that the trend in the counts (in particular, the
#'   modes) can be differentiated from noise introduced by spotty missing
#'   values. The size of the kernel used in the roll_mean function is
#'   indicated by the nonNA_sum_smootherSize argument.
#'   @keywords internal

fpRadiusEstimation <- function(surfaceMat,
                               sum_over = "row",
                               nonNA_sum_smootherSize){
  if(sum_over == "col"){
    surfaceMat_summedNAs <- surfaceMat %>%
      is.na() %>%
      colSums(na.rm = TRUE)
  }
  else{
    surfaceMat_summedNAs <- surfaceMat %>%
      is.na() %>%
      rowSums(na.rm = TRUE)
  }

  firingPinDiameter_estim <- surfaceMat_summedNAs %>%
    data.frame(naSum = .) %>%
    dplyr::mutate(nonNA_sum = nrow(surfaceMat) - naSum,
                  smoothednonNA_sum =c(rep(NA,nonNA_sum_smootherSize/2),
                                       zoo::rollmean(nonNA_sum,k = nonNA_sum_smootherSize),
                                       rep(NA,nonNA_sum_smootherSize/2)),
                  rowNum = seq(length.out = length(nonNA_sum)),
                  smoothednonNA_sum_diff = c(NA,diff(smoothednonNA_sum)),
                  smoothednonNA_sum_diff_lag = c(smoothednonNA_sum_diff[2:length(smoothednonNA_sum_diff)],NA)) %>%
    dplyr::filter(smoothednonNA_sum_diff >= 0 & smoothednonNA_sum_diff_lag < 0) %>% #search for rows that change from increasing to decreasing (local maxima)
    dplyr::slice(c(1,nrow(.))) %>% #the first and last of such local maxima rows should be the start and end of the firing pin impression circle
    dplyr::pull(rowNum) %>%
    diff()

  return(floor(firingPinDiameter_estim/2))
}

#' @name fpRadiusGridSearch
#'
#' @param surfaceMat a surface matrix representing a breech face impression scan
#' @param smootherSize size of average smoother (to be passed to zoo::roll_mean)
#' @param aggregation_function function to select initial radius estimate from
#'   those calculated using fpRadiusGridSearch
#'
#' @description This function effectively just repeatedly calls the
#'   fpRadiusEstimation function to obtain multiple estimates
#'   of the firing pin radius. It calculates row/col-wise radius estimates for 6
#'   different rotations of the surface matrix (0, 15, 30, 45, 60, and 75
#'   degrees). Based on the argument passed to aggregation_function, these 12
#'   radii estimates are reduced to a single, rough radius estimate (e.g.,
#'   minimum was determined to be an effective aggregation function in
#'   preliminary tests).
#'
#' @keywords internal

fpRadiusGridSearch <- function(surfaceMat,
                               smootherSize,
                               aggregation_function){
  if(!is.matrix(surfaceMat)){
    surfaceMat <- as.matrix(surfaceMat)
  }

  estim_rotated <- purrr::map(c(0,15,30,45,60,75),
                              function(theta){
                                if(theta != 0){
                                  surfaceMatRotated <- rotateSurfaceMatrix(surfaceMat,
                                                                           theta)
                                }
                                else{
                                  surfaceMatRotated <- surfaceMat
                                }

                                estimByRow_rotated <-
                                  cmcR:::fpRadiusEstimation(surfaceMatRotated,
                                                     sum_over = "row",
                                                     nonNA_sum_smootherSize = smootherSize)

                                estimByCol_rotated <-
                                  cmcR:::fpRadiusEstimation(surfaceMatRotated,
                                                     sum_over = "col",
                                                     nonNA_sum_smootherSize = smootherSize)

                                return(c("byRow" = estimByRow_rotated,
                                         "byCol" = estimByCol_rotated))
                              }) %>%
    setNames(c("0","15","30","45","60","75"))

  radiusEstim <- estim_rotated %>%
    unlist() %>%
    aggregation_function(na.rm = TRUE) %>%
    floor()

  return(list("radiusEstim" = radiusEstim,
              "individualEstims" = unlist(estim_rotated)))
}

#' Non-maxima suppression
#'
#' @name nms
#'
#' @seealso #https://www.rdocumentation.org/packages/imager/versions/0.41.2/topics/hough_circle
#' @keywords internal

nms <- function(im,sigma){
  #non-maxima suppression

  im[imager::dilate_square(im,sigma) != im] <- 0

  return(im)
}


