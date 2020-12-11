# Apply CMC classification logic of the original method of Song (2013)
#
# @name decision_originalMethod_classifyCMCs
# @param cellIndex vector/tibble column containing cell indices corresponding
#   to a reference cell
# @param x vector/tibble column containing x horizontal translation values
# @param y vector/tibble column containing y vertical translation values
# @param theta vector/tibble column containing theta rotation values
# @param corr vector/tibble column containing correlation similarity scores
#   between a reference cell and its associated target region
# @param xThresh used to classify particular x values "congruent" if they are
#   within xThresh of the median x value
# @param yThresh used to classify particular y values "congruent" if they are
#   within yThresh of the median y value
# @param thetaThresh used to classify particular theta values "congruent" if
#   they are within thetaThresh of the median theta value
# @param corrThresh to classify particular correlation values "congruent" if
#   they are at least corrThresh
#
# @seealso \url{https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193}
#@note the decision_CMC function internally calls this function if a value of
#   tau is not provided
#
# @keywords internal

decision_originalMethod_classifyCMCs <- function(cellIndex,
                                                 x,
                                                 y,
                                                 theta,
                                                 corr,
                                                 xThresh = 20,
                                                 yThresh = xThresh,
                                                 thetaThresh = 6,
                                                 corrThresh = .5){

  comparisonFeaturesDF <- data.frame(cellIndex = cellIndex,
                                     x = x,
                                     y = y,
                                     theta = theta,
                                     corr = corr)

  originalMethodCMCs <- comparisonFeaturesDF %>%
    dplyr::group_by(cellIndex) %>%
    dplyr::top_n(n = 1,wt = corr) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(originalMethodClassif = ifelse(abs(x - median(x)) <= xThresh &
                                                   abs(y - median(y)) <= yThresh &
                                                   abs(theta - median(theta)) <= thetaThresh &
                                                   corr >= corrThresh,"CMC","non-CMC")) %>%
    dplyr::filter(originalMethodClassif == "CMC") %>%
    dplyr::select(cellIndex,theta,originalMethodClassif)

  originalMethodClassif <- comparisonFeaturesDF %>%
    dplyr::left_join(originalMethodCMCs,by = c("cellIndex","theta")) %>%
    dplyr::mutate(originalMethodClassif = ifelse(is.na(originalMethodClassif),"non-CMC","CMC")) %>%
    dplyr::pull(originalMethodClassif)

  return(originalMethodClassif)
}

#'Compute CMC-theta distribution for a set of comparison features
#'
#'@name decision_highCMC_cmcThetaDistrib
#'@param cellIndex vector/tibble column containing cell indices corresponding to
#'  a reference cell
#'@param x vector/tibble column containing x horizontal translation values
#'@param y vector/tibble column containing y vertical translation values
#'@param theta vector/tibble column containing theta rotation values
#'@param corr vector/tibble column containing correlation similarity scores
#'  between a reference cell and its associated target region
#'@param xThresh used to classify particular x values "congruent" (conditional
#'  on a particular theta value) if they are within xThresh of the
#'  theta-specific median x value
#'@param yThresh used to classify particular y values "congruent" (conditional
#'  on a particular theta value) if they are within yThresh of the
#'  theta-specific median y value
#'@param corrThresh to classify particular correlation values "congruent"
#'  (conditional on a particular theta value) if they are at least corrThresh
#'@return a vector of the same length as the input containing a "CMC Candidate"
#'  or "Non-CMC Candidate" classification based on whether the particular
#'  cellIndex has congruent x,y, and theta features.
#'@note This function is a helper internally called in the decision_CMC
#'  function. It is exported to be used as a diagnostic tool for the High CMC
#'  method
#'
#'@examples
#'\dontrun{
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'comparisonDF <- purrr::map_dfr(seq(-30,30,by = 3),
#'                               ~ comparison_allTogether(fadul1.1_processed,
#'                                                        fadul1.2_processed,
#'                                                        theta = .))
#'
#'comparisonDF <- comparisonDF %>%
#' dplyr::mutate(cmcThetaDistribClassif = decision_highCMC_cmcThetaDistrib(cellIndex = cellIndex,
#'                                                                  x = x,
#'                                                                  y = y,
#'                                                                  theta = theta,
#'                                                                  corr = pairwiseCompCor))
#'
#'comparisonDF %>%
#' dplyr::filter(cmcThetaDistribClassif == "CMC Candidate") %>%
#' ggplot2::ggplot(ggplot2::aes(x = theta)) +
#' ggplot2::geom_bar(stat = "count")
#' }
#'@importFrom rlang .data
#'@export

decision_highCMC_cmcThetaDistrib <- function(cellIndex,
                                             x,
                                             y,
                                             theta,
                                             corr,
                                             xThresh = 20,
                                             yThresh = xThresh,
                                             corrThresh = .5){

  comparisonFeaturesDF <- data.frame(cellIndex = cellIndex,
                                     x = x,
                                     y = y,
                                     theta = theta,
                                     corr = corr)

  highCMC_candidates <- comparisonFeaturesDF %>%
    dplyr::group_by(theta) %>%
    dplyr::mutate(cmcThetaDistribClassif = ifelse(abs(x - median(x)) <= xThresh &
                                                    abs(y - median(y)) <= yThresh &
                                                    corr >= corrThresh,
                                                  "CMC Candidate","not CMC Candidate")) %>%
    dplyr::ungroup() %>%
    dplyr::select(cellIndex,theta,cmcThetaDistribClassif)


  cmcThetaDistribClassif <- comparisonFeaturesDF %>%
    dplyr::left_join(highCMC_candidates,by = c("cellIndex","theta")) %>%
    dplyr::pull(cmcThetaDistribClassif)

  return(cmcThetaDistribClassif)
}

#'Classify theta values in CMC-theta distribution as having "High" or "Low" CMC
#'candidate counts
#'
#'@name decision_highCMC_identifyHighCMCThetas
#'
#'@param cmcThetaDistrib output of the decision_highCMC_cmcThetaDistrib function
#'@param tau constant used to define a "high" CMC count. This number is
#'  subtracted from the maximum CMC count achieved in the CMC-theta
#'  distribution. Theta values with CMC counts above this value are considered
#'  to have "high" CMC counts.
#'@return A vector of the same length as the input containing "High" or "Low"
#'  classification based on whether the associated theta value has a High CMC
#'  Candidate count.
#'@note This function is a helper internally called in the decision_CMC
#'  function. It is exported to be used as a diagnostic tool for the High CMC
#'  method
#'@examples
#'\dontrun{
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'comparisonDF <- purrr::map_dfr(seq(-30,30,by = 3),
#'                               ~ comparison_allTogether(fadul1.1_processed,
#'                                                        fadul1.2_processed,
#'                                                        theta = .))
#'
#'highCMCthetas <- comparisonDF %>%
#' dplyr::mutate(cmcThetaDistribClassif = decision_highCMC_cmcThetaDistrib(cellIndex = cellIndex,
#'                                                                  x = x,
#'                                                                  y = y,
#'                                                                  theta = theta,
#'                                                                  corr = pairwiseCompCor)) %>%
#' decision_highCMC_identifyHighCMCThetas(tau = 1)
#'
#'
#'highCMCthetas %>%
#' dplyr::filter(cmcThetaDistribClassif == "CMC Candidate") %>%
#' ggplot2::ggplot(ggplot2::aes(x = theta,fill = thetaCMCIdentif)) +
#' ggplot2::geom_bar(stat = "count")
#' }
#'@importFrom rlang .data
#'@export

decision_highCMC_identifyHighCMCThetas <- function(cmcThetaDistrib,
                                                   tau = 1){

  thetaClassifications <- cmcThetaDistrib %>%
    dplyr::filter(.data$cmcThetaDistribClassif == "CMC Candidate") %>%
    dplyr::group_by(.data$theta) %>%
    dplyr::tally(name = "cmcCandidateCount") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(thetaCMCIdentif = ifelse(.data$cmcCandidateCount >= max(.data$cmcCandidateCount) - tau,"High","Low"))

  cmcThetaDistrib %>%
    dplyr::left_join(thetaClassifications,by = "theta")
}

#Apply CMC classification logic of the Tong et al. (2015) to the CMC-theta
#distribution returned by the decision_highCMC_cmcThetaDistrib function
#
#@name decision_highCMC_classifyCMCs
#@param cellIndex vector/tibble column containing cell indices corresponding to
#  a reference cell
#@param x vector/tibble column containing x horizontal translation values
#@param y vector/tibble column containing y vertical translation values
#@param theta vector/tibble column containing theta rotation values
#@param corr vector/tibble column containing correlation similarity scores
#  between a reference cell and its associated target region
#@param xThresh used to classify particular x values "congruent" (conditional
#  on a particular theta value) if they are within xThresh of the
#  theta-specific median x value
#@param yThresh used to classify particular y values "congruent" (conditional
#  on a particular theta value) if they are within yThresh of the
#  theta-specific median y value
#@param thetaThresh defines how wide a High CMC mode is allowed to be in the
#  CMC-theta distribution before it's considered too diffuse
#@param corrThresh to classify particular correlation values "congruent"
#  (conditional on a particular theta value) if they are at least corrThresh
#@param tau constant used to define a "high" CMC count. This number is
#  subtracted from the maximum CMC count achieved in the CMC-theta
#  distribution. Theta values with CMC counts above this value are considered
#  to have "high" CMC counts.
#@seealso
#\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#
#@note the decision_CMC function internally calls this function if a value of
#  tau is provided
#@importFrom rlang .data
# @keywords internal

decision_highCMC_classifyCMCs <- function(cellIndex,
                                          x,
                                          y,
                                          theta,
                                          corr,
                                          xThresh = 20,
                                          yThresh = xThresh,
                                          thetaThresh = 6,
                                          corrThresh = .5,
                                          tau = 1){

  comparisonFeaturesDF <- data.frame(cellIndex = cellIndex,
                                     x = x,
                                     y = y,
                                     theta = theta,
                                     corr = corr)

  cmcThetaDistrib <- comparisonFeaturesDF %>%
    dplyr::mutate(cmcThetaDistribClassif = decision_highCMC_cmcThetaDistrib(cellIndex = cellIndex,
                                                                            x = x,
                                                                            y = y,
                                                                            theta = theta,
                                                                            corr = corr,
                                                                            xThresh = xThresh,
                                                                            yThresh = yThresh,
                                                                            corrThresh = corrThresh)) %>%
    dplyr::filter(.data$cmcThetaDistribClassif == "CMC Candidate")

  if(nrow(cmcThetaDistrib) == 0){
    highCMCClassif <- comparisonFeaturesDF %>%
      dplyr::mutate(highCMCClassif = "non-CMC (failed)") %>%
      dplyr::pull(highCMCClassif)

    return(highCMCClassif)
  }

  cmcThetaDistrib_classified <- decision_highCMC_identifyHighCMCThetas(cmcThetaDistrib,
                                                                       tau = tau)

  passesHighCMCCriterion <- cmcThetaDistrib_classified %>%
    dplyr::filter(.data$thetaCMCIdentif == "High") %>%
    dplyr::select(theta) %>%
    dplyr::distinct() %>%
    dplyr::summarise(distance = abs(max(theta) - min(theta))) %>%
    dplyr::pull(.data$distance) %>%
    {. <= thetaThresh}

  if(passesHighCMCCriterion){
    highCMCs <- cmcThetaDistrib_classified %>%
      dplyr::mutate(highCMCClassif = ifelse(.data$thetaCMCIdentif == "High","CMC","non-CMC (passed)")) %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter(highCMCClassif == "CMC") %>%
      dplyr::filter(corr == max(corr)) %>%
      dplyr::select(c(cellIndex,theta,highCMCClassif))

    highCMCClassif <- comparisonFeaturesDF %>%
      dplyr::left_join(highCMCs,by = c("cellIndex","theta")) %>%
      dplyr::mutate(highCMCClassif = ifelse(highCMCClassif == "non-CMC (passed)" | is.na(highCMCClassif),"non-CMC (passed)","CMC")) %>%
      dplyr::pull(highCMCClassif)
  }
  else{
    highCMCClassif <- comparisonFeaturesDF %>%
      dplyr::mutate(highCMCClassif = "non-CMC (failed)") %>%
      dplyr::pull(highCMCClassif)
  }

  return(highCMCClassif)
}

#'Applies the decision rules of the original method of Song (2013) or the High
#'CMC method of Tong et al. (2015)
#'
#'@name decision_CMC
#'@param cellIndex vector/tibble column containing cell indices corresponding to
#'  a reference cell
#'@param x vector/tibble column containing x horizontal translation values
#'@param y vector/tibble column containing y vertical translation values
#'@param theta vector/tibble column containing theta rotation values
#'@param corr vector/tibble column containing correlation similarity scores
#'  between a reference cell and its associated target region
#'@param xThresh used to classify particular x values "congruent" (conditional
#'  on a particular theta value) if they are within xThresh of the
#'  theta-specific median x value
#'@param yThresh used to classify particular y values "congruent" (conditional
#'  on a particular theta value) if they are within yThresh of the
#'  theta-specific median y value
#'@param thetaThresh (original method of Song (2013)) used to classify
#'  particular theta values "congruent" if they are within thetaThresh of the
#'  median theta value. (High CMC) defines how wide a High CMC mode is allowed
#'  to be in the CMC-theta distribution before it's considered too diffuse
#'@param corrThresh to classify particular correlation values "congruent"
#'  (conditional on a particular theta value) if they are at least corrThresh
#'@param tau (optional) parameter required to apply the High CMC method of Tong
#'  et al. (2015). If not given, then the decision rule of the original method
#'  of Song (2013) is applied. This number is subtracted from the maximum CMC
#'  count achieved in the CMC-theta distribution. Theta values with CMC counts
#'  above this value are considered to have "high" CMC counts.
#'
#'@return A vector of the same length as the input containing the CMC
#'  classification under one of the two decision rules.
#'@seealso \url{https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193}
#'@seealso
#'\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'@examples
#'\dontrun{
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'comparisonDF <- purrr::map_dfr(seq(-30,30,by = 3),
#'                               ~ comparison_allTogether(fadul1.1_processed,
#'                                                        fadul1.2_processed,
#'                                                        theta = .))
#'
#'
#'comparisonDF <- comparisonDF %>%
#' dplyr::mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
#'                                                    x = x,
#'                                                    y = y,
#'                                                    theta = theta,
#'                                                    corr = pairwiseCompCor),
#'               highCMCClassif = decision_CMC(cellIndex = cellIndex,
#'                                            x = x,
#'                                            y = y,
#'                                            theta = theta,
#'                                            corr = pairwiseCompCor,
#'                                            tau = 1))
#'
#'
#'comparisonDF %>%
#' dplyr::filter(originalMethodClassif == "CMC" | highCMCClassif == "CMC")
#' }
#'@export

decision_CMC <- function(cellIndex,
                         x,
                         y,
                         theta,
                         corr,
                         xThresh = 20,
                         yThresh = xThresh,
                         thetaThresh = 6,
                         corrThresh = .5,
                         tau = NULL){

  comparisonFeaturesDF <- data.frame(cellIndex = cellIndex,
                                     x = x,
                                     y = y,
                                     theta = theta,
                                     corr = corr)

  if(is.null(tau)){
    comparisonFeaturesDF <- comparisonFeaturesDF %>%
      dplyr::mutate(originalMethodClassif = decision_originalMethod_classifyCMCs(cellIndex = cellIndex,
                                                                                 x = x,
                                                                                 y = y,
                                                                                 theta = theta,
                                                                                 corr = corr,
                                                                                 xThresh = xThresh,
                                                                                 yThresh = yThresh,
                                                                                 thetaThresh = thetaThresh,
                                                                                 corrThresh = corrThresh))
    originalMethodClassif <- comparisonFeaturesDF %>%
      dplyr::pull(originalMethodClassif)

    return(originalMethodClassif)
  }

  if(is.numeric(tau)){
    highCMCClassif <- comparisonFeaturesDF %>%
      dplyr::mutate(highCMCClassif = decision_highCMC_classifyCMCs(cellIndex = cellIndex,
                                                                   x = x,
                                                                   y = y,
                                                                   theta = theta,
                                                                   corr = corr,
                                                                   xThresh = xThresh,
                                                                   yThresh = yThresh,
                                                                   thetaThresh = thetaThresh,
                                                                   corrThresh = corrThresh,
                                                                   tau = tau)) %>%
      dplyr::pull(highCMCClassif)

    return(highCMCClassif)
  }
}

# @name calcMaxCMCTheta
#
# @param distanceToCMCMaxTieBreaker decides what to do in cases where there are
#   consecutive theta values all tied for the max CMC count (no guidance is
#   given by Tong et al. (2015) of what to do in this situation). The default
#   is to determine the distance between any high CMC theta value and its
#   closest CMC max theta value.
#
#   For example, suppose 3 consecutive theta values each have the max CMC count
#   of 17. So the CMC counts around these theta values may look like ...10, 16,
#   17, 17, 17, 15, 12... For a highCMC_thresh = 1, note that there is one
#   theta value with associated count equal to the high CMC count of 17 - 1 =
#   16. The default setting of distanceToCMCMaxTieBreaker = "minDistance" will
#   consider the minimum distance between a max CMC theta value and this high
#   CMC theta. So in this case, because a theta value with CMC count 17 is
#   right next to this high CMC theta value, we would say that this particular
#   cartridge case pair "passes" the high CMC criterion, assuming the
#   thetaThresh is set to be equal to the grid spacing of the theta values
#   considered (3 degrees by default).
#
#   A slightly more "conservative" option would be to set
#   distanceToCMCMaxTieBreaker = "medDistance" take the median of these max CMC
#   theta values and apply the high CMC criterion. The example under
#   consideration *wouldn't* pass the high CMC criterion due to the particular
#   choice of theta grid spacing, thetaThresh, and highCMC_thresh. However,
#   different combinations of these arguments may lead to passing the
#   criterion.
#
#   Finally, the most "conservative" option would be to set
#   distanceToCMCMaxTieBreaker = "failCriterion" in which we immediately "fail"
#   a cartridge case pair if there are any ties for max CMC theta (i.e., "there
#   can only be one!" - Highlander).
#
# @keywords internal
#
# @importFrom stats median

utils::globalVariables(c("theta","n","distanceToCMCMax"))

calcMaxCMCTheta <- function(cmcPerTheta,
                            highCMC_thresh = 1,
                            thetaThresh = 3,
                            distanceToCMCMaxTieBreaker = "minDistance"){

  cmcCountPerTheta <- cmcPerTheta %>%
    dplyr::group_by(theta) %>%
    dplyr::tally()

  if(nrow(cmcCountPerTheta) == 0){
    return(NA)
  }

  cmcMax <- cmcCountPerTheta %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n))

  if(purrr::is_empty(cmcMax$theta) | nrow(cmcMax) == 0){
    return(NA)
  }

  #if there may be multiple theta values tied for cmcMax that are within
  #thetaThresh of each other, then we need to somehow "break" the ties. The
  #distanceToCMCMaxTieBreaker argument (see documentation) determines how to do
  #so.

  if(distanceToCMCMaxTieBreaker == "failCriterion" & nrow(cmcMax) > 1){
    return(NA)
  }

  #there may be more than one theta tied with the maximum CMC - in which case
  #we should determine if such thetas are "close" to each other, where
  #proximity is defined based on the theta_tresh set. If not, then we can rule
  #out that the comparison under scrutiny is a known match.
  if(any(diff(cmcMax$theta) > thetaThresh)){
    return(NA)
  }

  if(distanceToCMCMaxTieBreaker == "medDistance"){
    tieBreaker <- median
  }
  else if(distanceToCMCMaxTieBreaker == "minDistance"){
    tieBreaker <- min
  }

  maxDistancetoCMCMax <- cmcCountPerTheta %>%
    dplyr::filter(n >= unique(cmcMax$n) - highCMC_thresh) %>%
    dplyr::group_by(theta) %>%
    dplyr::summarise(distanceToCMCMax = tieBreaker(abs(cmcMax$theta - theta))) %>%
    dplyr::pull(distanceToCMCMax) %>%
    max()

  if(all(maxDistancetoCMCMax > thetaThresh)){
    return(NA)
  }
  else{
    if(distanceToCMCMaxTieBreaker == "medDistance"){
      return(median(cmcMax$theta))
    }
    return(median(cmcMax$theta))
  }
}

# Implements "improved" CMC logic on a list of CCF results for a comparison
# between two cartridge case scans as proposed by Tong et al. (2015)
#
# @name cmcFilter_improved
#
# @description Implements "improved' Congruent Matching Cells logic, as
#   proposed by Tong et al. (2015), to the CCF results of a comparison between
#   two cartridge case scans.
#
# @param consensus_function function to aggregate the translation (x and y)
#   and rotation (theta) values in the ccfDF data frame to determine
#   "consensus" values
# @param ccf_thresh minimum correlation threshold to call a cell pair
#   "congruent matching"
# @param dx_thresh maximum distance from the consensus x value that a cell
#   pair can be to be called "congruent matching"
# @param dy_thresh  maximum distance from the consensus y value that a cell
#   pair can be to be called "congruent matching"
# @param thetaThresh maximum distance from the consensus theta value that a
#   cell pair can be to be called "congruent matching"
# @param missingThetaDecision dictates how function should handle situations
#   in which one direction passes the high CMC criterion while another
#   direction does not. "replace": replaces theta value in failed direction
#   with opposite of theta value in successful direction. "dismiss": only
#   counts the initial CMCs in failed direction and high CMCs in successful
#   direction. "fail": only counts the initial CMCs in either direction.
# @param compareThetas dictates if the consensus theta values
#   determined under the initially proposed method should be compared to the
#   consensus theta values determined under the High CMC method. In particular,
#   determines for each direction whether the consensus theta values determined
#   under the two methods are within thetaThresh of each other. It is often
#   the case that non-matching cartridge cases, even if they pass the High CMC
#   criterion, will have differing consensus theta values under the two
#   methods. If this isn't taken into account, non-matches tend to be assigned
#   a lot of false positive CMCs under the High CMC method.
# @param consensus_function_theta *(OPTIONAL)* function (separate from
#   consensus_function) to aggregate the rotation (theta) values in the ccfDF
#   data frame to determine "consensus" value
#
# @seealso
# \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#
#@keywords internal
# @importFrom stats median
# @importFrom rlang .data

utils::globalVariables(c(".","cellIndex","comparison","theta"))

cmcFilter_improved <- function(reference_v_target_CMCs,
                               target_v_reference_CMCs,
                               corColName = "pairwiseCompCor",
                               thetaThresh = 6,
                               missingThetaDecision = "fail",
                               compareThetas = TRUE){

  originalCMCs_reference_v_target <- reference_v_target_CMCs %>%
    dplyr::filter(.data$originalMethodClassif == "CMC") %>%
    dplyr::select(-c(.data$originalMethodClassif,.data$highCMCClassif)) %>%
    dplyr::mutate(direction = "reference_v_target")

  originalCMCs_target_v_reference <- target_v_reference_CMCs %>%
    dplyr::filter(.data$originalMethodClassif == "CMC") %>%
    dplyr::select(-c(.data$originalMethodClassif,.data$highCMCClassif)) %>%
    dplyr::mutate(direction = "target_v_reference")

  smallerDirection <- which.min(c(nrow(originalCMCs_reference_v_target),nrow(originalCMCs_target_v_reference)))

  originalMethodCMCs <- list(originalCMCs_reference_v_target,originalCMCs_target_v_reference)

  highCMC_initial <- list(originalCMCs_reference_v_target,originalCMCs_target_v_reference)[[smallerDirection]]

  cmcPerTheta <- list(reference_v_target_CMCs,
                      target_v_reference_CMCs)

  #Important note: there may be multiple, consecutive theta values that tie for
  #the CMC max count. The default of the cmcR package is to take the median of
  #these theta values as the thetaMax value
  thetaMax <- purrr::map(cmcPerTheta,
                         function(cmc_direction_df){

                           if(suppressWarnings(
                             any(stringr::str_detect(cmc_direction_df,"(failed)"))
                           )){
                             return(NA)
                           }
                           else{
                             cmc_direction_df %>%
                               dplyr::filter(.data$highCMCClassif == "CMC") %>%
                               calcMaxCMCTheta(highCMC_thresh = 1,
                                               thetaThresh = thetaThresh,
                                               distanceToCMCMaxTieBreaker = "minDistance") %>%
                               return()
                           }
                         }) %>%
    purrr::set_names(c("reference_v_target","target_v_reference"))

  if(purrr::is_empty(thetaMax$reference_v_target) & purrr::is_empty(thetaMax$target_v_reference)){
    thetaMax$target_v_reference <- NA
    thetaMax$reference_v_target <- NA
  }

  #if neither direction passed the high CMC criterion...
  if(all(is.na(thetaMax))){
    return(list("originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = highCMC_initial))
  }

  #if one direction didn't pass the high CMC criterion...
  if(any(is.na(thetaMax))){

    #Note: getting the missingThetaDecision = "replace" option working would
    #require re-calculating the CMCs after replacing the missing theta value.
    #However, this would also require specifying the thresholds used to
    #determine the CMCs all over again, which is not required for this function
    #(nor do I think it should be). Perhaps we could use the ellipsis ... to
    #allow the user to input the thresholds if they set missingThetaDecision =
    #"replace", but for now we'll keep it commented-out

    # Most "liberal" decision is to replace the missing theta value with the
    # opposite of the other theta value
    # if(missingThetaDecision == "replace"){
    #   thetaMax[[which(is.na(thetaMax))]] <- -1*thetaMax[[which(!is.na(thetaMax))]]
    # }

    # More "moderate" decision is to dismiss the missing direction and only take
    # the initial CMCs defined for that direction
    # else if(missingThetaDecision == "dismiss"){
    if(missingThetaDecision == "dismiss"){
      #the direction that didn't pass gets assigned initial CMCs:
      highCMCs_failing <- originalMethodCMCs[[which(is.na(thetaMax))]] %>%
        dplyr::ungroup() %>%
        dplyr::mutate(direction = names(thetaMax)[which(is.na(thetaMax))]) %>%
        dplyr::mutate(originalMethodClassif = "CMC",
                      highCMCClassif = "CMC")

      #the direction that passed the high CMC criterion gets all of its high
      #CMCs
      highCMCs_passing <- cmcPerTheta[[which(!is.na(thetaMax))]] %>%
        dplyr::filter(.data$theta >= thetaMax[[which(!is.na(thetaMax))]] - thetaThresh & .data$theta <= thetaMax[[which(!is.na(thetaMax))]] + thetaThresh & .data$highCMCClassif == "CMC") %>%
        dplyr::mutate(direction = names(thetaMax)[[which(!is.na(thetaMax))]])

      highCMCs <- highCMCs_passing %>%
        dplyr::bind_rows() %>%
        dplyr::bind_rows(highCMCs_failing) %>%
        dplyr::distinct() %>%
        dplyr::group_by(.data$cellIndex) %>% #we don't want a cell being double-counted between the two comparisons
        dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
        dplyr::ungroup() %>%
        dplyr::select(-c(.data$originalMethodClassif,.data$highCMCClassif))

      #we want to make sure that the modal theta value in one direction is the
      #opposite (or close to the opposite) of the modal theta value in the other
      #direction
      thetaMax_dismissed <- highCMCs %>%
        dplyr::group_by(.data$direction,.data$theta) %>%
        dplyr::tally() %>%
        dplyr::filter(.data$n == max(.data$n))

      #it's theoretically possible, albeit improbable, that one direction will
      #pass the high CMC criterion while the other direction fails *and*
      #produces 0 initial CMCs. In this case, we would only take the high CMCs
      #in the one direction. Not including this if statement first would throw
      #an error in the next if statement
      if(nrow(thetaMax_dismissed) == 1){
        return(list("originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = highCMCs))
      }

      #another possibility is that more than one theta in one direction ties for
      #the CMC max count -- we can again determine whether these theta values
      #are "close" to each other (if not, then fail the criterion) and otherwise
      #take their median. Note that if there are three
      else if(nrow(thetaMax_dismissed) > 2){

        thetaMax_reference_v_target_diff <- thetaMax_dismissed %>%
          dplyr::filter(.data$direction =="reference_v_target") %>%
          dplyr::arrange(.data$theta) %>%
          dplyr::pull(.data$theta) %>%
          diff()

        thetaMax_target_v_reference_diff <- thetaMax_dismissed %>%
          dplyr::filter(.data$direction =="target_v_reference") %>%
          dplyr::arrange(.data$theta) %>%
          dplyr::pull(.data$theta) %>%
          diff()

        reference_v_target_failure <- any(thetaMax_reference_v_target_diff > thetaThresh)

        target_v_reference_failure <- any(thetaMax_target_v_reference_diff > thetaThresh)

        #in the event of a failure, don't assign any high CMCs

        if(reference_v_target_failure | target_v_reference_failure){
          return(list("originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = highCMC_initial))
        }
        #otherwise, take the median of the theta values within their respective
        #directions
        else{
          thetaMax_dismissed <- thetaMax_dismissed %>%
            dplyr::group_by(.data$direction) %>%
            dplyr::summarise(theta = median(.data$theta))
        }
      }
      #if thetaMax_dismissed has length 2, then we want to make sure that these
      #are actually opposites or close to opposites of each other. If not, then
      #we won't assign any high CMCs to the comparison
      thetaMax_dismissed <- thetaMax_dismissed %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$direction) #%>%
      # dplyr::pull(theta)

      if(nrow(thetaMax_dismissed) == 1){

        if(compareThetas){

          intitialCMCs_nonMissingDirection <- originalMethodCMCs[[which(names(thetaMax) == thetaMax_dismissed$direction)]]

          thetaCompareBool_dismissed <- abs(thetaMax_dismissed$theta - median(intitialCMCs_nonMissingDirection$theta,na.rm = TRUE)) > thetaThresh

          if(thetaCompareBool_dismissed | is.na(thetaCompareBool_dismissed) | purrr::is_empty(thetaCompareBool_dismissed)){
            return(list("originalMethodCMCs" = originalMethodCMCs,
                        "highCMCs" = highCMC_initial))

          }
          else{
            return(list("originalMethodCMCs" = originalMethodCMCs,
                        "highCMCs" = highCMCs))
          }
        }
        else{
          return(list("originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = highCMCs))
        }
      }

      if(compareThetas){
        thetaCompareBool_dismissed <- (((abs((thetaMax_dismissed[1,"theta"] - median(originalCMCs_reference_v_target$theta,na.rm = TRUE)))) > thetaThresh) |
                                         (abs((thetaMax_dismissed[2,"theta"] - median(originalCMCs_target_v_reference$theta,na.rm = TRUE))) > thetaThresh))

        if(thetaCompareBool_dismissed | is.na(thetaCompareBool_dismissed) | purrr::is_empty(thetaCompareBool_dismissed)){
          return(list("originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = highCMC_initial))
        }
      }

      sign(thetaMax$reference_v_target) == sign(thetaMax$target_v_reference) & abs(thetaMax$reference_v_target - -1*thetaMax$target_v_reference) > thetaThresh

      if((sign(thetaMax_dismissed[1,"theta"]) == sign(thetaMax_dismissed[2,"theta"]) & sign(thetaMax_dismissed[1,"theta"]) != 0 & sign(thetaMax_dismissed[2,"theta"]) != 0) |
         (abs((abs(thetaMax_dismissed[1,"theta"]) - abs(thetaMax_dismissed[2,"theta"]))) > thetaThresh)){
        return(list("originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = highCMC_initial))
      }
      # if we've made it this far, then the theta values should be at least to
      # within thetaThresh of being opposites of each other, so we can return
      # the highCMCs without worrying about a disagreement
      else{
        return(list("originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = highCMCs))
      }
    }

    # Most "conservative" decision is to flat-out fail the whole comparison if
    # one direction doesn't pass the high CMC criterion
    #else if(missingThetaDecision == "fail"){
    if(missingThetaDecision == "fail"){
      return(list("originalMethodCMCs" = originalMethodCMCs,
                  "highCMCs" = highCMC_initial))
    }
  }

  #if we've made it this far in the function, then both directions pass the high
  #CMC criterion. That is a CMC count "mode" has been identified at a particular
  #theta value in both directions.

  #It is often the case for known non-matches, even if they pass the high CMC
  #criterion, that the theta values obtained from the initially proposed method
  #disagree considerably with the theta values obtained from the high CMC
  #method. This leads to assigning a lot of false positive high CMCs when
  #there's already evidence that they are in fact NOT matching. We can use this
  #to our advantage for teasing out non-matches. If there is a disagreement in
  #theta values in *either* direction between the initial and High CMC methods,
  #then the pair will not be assigned any high CMCs. Note: this might be changed
  #to determine whether there is a disagreement in *both* directions later on if
  #this is determined to be too strict

  if(compareThetas){

    thetaCompareBool <- (((abs((thetaMax$reference_v_target - median(originalCMCs_reference_v_target$theta,na.rm = TRUE)))) > thetaThresh) |
                           (abs((thetaMax$target_v_reference - median(originalCMCs_target_v_reference$theta,na.rm = TRUE))) > thetaThresh))

    if(thetaCompareBool | is.na(thetaCompareBool) | purrr::is_empty(thetaCompareBool)){
      return(list("originalMethodCMCs" = originalMethodCMCs,
                  "highCMCs" = highCMC_initial))

    }
  }

  #The last contingency is making sure that these theta modes are opposites of
  #each other (or within thetaThresh of being opposites).

  #The following determines whether (1) the theta modes in either direction are
  #the same sign as each other (which they shouldn't be for matches) AND (2)
  #whether they are far away from 0. An older version of this function didn't
  #cover cases in which, for example, the thetaThresh was 6 and both
  #comparisons voted for a consensual theta value of 3. Given the thetaThresh,
  #this should be admissible since they are "close enough" (according to
  #thetaThresh) to each other.
  if(sign(thetaMax$reference_v_target) == sign(thetaMax$target_v_reference) & abs(thetaMax$reference_v_target - -1*thetaMax$target_v_reference) > thetaThresh){

    return(list("originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = highCMC_initial))
  }
  #Even if the two theta values are of opposite signs, it's possible that they
  #don't agree with each other. For example, one direction might vote for -27
  #degrees as the consensual theta value while the other votes for 12. Such a
  #comparison shouldn't pass the High CMC criterion.
  else if((abs((abs(thetaMax$reference_v_target) - abs(thetaMax$target_v_reference))) > thetaThresh)){
    return(list("originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = highCMC_initial))
  }
  else{
    highCMCs <- purrr::pmap(.l = list(cmcPerTheta,
                                      c("reference_v_target","target_v_reference"),
                                      thetaMax),
                            function(cmcs,compDirection,th){
                              cmcs %>%
                                dplyr::filter(.data$theta >= th - thetaThresh & .data$theta <= th + thetaThresh & .data$highCMCClassif == "CMC") %>%
                                dplyr::mutate(direction = compDirection)
                            }) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct() %>%
      dplyr::group_by(.data$cellIndex) %>% #we don't want a cell being double-counted between the two comparisons
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)), na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(.data$originalMethodClassif,.data$highCMCClassif))

    return(list("originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = highCMCs))
  }
}

#'Combine data frames containing CMC results from 2 comparison directions
#'
#'@name decision_combineDirections
#'
#'@description Combines CMC results from two comparison directions of a single
#'  cartridge case pair (i.e., where each cartridge case scan has been treated
#'  as both the reference and target scan). This function assumes that the CMC
#'  results are data frames withcolumns called "originalMethodClassif" and
#'  "highCMCClassif" containing CMCs identified under the original method of
#'  Song (2013) and the High CMC method of Tong et al. (2015) (see example).
#'@param reference_v_target_CMCs CMCs for the comparison between the reference
#'  scan and the target scan.
#'@param target_v_reference_CMCs (optional) CMCs for the comparison between the
#'  target scan and the reference scan. If this is missing, then only the
#'  original method CMCs will be plotted
#'@param corColName name of correlation similarity score column used to identify
#'  the CMCs in the two comparison_*_df data frames (e.g., pairwiseCompCor)
#'@param missingThetaDecision dictates how function should handle situations in
#'  which one direction passes the high CMC criterion while another direction
#'  does not. "dismiss": only counts the initial CMCs in failed direction and
#'  high CMCs in successful direction. "fail": only counts the initial CMCs in
#'  either direction and returns the minimum of these two numbers.
#'@param compareThetas dictates if the consensus theta values determined under
#'  the initially proposed method should be compared to the consensus theta
#'  values determined under the High CMC method. In particular, determines for
#'  each direction whether the consensus theta values determined under the two
#'  methods are within theta_thresh of each other. It is often the case that
#'  non-matching cartridge cases, even if they pass the High CMC criterion, will
#'  have differing consensus theta values under the two methods. If this isn't
#'  taken into account, non-matches tend to be assigned a lot of false positive
#'  CMCs under the High CMC method.
#'@param thetaThresh (original method of Song (2013)) used to classify
#'  particular theta values "congruent" if they are within thetaThresh of the
#'  median theta value. (High CMC) defines how wide a High CMC mode is allowed
#'  to be in the CMC-theta distribution before it's considered too diffuse. This
#'  is also used in this function to determine whether the estimated alignment
#'  theta values from the two comparison directions are "approximately" opposite
#'  (i.e., within thetaThresh of each other in absolute value), which they
#'  should be if the cartridge case pair is a known match.
#'
#'@return a list of 2 elements: (1) the CMCs identified under the original
#'  method of Song (2013) for both comparison directions since Song (2013) does
#'  not indicate whether/how results are combined and (2) the combined CMC
#'  results under the High CMC method.
#'@examples
#'\dontrun{
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'comparisonDF_1to2 <- purrr::map_dfr(seq(-30,30,by = 3),
#'                                    ~ comparison_allTogether(fadul1.1_processed,
#'                                                        fadul1.2_processed,
#'                                                        theta = .))
#'comparisonDF_2to1 <- purrr::map_dfr(seq(-30,30,by = 3),
#'                                    ~ comparison_allTogether(fadul1.2_processed,
#'                                                        fadul1.1_processed,
#'                                                        theta = .))
#'
#'comparisonDF_1to2 <- comparisonDF_1to2 %>%
#' dplyr::mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
#'                                                    x = x,
#'                                                    y = y,
#'                                                    theta = theta,
#'                                                    corr = pairwiseCompCor),
#'               highCMCClassif = decision_CMC(cellIndex = cellIndex,
#'                                            x = x,
#'                                            y = y,
#'                                            theta = theta,
#'                                            corr = pairwiseCompCor,
#'                                            tau = 1))
#'
#'
#'comparisonDF_2to1 <- comparisonDF_2to1 %>%
#' dplyr::mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
#'                                                    x = x,
#'                                                    y = y,
#'                                                    theta = theta,
#'                                                    corr = pairwiseCompCor),
#'               highCMCClassif = decision_CMC(cellIndex = cellIndex,
#'                                            x = x,
#'                                            y = y,
#'                                            theta = theta,
#'                                            corr = pairwiseCompCor,
#'                                            tau = 1))
#'
#'decision_combineDirections(comparisonDF_1to2,comparisonDF_2to1)
#'}
#'@export

decision_combineDirections <- function(reference_v_target_CMCs,
                                       target_v_reference_CMCs,
                                       corColName = "pairwiseCompCor",
                                       missingThetaDecision = "fail",
                                       compareThetas = TRUE,
                                       thetaThresh = 6){

  #internally, this function will basically just call a legacy function (for a
  #version of the package that wasn't pipe-friendly) that implements the
  #necessary logic to combine the data from two comparison directions.

  cmcResults <- cmcFilter_improved(reference_v_target_CMCs = reference_v_target_CMCs,
                                   target_v_reference_CMCs = target_v_reference_CMCs,
                                   corColName = corColName,
                                   missingThetaDecision = missingThetaDecision,
                                   compareThetas = compareThetas,
                                   thetaThresh = thetaThresh)

  return(cmcResults)
}

# Calculates the mode of a vector of numbers
#
# @name getMode
#
# @description Calculates the mode of a vector. Can be used as a consensus
#   function in cmcR::cmcFilter or cmcR::cmcFilter_improved.
#
# @param x a numeric vector
#
# @examples
# \dontrun{
# #x3p1 and x3p2 are two x3p objects containing processed surface matrices
#
# #calculate "initial" cmcs
# comparison1$ccfResults %>%
#   cmcR::topResultsPerCell() %>%
#   cmcR::cmcFilter(consensus_function = median,
#                   ccf_thresh = .4,
#                   dx_thresh = 20,
#                   dy_thresh = dx_thresh,
#                   theta_thresh = 3,
#                   consensus_function_theta = getMode)
# }
#@export

getMode <- function(x){
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
