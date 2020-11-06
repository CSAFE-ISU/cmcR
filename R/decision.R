#' Apply CMC classification logic of the original method of Song (2013)
#'
#' @name decision_originalMethod_classifyCMCs
#' @param cellIndex vector/tibble column containing cell indices corresponding
#'   to a reference cell
#' @param x vector/tibble column containing x horizontal translation values
#' @param y vector/tibble column containing y vertical translation values
#' @param theta vector/tibble column containing theta rotation values
#' @param corr vector/tibble column containing correlation similarity scores
#'   between a reference cell and its associated target region
#' @param xThresh used to classify particular x values "congruent" if they are
#'   within xThresh of the median x value
#' @param yThresh used to classify particular y values "congruent" if they are
#'   within yThresh of the median y value
#' @param thetaThresh used to classify particular theta values "congruent" if
#'   they are within thetaThresh of the median theta value
#' @param corrThresh to classify particular correlation values "congruent" if
#'   they are at least corrThresh
#'
#' @seealso \url{https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193}
#'@note the decision_CMC function internally calls this function if a value of
#'   tau is not provided
#'
#' @keywords internal

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

#' Compute CMC-theta distribution for a set of comparison features
#'
#' @name decision_highCMC_cmcThetaDistrib
#' @param cellIndex vector/tibble column containing cell indices corresponding
#'   to a reference cell
#' @param x vector/tibble column containing x horizontal translation values
#' @param y vector/tibble column containing y vertical translation values
#' @param theta vector/tibble column containing theta rotation values
#' @param corr vector/tibble column containing correlation similarity scores
#'   between a reference cell and its associated target region
#' @param xThresh used to classify particular x values "congruent" (conditional
#'   on a particular theta value) if they are within xThresh of the
#'   theta-specific median x value
#' @param yThresh used to classify particular y values "congruent" (conditional
#'   on a particular theta value) if they are within yThresh of the
#'   theta-specific median y value
#' @param corrThresh to classify particular correlation values "congruent"
#'   (conditional on a particular theta value) if they are at least corrThresh
#'
#' @note This function is a helper internally called in the decision_CMC
#'   function. It is exported to be used as a diagnostic tool for the High CMC
#'   method
#'
#'@examples
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
#' @importFrom rlang .data
#' @export

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

#' Classify theta values in CMC-theta distribution as having "High" or "Low" CMC
#' candidate counts
#'
#' @name decision_highCMC_identifyHighCMCThetas
#'
#' @param cmcThetaDistrib output of the decision_highCMC_cmcThetaDistrib
#'   function
#' @param tau constant used to define a "high" CMC count. This number is
#'   subtracted from the maximum CMC count achieved in the CMC-theta
#'   distribution. Theta values with CMC counts above this value are considered
#'   to have "high" CMC counts.
#'
#' @note This function is a helper internally called in the decision_CMC
#'   function. It is exported to be used as a diagnostic tool for the High CMC
#'   method
#'@examples
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
#' @importFrom rlang .data
#' @export

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

#'Apply CMC classification logic of the Tong et al. (2015) to the CMC-theta
#'distribution returned by the decision_highCMC_cmcThetaDistrib function
#'
#'@name decision_highCMC_classifyCMCs
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
#'@param thetaThresh defines how wide a High CMC mode is allowed to be in the
#'  CMC-theta distribution before it's considered too diffuse
#'@param corrThresh to classify particular correlation values "congruent"
#'  (conditional on a particular theta value) if they are at least corrThresh
#'@param tau constant used to define a "high" CMC count. This number is
#'  subtracted from the maximum CMC count achieved in the CMC-theta
#'  distribution. Theta values with CMC counts above this value are considered
#'  to have "high" CMC counts.
#'@seealso
#'\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'
#'@note the decision_CMC function internally calls this function if a value of
#'  tau is provided
#'@importFrom rlang .data
#' @keywords internal

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
#'@seealso \url{https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193}
#'@seealso
#'\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'@examples
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

#' Calculates the mode of a vector of numbers
#'
#' @name getMode
#'
#' @description Calculates the mode of a vector. Can be used as a consensus
#'   function in cmcR::cmcFilter or cmcR::cmcFilter_improved.
#'
#' @param x a numeric vector
#'
#' @examples
#' \dontrun{
#' #x3p1 and x3p2 are two x3p objects containing processed surface matrices
#' comparison1 <- cellCCF(x3p1,x3p2) #defaults assumed
#'
#' #calculate "initial" cmcs
#' comparison1$ccfResults %>%
#'   cmcR::topResultsPerCell() %>%
#'   cmcR::cmcFilter(consensus_function = median,
#'                   ccf_thresh = .4,
#'                   dx_thresh = 20,
#'                   dy_thresh = dx_thresh,
#'                   theta_thresh = 3,
#'                   consensus_function_theta = getMode)
#' }
#'
#' @export

getMode <- function(x){
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
