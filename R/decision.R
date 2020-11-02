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
#' dplyr::mutate(originalMethodClassif = decision_originalMethod_classifyCMCs(cellIndex = cellIndex,
#'                                                                            x = x,
#'                                                                            y = y,
#'                                                                            theta = theta,
#'                                                                            corr = pairwiseCompCor))
#'
#'comparisonDF %>%
#' dplyr::filter(originalMethodClassif == "CMC")
#' @export

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
#' @export

decision_highCMC_identifyHighCMCThetas <- function(cmcThetaDistrib,
                                                   tau = 1){

  thetaClassifications <- cmcThetaDistrib %>%
    dplyr::filter(cmcThetaDistribClassif == "CMC Candidate") %>%
    dplyr::group_by(theta) %>%
    dplyr::tally(name = "cmcCandidateCount") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(thetaCMCIdentif = ifelse(cmcCandidateCount >= max(cmcCandidateCount) - tau,"High","Low"))

  cmcThetaDistrib %>%
    dplyr::left_join(thetaClassifications,by = "theta")
}

#' Apply CMC classification logic of the Tong et al. (2015) to the CMC-theta
#' distribution returned by the decision_highCMC_cmcThetaDistrib function
#'
#' @name decision_highCMC_classifyCMCs
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
#' @param tau constant used to define a "high" CMC count. This number is
#'   subtracted from the maximum CMC count achieved in the CMC-theta
#'   distribution. Theta values with CMC counts above this value are considered
#'   to have "high" CMC counts.
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'
#' @note the decision_CMC function internally calls this function if a value of
#'   tau is provided
#'@examples
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'comparisonDF <- purrr::map_dfr(seq(-30,30,by = 3),
#'                               ~ comparison_allTogether(fadul1.1_processed,
#'                                                        fadul1.2_processed,
#'                                                        theta = .))
#'
#'comparisonDF <- comparisonDF %>%
#' dplyr::mutate(highCMCClassif = decision_highCMC_classifyCMCs(cellIndex = cellIndex,
#'                                                              x = x,
#'                                                              y = y,
#'                                                              theta = theta,
#'                                                              corr = pairwiseCompCor))
#'
#'comparisonDF %>%
#' dplyr::filter(highCMCClassif == "CMC")
#' @export

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
    dplyr::filter(cmcThetaDistribClassif == "CMC Candidate")

  cmcThetaDistrib_classified <- decision_highCMC_identifyHighCMCThetas(cmcThetaDistrib,
                                                                       tau = tau)

  passesHighCMCCriterion <- cmcThetaDistrib_classified %>%
    dplyr::filter(thetaCMCIdentif == "High") %>%
    dplyr::select(theta) %>%
    dplyr::distinct() %>%
    dplyr::summarise(distance = abs(max(theta) - min(theta))) %>%
    dplyr::pull(distance) %>%
    {. <= thetaThresh}

  if(passesHighCMCCriterion){
    highCMCs <- cmcThetaDistrib_classified %>%
      dplyr::mutate(highCMCClassif = ifelse(thetaCMCIdentif == "High","CMC","non-CMC (passed)")) %>%
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

#' Combine CMC data from two comparison directions
#'
#' @name decision_combineCMCDirections
#'
#' @param comparison_1to2_df a data frame containing a column
#'   "originalMethodClassif" or "highCMCClassif" (or both) representing CMCs
#'   identified for a comparison between two cartridge cases in which one is
#'   considered the "reference scan" and the other the "target scan". The
#'   decision_CMC function can be used to add CMC column(s) to a data frame.
#' @param comparison_2to1_df a data frame containing a column
#'   "originalMethodClassif" or "highCMCClassif" (or both) representing CMCs
#'   identified for a comparison between two cartridge cases where the roles of
#'   "reference" and "target" are reversed relative to the comparison_1to2_df
#'   CMCs
#' @param corColName name of correlation similarity score column used to
#'   identify the CMCs in the two comparison_*_df data frames (e.g.,
#'   pairwiseCompCor)
#' @param method whether to apply the combination logic of the High CMC method
#'   or the original method of Song (2013)
#' @param compareThetas controls whether the estimated "true" rotation values
#'   calculated under the original method of Song (2013) and the High CMC method
#'   are compared to each other (they should be similar for a truly matching
#'   scan).
#' @param thetaThresh threshold to use when comparing the estimated theta values
#'   under the original method of Song (2013) and the High CMC method. Also used
#'   to compare the estimated theta values in both comparison directions under
#'   the High CMC method (as they should be approximate opposites of each
#'   other).
#'@examples
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
#'decision_combineCMCDirections(comparisonDF_1to2,
#'                              comparisonDF_2to1,
#'                              corColName = "pairwiseCompCor",
#'                              method = "highCMC")
#' @export

decision_combineCMCDirections <- function(comparison_1to2_df,
                                          comparison_2to1_df,
                                          corColName = "pairwiseCompCor",
                                          method = "highCMC",
                                          compareThetas = TRUE,
                                          thetaThresh = 6){
  #missingTheta_decision is only applicable when method == "highCMC" or "highCMC".
  #compareThetas argument is only applicable when method == "highCMC". Test for
  #this.

  #For the original method, it's unclear which of the two directions should be
  #taken as the "official" CMC count. To be conservative, assign the pair the
  #minimum of the two CMC counts.
  if(method == "original"){

    if(any(stringr::str_detect(names(comparison_1to2_df),"high")) |
       any(stringr::str_detect(names(comparison_2to1_df),"high"))){
      stop("Cannot have High CMC method information for method = 'original'")
    }

    originalMethodCMCs <- dplyr::bind_rows(comparison_1to2_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_AtoB"),
                                           comparison_2to1_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_BtoA"))

    if(nrow(originalMethodCMCs) == 0){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }

    originalMethodCMCs <- originalMethodCMCs %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(direction) %>%
      dplyr::group_split()

    if(purrr::is_empty(originalMethodCMCs)){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }
    else{
      minIndex <- originalMethodCMCs %>%
        purrr::map_int(nrow) %>%
        which.min()

      highCMCs<- originalMethodCMCs[[minIndex]]
    }

    return(originalMethodCMCs)
  }

  #If both directions fail the High CMC criterion, assign them the minimum of
  #the two original method CMC counts (which may be 0)
  if(any(stringr::str_detect(comparison_1to2_df$highCMCClassif,"failed")) |
     any(stringr::str_detect(comparison_2to1_df$highCMCClassif,"failed"))){

    originalMethodCMCs <- dplyr::bind_rows(comparison_1to2_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_AtoB"),
                                           comparison_2to1_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_BtoA"))

    if(nrow(originalMethodCMCs) == 0){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }

    originalMethodCMCs <- originalMethodCMCs %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(direction) %>%
      dplyr::group_split()

    if(purrr::is_empty(originalMethodCMCs)){
      highCMCs <- data.frame(cellIndex = vector(length = 0),
                             x = vector(length = 0),
                             y = vector(length = 0),
                             fft_ccf = vector(length = 0),
                             pairwiseCompCor = vector(length = 0),
                             theta = vector(length = 0),
                             highCMCCriterion = vector(length = 0),
                             direction = vector(length = 0),
                             originalMethodClassif = vector(length = 0),
                             highCMCClassif = vector(length = 0))
    }
    else{
      minIndex <- originalMethodCMCs %>%
        purrr::map_int(nrow) %>%
        which.min()

      highCMCs<- originalMethodCMCs[[minIndex]]
    }

    return(highCMCs)

  }

  #Conservative: if one direction fails the High CMC criterion, then behave as
  #if both directions failed. Might implement a "moderate" or "liberal" option
  #(where we, e.g., replace the missing theta value in one direction and carry
  #on as if neither direction failed), but the conservative option currently
  #seems to yield best results.

  #Checking if either direction failed
  if(any(stringr::str_detect(comparison_1to2_df$highCMCClassif,"failed")) |
     any(stringr::str_detect(comparison_2to1_df$highCMCClassif,"failed"))){

    originalMethodCMCs <- dplyr::bind_rows(comparison_1to2_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_AtoB"),
                                           comparison_2to1_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_BtoA"))

    if(nrow(originalMethodCMCs) == 0){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }

    originalMethodCMCs <- originalMethodCMCs %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(direction) %>%
      dplyr::group_split()

    if(purrr::is_empty(originalMethodCMCs)){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }
    else{
      minIndex <- originalMethodCMCs %>%
        purrr::map_int(nrow) %>%
        which.min()

      highCMCs<- originalMethodCMCs[[minIndex]]
    }

    return(highCMCs)

  }

  #If we've gotten this far, then both directions pass the High CMC criterion.
  #However, non-matches may pass the High CMC criterion in both directions yet
  #vote for wildly different theta values in either direction (e.g., -30 degrees
  #in one direction vs. 3 degrees in the other). We need to make sure that the
  #theta values are approximately each other's opposites. Note that if multiple
  #theta values tie for the mode, we take their median.

  thetaMode_1to2 <- comparison_1to2_df %>%
    dplyr::filter(highCMCClassif == "CMC") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(theta) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::pull(theta) %>%
    median()

  thetaMode_2to1 <- comparison_2to1_df %>%
    dplyr::filter(highCMCClassif == "CMC") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(theta) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::pull(theta) %>%
    median()

  #The first boolean expression determines whether (1) the theta modes in either
  #direction are the same sign as each other (which they shouldn't be for
  #matches) AND (2) whether they are far away from 0 (e.g., if both directions
  #vote for a 0 degree rotation, then we shouldn't be failing the pair). The
  #expression after the "|" checks if they are far from each other in
  #absolute value (e.g., if one direction votes for -30 degrees while the other
  #for 12 degrees)
  if((sign(thetaMode_1to2) == sign(thetaMode_2to1) & abs(thetaMode_1to2 - -1*thetaMode_2to1) > thetaThresh) |
     (abs(abs(thetaMode_1to2) - abs(thetaMode_2to1)) > thetaThresh)){

    originalMethodCMCs <- dplyr::bind_rows(comparison_1to2_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_AtoB"),
                                           comparison_2to1_df %>%
                                             dplyr::filter(originalMethodClassif == "CMC") %>%
                                             dplyr::mutate(direction = "comparison_BtoA"))

    if(nrow(originalMethodCMCs) == 0){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }

    originalMethodCMCs <- originalMethodCMCs %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(direction) %>%
      dplyr::group_split()

    if(purrr::is_empty(originalMethodCMCs)){
      return(data.frame(cellIndex = vector(length = 0),
                        x = vector(length = 0),
                        y = vector(length = 0),
                        fft_ccf = vector(length = 0),
                        pairwiseCompCor = vector(length = 0),
                        theta = vector(length = 0),
                        highCMCCriterion = vector(length = 0),
                        direction = vector(length = 0),
                        originalMethodClassif = vector(length = 0),
                        highCMCClassif = vector(length = 0)))
    }
    else{
      minIndex <- originalMethodCMCs %>%
        purrr::map_int(nrow) %>%
        which.min()

      highCMCs<- originalMethodCMCs[[minIndex]]
    }

    return(highCMCs)

  }

  #Optionally, we can compare the median theta values from the original method
  #and the modal theta values from the High CMC method. Again, for a matching
  #pair, we would expect agreement no matter what method we use.

  if(compareThetas){
    medianTheta_1to2 <- comparison_1to2_df %>%
      dplyr::filter(originalMethodClassif == "CMC") %>%
      dplyr::ungroup() %>%
      dplyr::pull(theta) %>%
      median()

    medianTheta_2to1 <- comparison_2to1_df %>%
      dplyr::filter(highCMCClassif == "CMC") %>%
      dplyr::ungroup() %>%
      dplyr::pull(theta) %>%
      median()

    #compare the modal and median theta values estimated under the two methods
    #in both directions. If either of the directions are too large (larger than
    #thetaThresh), then say that pair fails the High CMC criterion.
    if(abs(thetaMode_1to2 - medianTheta_1to2) > thetaThresh | abs(thetaMode_2to1 - medianTheta_2to1) > thetaThresh){

      originalMethodCMCs <- dplyr::bind_rows(comparison_1to2_df %>%
                                               dplyr::filter(originalMethodClassif == "CMC") %>%
                                               dplyr::mutate(direction = "comparison_AtoB"),
                                             comparison_2to1_df %>%
                                               dplyr::filter(originalMethodClassif == "CMC") %>%
                                               dplyr::mutate(direction = "comparison_BtoA"))

      if(nrow(originalMethodCMCs) == 0){
        return(data.frame(cellIndex = vector(length = 0),
                          x = vector(length = 0),
                          y = vector(length = 0),
                          fft_ccf = vector(length = 0),
                          pairwiseCompCor = vector(length = 0),
                          theta = vector(length = 0),
                          highCMCCriterion = vector(length = 0),
                          direction = vector(length = 0),
                          originalMethodClassif = vector(length = 0),
                          highCMCClassif = vector(length = 0)))
      }

      originalMethodCMCs <- originalMethodCMCs %>%
        dplyr::group_by(cellIndex) %>%
        dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(direction) %>%
        dplyr::group_split()

      if(purrr::is_empty(originalMethodCMCs)){
        return(data.frame(cellIndex = vector(length = 0),
                          x = vector(length = 0),
                          y = vector(length = 0),
                          fft_ccf = vector(length = 0),
                          pairwiseCompCor = vector(length = 0),
                          theta = vector(length = 0),
                          highCMCCriterion = vector(length = 0),
                          direction = vector(length = 0),
                          originalMethodClassif = vector(length = 0),
                          highCMCClassif = vector(length = 0)))
      }
      else{
        minIndex <- originalMethodCMCs %>%
          purrr::map_int(nrow) %>%
          which.min()

        highCMCs <- originalMethodCMCs[[minIndex]]
      }

      return(highCMCs)
    }
  }

  #If all contingencies pass, assign all of the appropriate CMCs, excluding
  #replicates:
  highCMCs <- bind_rows(comparison_1to2_df %>%
                          dplyr::filter(originalMethodClassif == "CMC" | highCMCClassif == "CMC") %>%
                          dplyr::mutate(direction = "comparison_AtoB"),
                        comparison_2to1_df %>%
                          dplyr::filter(originalMethodClassif == "CMC" | highCMCClassif == "CMC") %>%
                          dplyr::mutate(direction = "comparison_BtoA")) %>%
    tidyr::pivot_longer(cols = c(originalMethodClassif,highCMCClassif),
                        names_to = "method",
                        values_to = "classification") %>%
    dplyr::group_by(cellIndex,method) %>%
    filter(classification == "CMC") %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName)))) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = 1:(ncol(.)-2),
                       names_from = method,
                       values_from = classification) %>%
    dplyr::mutate(highCMCClassif = ifelse(is.na(highCMCClassif),"non-CMC","CMC"),
                  originalMethodClassif = ifelse(is.na(originalMethodClassif),"non-CMC","CMC"),
                  highCMCCriterion = "passed")

  return(highCMCs)
}