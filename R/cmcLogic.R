#' Apply CMC classification logic of the original method of Song (2013)
#'
#' @name decision_originalMethod
#'
#' @seealso \url{https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193}
#'
#' @export

decision_originalMethod_classifyCMCs <- function(comparisonFeaturesDF,
                                                 corColName = "pairwiseCompCor",
                                                 xThresh = 20,
                                                 yThresh = xThresh,
                                                 thetaThresh = 6,
                                                 corThresh = .5){

  comparisonFeaturesDF %>%
    dplyr::group_by(cellIndex) %>%
    dplyr::top_n(n = 1,wt = (!!as.name(corColName))) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(originalMethodClassif = ifelse(abs(x - median(x)) <= xThresh &
                                                   abs(y - median(y)) <= yThresh &
                                                   abs(theta - median(theta)) <= thetaThresh &
                                                   (!!as.name(corColName)) >= corThresh,"CMC","non-CMC")) %>%
    dplyr::select(cellIndex,x,y,(!!as.name(corColName)),originalMethodClassif)
}

#' Compute CMC-theta distribution for a set of comparison features
#'
#' @name decision_highCMC_cmcThetaDistrib
#'
#' @export

decision_highCMC_cmcThetaDistrib <- function(comparisonFeaturesDF,
                                             corColName = "pairwiseCompCor",
                                             xThresh = 20,
                                             yThresh = xThresh,
                                             corThresh = .5){

  comparisonFeaturesDF %>%
    dplyr::group_by(theta) %>%
    dplyr::filter(abs(x - median(x)) <= xThresh &
                    abs(y - median(y)) <= yThresh &
                    (!!as.name(corColName)) >= corThresh) %>%
    dplyr::ungroup() %>%
    dplyr::select(cellIndex,x,y,(!!as.name(corColName)),theta)

}

#' Classify theta values in CMC-theta distribution as having "High" or "Low" CMC
#' candidate counts
#'
#' @name decision_highCMC_identifyHighCMCThetas
#'
#' @export

decision_highCMC_identifyHighCMCThetas <- function(cmcThetaDistrib,
                                                   tau = 1){

  thetaClassifications <- cmcThetaDistrib %>%
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
#'
#' @seealso
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'
#'
#' @export

decision_highCMC_classifyCMCs <- function(cmcThetaDistrib,
                                          tau = 1,
                                          thetaThresh = 6){

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
      dplyr::mutate(highCMCClassif = ifelse(thetaCMCIdentif == "High","CMC","non-CMC")) %>%
      dplyr::select(-c(thetaCMCIdentif,cmcCandidateCount))
  }
  else{
    highCMCs <- cmcThetaDistrib_classified %>%
      dplyr::mutate(highCMCClassif = "non-CMC") %>%
      dplyr::select(-c(thetaCMCIdentif,cmcCandidateCount))
  }

  return(highCMCs)
}

#' Calculate x, y, theta values at which top similarity is achieved in a cell
#' pair.
#'
#' @name topResultsPerCell
#'
#' @description Given a list of data frames like the one returned in the
#'   `ccfResults` element of the list returned from cmcR::cellCCF, returns a
#'   data frame of the theta, x, y, ccf, and ccf values at which each cell
#'   pair attains maximum ccf.
#'
#' @param ccfResults list of data frames, like the one returned by
#'   cmcR::cellCCF, containing breech face cell comparison results
#'
#' @examples
#' \dontrun{
#' #x3p1 and x3p2 are two x3p objects containing processed surface matrices
#' comparison1 <- cellCCF(x3p1,x3p2) #defaults assumed
#'
#' comparison1$ccfResults %>%
#'   cmcR::topResultsPerCell()
#' }
#'
#' @seealso cmcR::cellCCF
#' @export
#'

utils::globalVariables(c(".","cellRange","ccf"))

topResultsPerCell <- function(ccfResults){

  ccfResults  %>%
    purrr::map2_dfr(.x = .,
                    .y = names(.),
                    ~ .x %>%
                      dplyr::mutate(theta = as.numeric(rep(.y,times = nrow(.))))) %>%
    dplyr::group_by(cellRange) %>%
    dplyr::filter(!is.na(ccf)) %>%
    dplyr::filter(ccf == max(ccf,na.rm = TRUE)) %>%
    dplyr::arrange(cellRange)
}

#' Implements CMC logic on a single data frame as proposed by Song (2013).
#'
#' @name cmcFilter
#'
#' @description Applies initially proposed Congruent Matching Cells method logic
#'   to the CCF results of a comparison between two cartridge case scans.
#'
#' @param ccfDF data frame containing ccf results from a comparison between two
#'   cartridge case scans
#' @param consensus_function function to aggregate the translation (x and y)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus x value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus y value that a cell
#'   pair can be to be called "congruent matching"
#' @param theta_thresh maximum distance from the consensus theta value that a
#'   cell pair can be to be called "congruent matching"
#' @param consensus_function_theta *(OPTIONAL)* function (separate from
#'   consensus_function) to aggregate the rotation (theta) values in the ccfDF
#'   data frame to determine "consensus" values
#' @param ... arguments to be passed to consensus_function or
#'   consensus_function_theta (e.g., na.rm = TRUE) if necessary
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
#'                   consensus_function_theta = median)
#' }
#'
#' @seealso
#' \url{https://pdfs.semanticscholar.org/4bf3/0b3a23c38d8396fa5e0d116cba63a3681494.pdf}
#'
#' @export
#'
#' @importFrom stats median

utils::globalVariables(c("x","y","theta","ccf"))

cmcFilter <- function(ccfDF,
                      consensus_function = median,
                      ccf_thresh = .6,
                      dx_thresh = 10,
                      dy_thresh = dx_thresh,
                      theta_thresh = 3,
                      consensus_function_theta = consensus_function,...){
  #Required tests:
  # ccfResults needs to be a dataframe containing columns: ccf,theta,x,y
  # ccf_thresh should be between 0 and 1
  # dx_thresh should positive
  # consensus_function and consensus_function_theta should be a function name that exists

  consensus_dx <- consensus_function(ccfDF$x,...)
  consensus_dy <- consensus_function(ccfDF$y,...)
  consensus_theta <- consensus_function_theta(ccfDF$theta,...)

  ccfDF %>%
    dplyr::filter(ccf >= ccf_thresh &
                    x >= consensus_dx - dx_thresh & x <= consensus_dx + dx_thresh &
                    y >= consensus_dy - dy_thresh & y <= consensus_dy + dy_thresh &
                    theta >= consensus_theta - theta_thresh & theta <= consensus_theta + theta_thresh)
}

#' @name cmcFilterPerTheta
#'
#' @param ccfResults list of data frames, like the one returned by
#'   cmcR::cellCCF, containing breech face cell comparison results
#' @param consensus_function function to aggregate the translation (x and y)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus x value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus y value that a cell
#'   pair can be to be called "congruent matching"
#' @param theta_thresh maximum distance from the consensus theta value that a
#'   cell pair can be to be called "congruent matching"
#' @param consensus_function_theta *(OPTIONAL)* function (separate from
#'   consensus_function) to aggregate the rotation (theta) values in the ccfDF
#'   data frame to determine "consensus" values
#' @param ... arguments to be passed to consensus_function or
#'   consensus_function_theta (e.g., na.rm = TRUE) if necessary
#'
#' @keywords internal
#'
#' @importFrom stats median setNames

utils::globalVariables(c(".","ccf"))

cmcFilterPerTheta <- function(ccfResults,
                              consensus_function = median,
                              ccf_thresh = .6,
                              dx_thresh = 10,
                              dy_thresh = dx_thresh,
                              theta_thresh = 3,
                              consensus_function_theta = consensus_function,...){

  ccfResults <- ccfResults %>%
    purrr::map2(.x = .,
                .y = names(.),
                function(thetaSpecificResults,theta){
                  thetaSpecificResults %>%
                    dplyr::select(-theta) %>%
                    dplyr::mutate(theta = rep(as.numeric(theta),
                                              times = nrow(.)))
                }) %>%
    setNames(names(ccfResults))

  ccfResults %>%
    purrr::map(~ cmcFilter(ccfDF = .,
                           consensus_function = consensus_function,
                           ccf_thresh = ccf_thresh,
                           dx_thresh = dx_thresh,
                           dy_thresh = dy_thresh,
                           theta_thresh = theta_thresh,
                           consensus_function_theta = median)) %>%
    dplyr::bind_rows()
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

#' @name calcMaxCMCTheta
#'
#' @param distanceToCMCMaxTieBreaker decides what to do in cases where there are
#'   consecutive theta values all tied for the max CMC count (no guidance is
#'   given by Tong et al. (2015) of what to do in this situation). The default
#'   is to determine the distance between any high CMC theta value and its
#'   closest CMC max theta value.
#'
#'   For example, suppose 3 consecutive theta values each have the max CMC count
#'   of 17. So the CMC counts around these theta values may look like ...10, 16,
#'   17, 17, 17, 15, 12... For a highCMC_thresh = 1, note that there is one
#'   theta value with associated count equal to the high CMC count of 17 - 1 =
#'   16. The default setting of distanceToCMCMaxTieBreaker = "minDistance" will
#'   consider the minimum distance between a max CMC theta value and this high
#'   CMC theta. So in this case, because a theta value with CMC count 17 is
#'   right next to this high CMC theta value, we would say that this particular
#'   cartridge case pair "passes" the high CMC criterion, assuming the
#'   theta_thresh is set to be equal to the grid spacing of the theta values
#'   considered (3 degrees by default).
#'
#'   A slightly more "conservative" option would be to set
#'   distanceToCMCMaxTieBreaker = "medDistance" take the median of these max CMC
#'   theta values and apply the high CMC criterion. The example under
#'   consideration *wouldn't* pass the high CMC criterion due to the particular
#'   choice of theta grid spacing, theta_thresh, and highCMC_thresh. However,
#'   different combinations of these arguments may lead to passing the
#'   criterion.
#'
#'   Finally, the most "conservative" option would be to set
#'   distanceToCMCMaxTieBreaker = "failCriterion" in which we immediately "fail"
#'   a cartridge case pair if there are any ties for max CMC theta (i.e., "there
#'   can only be one!" - Highlander).
#'
#' @keywords internal
#'
#' @importFrom stats median

utils::globalVariables(c("theta","n","distanceToCMCMax"))

calcMaxCMCTheta <- function(cmcPerTheta,
                            highCMC_thresh = 1,
                            theta_thresh = 3,
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
  #theta_thresh of each other, then we need to somehow "break" the ties. The
  #distanceToCMCMaxTieBreaker argument (see documentation) determines how to do
  #so.

  if(distanceToCMCMaxTieBreaker == "failCriterion" & nrow(cmcMax) > 1){
    return(NA)
  }

  #there may be more than one theta tied with the maximum CMC - in which case
  #we should determine if such thetas are "close" to each other, where
  #proximity is defined based on the theta_tresh set. If not, then we can rule
  #out that the comparison under scrutiny is a known match.
  if(any(diff(cmcMax$theta) > theta_thresh)){
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

  if(all(maxDistancetoCMCMax > theta_thresh)){
    return(NA)
  }
  else{
    if(distanceToCMCMaxTieBreaker == "medDistance"){
      return(median(cmcMax$theta))
    }
    return(median(cmcMax$theta))
  }
}

#' Implements "improved" CMC logic on a list of CCF results for a comparison
#' between two cartridge case scans as proposed by Tong et al. (2015)
#'
#' @name cmcFilter_improved
#'
#' @description Implements "improved' Congruent Matching Cells logic, as
#'   proposed by Tong et al. (2015), to the CCF results of a comparison between
#'   two cartridge case scans.
#'
#' @param cellCCF_bothDirections_output list returned by the function
#'   cmcR::cellCCF_bothdirections
#' @param consensus_function function to aggregate the translation (x and y)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus x value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus y value that a cell
#'   pair can be to be called "congruent matching"
#' @param theta_thresh maximum distance from the consensus theta value that a
#'   cell pair can be to be called "congruent matching"
#' @param missingTheta_decision dictates how function should handle situations
#'   in which one direction passes the high CMC criterion while another
#'   direction does not. "replace": replaces theta value in failed direction
#'   with opposite of theta value in successful direction. "dismiss": only
#'   counts the initial CMCs in failed direction and high CMCs in successful
#'   direction. "fail": only counts the initial CMCs in either direction.
#' @param compareInitialAndHighThetas dictates if the consensus theta values
#'   determined under the initially proposed method should be compared to the
#'   consensus theta values determined under the High CMC method. In particular,
#'   determines for each direction whether the consensus theta values determined
#'   under the two methods are within theta_thresh of each other. It is often
#'   the case that non-matching cartridge cases, even if they pass the High CMC
#'   criterion, will have differing consensus theta values under the two
#'   methods. If this isn't taken into account, non-matches tend to be assigned
#'   a lot of false positive CMCs under the High CMC method.
#' @param consensus_function_theta *(OPTIONAL)* function (separate from
#'   consensus_function) to aggregate the rotation (theta) values in the ccfDF
#'   data frame to determine "consensus" values
#'
#' @examples
#' \dontrun{
#' comparison1 <- cellCCF_bothDirections(x3p1,x3p2)
#'
#' cmcs <- cmcFilter_improved(comparison1,
#'                            consensus_function = median,
#'                            ccf_thresh = .4,
#'                            dx_thresh = 20,
#'                            dy_thresh = dx_thresh,
#'                            theta_thresh = 3,
#'                            missingTheta_decision = "replace",
#'                            consensus_function_theta = getMode)
#' }
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'
#' @export
#'
#' @importFrom stats median

utils::globalVariables(c(".","cellNum","comparison","theta"))

cmcFilter_improved <- function(cellCCF_bothDirections_output,
                               consensus_function = median,
                               ccf_thresh = .6,
                               dx_thresh = 10,
                               dy_thresh = dx_thresh,
                               theta_thresh = 3,
                               missingTheta_decision = "replace",
                               compareInitialAndHighThetas = FALSE,
                               consensus_function_theta = consensus_function){

  originalMethodCMCs <- cellCCF_bothDirections_output %>%
    purrr::map(~ topResultsPerCell(.$ccfResults) %>%
                 cmcFilter(consensus_function = consensus_function,
                           ccf_thresh = ccf_thresh,
                           dx_thresh = dx_thresh,
                           dy_thresh = dy_thresh,
                           theta_thresh = theta_thresh,
                           consensus_function_theta = consensus_function_theta) %>%
                 dplyr::ungroup())

  cmcPerTheta <-  cellCCF_bothDirections_output %>%
    purrr::map(~ cmcFilterPerTheta(ccfResults = .$ccfResults,
                                   consensus_function = consensus_function,
                                   ccf_thresh = ccf_thresh,
                                   dx_thresh = dx_thresh,
                                   dy_thresh = dy_thresh,
                                   theta_thresh = theta_thresh,
                                   consensus_function_theta = consensus_function_theta))

  #Important note: there may be multiple theta values that tie for the CMC max
  #count. The default of the cmcR package is to take the median of these theta
  #values as the thetaMax value
  thetaMax <- purrr::map(cmcPerTheta,~ calcMaxCMCTheta(cmcPerTheta = .,
                                                       highCMC_thresh = 1,
                                                       theta_thresh = theta_thresh,
                                                       distanceToCMCMaxTieBreaker = "minDistance"))

  if(purrr::is_empty(thetaMax$comparison_1to2) & purrr::is_empty(thetaMax$comparison_2to1)){
    thetaMax$comparison_2to1 <- NA
    thetaMax$comparison_1to2 <- NA
  }

  #if neither direction passed the high CMC criterion...
  if(all(is.na(thetaMax))){
    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = data.frame(cellNum = integer(0),
                                        cellRange = character(0),
                                        ccf = double(0),
                                        fft.ccf = double(0),
                                        x = integer(0),
                                        y = integer(0),
                                        theta = integer(0),
                                        comparison = character(0))))
  }

  #if one direction didn't pass the high CMC criterion...
  if(any(is.na(thetaMax))){

    #Most "liberal" decision is to replace the missing theta value with the
    #opposite of the other theta value
    if(missingTheta_decision == "replace"){
      thetaMax[[which(is.na(thetaMax))]] <- -1*thetaMax[[which(!is.na(thetaMax))]]
    }

    #More "moderate" decision is to dismiss the missing direction and only take
    #the initial CMCs defined for that direction
    else if(missingTheta_decision == "dismiss"){
      #the direction that didn't pass gets assigned initial CMCs:
      highCMCs1 <- originalMethodCMCs[[which(is.na(thetaMax))]] %>%
        dplyr::ungroup() %>%
        dplyr::mutate(comparison = rep(names(cmcPerTheta)[which(is.na(thetaMax))],times = nrow(.)))

      #the direction that passed the high CMC criterion gets all of its high
      #CMCs
      highCMCs2 <- purrr::pmap(.l = list(cmcPerTheta[which(!is.na(thetaMax))],
                                         names(cmcPerTheta)[which(!is.na(thetaMax))],
                                         thetaMax[which(!is.na(thetaMax))]),
                               function(cmcs,compName,th){
                                 purrr::map_dfr(th,~ dplyr::filter(cmcs,theta >= . - theta_thresh & theta <= . + theta_thresh)) %>%
                                   dplyr::mutate(comparison = rep(compName,times = nrow(.)))
                               })

      highCMCs <- highCMCs2 %>%
        dplyr::bind_rows() %>%
        dplyr::bind_rows(highCMCs1) %>%
        dplyr::distinct() %>%
        dplyr::group_by(cellNum) %>% #we don't want a cell being double-counted between the two comparisons
        dplyr::filter(ccf == max(ccf, na.rm = TRUE)) %>%
        dplyr::ungroup()

      #we want to make sure that the modal theta value in one direction is the
      #opposite (or close to the opposite) of the modal theta value in the other
      #direction
      thetaMax_dismissed <- highCMCs %>%
        dplyr::group_by(comparison,theta) %>%
        dplyr::tally() %>%
        dplyr::filter(n == max(n))

      #it's theoretically possible, albeit improbable, that one direction will
      #pass the high CMC criterion while the other direction fails *and*
      #produces 0 initial CMCs. In this case, we would only take the high CMCs
      #in the one direction. Not including this if statement first would throw
      #an error in the next if statement
      if(nrow(thetaMax_dismissed) == 1){
        return(list("params" = list(consensus_function = consensus_function,
                                    ccf_thresh = ccf_thresh,
                                    dx_thresh = dx_thresh,
                                    dy_thresh = dy_thresh,
                                    theta_thresh = theta_thresh,
                                    consensus_function_theta = consensus_function_theta),
                    "originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = highCMCs))
      }

      #another possibility is that more than one theta in one direction ties for
      #the CMC max count -- we can again determine whether these theta values
      #are "close" to each other (if not, then fail the criterion) and otherwise
      #take their median. Note that if there are three
      else if(nrow(thetaMax_dismissed) > 2){

        thetaMax_comparison_1to2_diff <- thetaMax_dismissed %>%
          dplyr::filter(comparison == "comparison_1to2") %>%
          dplyr::arrange(theta) %>%
          dplyr::pull(theta) %>%
          diff()

        thetaMax_comparison_2to1_diff <- thetaMax_dismissed %>%
          dplyr::filter(comparison == "comparison_2to1") %>%
          dplyr::arrange(theta) %>%
          dplyr::pull(theta) %>%
          diff()

        comparison_1to2_failure <- any(thetaMax_comparison_1to2_diff > theta_thresh)

        comparison_2to1_failure <- any(thetaMax_comparison_2to1_diff > theta_thresh)

        #in the event of a failure, don't assign any high CMCs

        if(comparison_1to2_failure | comparison_2to1_failure){
          return(list("params" = list(consensus_function = consensus_function,
                                      ccf_thresh = ccf_thresh,
                                      dx_thresh = dx_thresh,
                                      dy_thresh = dy_thresh,
                                      theta_thresh = theta_thresh,
                                      consensus_function_theta = consensus_function_theta),
                      "originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = data.frame(cellNum = integer(0),
                                              cellRange = character(0),
                                              ccf = double(0),
                                              fft.ccf = double(0),
                                              x = integer(0),
                                              y = integer(0),
                                              theta = integer(0),
                                              comparison = character(0))))
        }
        #otherwise, take the median of the theta values within their respective
        #directions
        else{
          thetaMax_dismissed <- thetaMax_dismissed %>%
            dplyr::group_by(comparison) %>%
            dplyr::summarise(theta = median(theta))
        }
      }
      #if thetaMax_dismissed has length 2, then we want to make sure that these
      #are actually opposites or close to opposites of each other. If not, then
      #we won't assign any high CMCs to the comparison
      thetaMax_dismissed <- thetaMax_dismissed %>%
        dplyr::ungroup() %>%
        dplyr::arrange(comparison) #%>%
      # dplyr::pull(theta)

      if(nrow(thetaMax_dismissed) == 1){

        if(compareInitialAndHighThetas){

          intitialCMCs_nonMissingDirection <- originalMethodCMCs[[which(names(originalMethodCMCs) == thetaMax_dismissed$comparison)]]

          thetaCompareBool_dismissed <- abs(thetaMax_dismissed$theta - median(intitialCMCs_nonMissingDirection$theta,na.rm = TRUE)) > theta_thresh

          if(thetaCompareBool_dismissed | is.na(thetaCompareBool_dismissed) | purrr::is_empty(thetaCompareBool_dismissed)){
            return(list("params" = list(consensus_function = consensus_function,
                                        ccf_thresh = ccf_thresh,
                                        dx_thresh = dx_thresh,
                                        dy_thresh = dy_thresh,
                                        theta_thresh = theta_thresh,
                                        consensus_function_theta = consensus_function_theta),
                        "originalMethodCMCs" = originalMethodCMCs,
                        "highCMCs" = data.frame(cellNum = integer(0),
                                                cellRange = character(0),
                                                ccf = double(0),
                                                fft.ccf = double(0),
                                                x = integer(0),
                                                y = integer(0),
                                                theta = integer(0),
                                                comparison = character(0))))

          }
          else{
            return(list("params" = list(consensus_function = consensus_function,
                                        ccf_thresh = ccf_thresh,
                                        dx_thresh = dx_thresh,
                                        dy_thresh = dy_thresh,
                                        theta_thresh = theta_thresh,
                                        consensus_function_theta = consensus_function_theta),
                        "originalMethodCMCs" = originalMethodCMCs,
                        "highCMCs" = highCMCs))
          }
        }
        else{
          return(list("params" = list(consensus_function = consensus_function,
                                      ccf_thresh = ccf_thresh,
                                      dx_thresh = dx_thresh,
                                      dy_thresh = dy_thresh,
                                      theta_thresh = theta_thresh,
                                      consensus_function_theta = consensus_function_theta),
                      "originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = highCMCs))
        }
      }

      if(compareInitialAndHighThetas){
        thetaCompareBool_dismissed <- (((abs((thetaMax_dismissed[1,"theta"] - median(originalMethodCMCs$comparison_1to2$theta,na.rm = TRUE)))) > theta_thresh) |
                                         (abs((thetaMax_dismissed[2,"theta"] - median(originalMethodCMCs$comparison_2to1$theta,na.rm = TRUE))) > theta_thresh))

        if(thetaCompareBool_dismissed | is.na(thetaCompareBool_dismissed) | purrr::is_empty(thetaCompareBool_dismissed)){
          return(list("params" = list(consensus_function = consensus_function,
                                      ccf_thresh = ccf_thresh,
                                      dx_thresh = dx_thresh,
                                      dy_thresh = dy_thresh,
                                      theta_thresh = theta_thresh,
                                      consensus_function_theta = consensus_function_theta),
                      "originalMethodCMCs" = originalMethodCMCs,
                      "highCMCs" = data.frame(cellNum = integer(0),
                                              cellRange = character(0),
                                              ccf = double(0),
                                              fft.ccf = double(0),
                                              x = integer(0),
                                              y = integer(0),
                                              theta = integer(0),
                                              comparison = character(0))))
        }
      }

      sign(thetaMax$comparison_1to2) == sign(thetaMax$comparison_2to1) & abs(thetaMax$comparison_1to2 - -1*thetaMax$comparison_2to1) > theta_thresh

      if((sign(thetaMax_dismissed[1,"theta"]) == sign(thetaMax_dismissed[2,"theta"]) & sign(thetaMax_dismissed[1,"theta"]) != 0 & sign(thetaMax_dismissed[2,"theta"]) != 0) |
         (abs((abs(thetaMax_dismissed[1,"theta"]) - abs(thetaMax_dismissed[2,"theta"]))) > theta_thresh)){
        return(list("params" = list(consensus_function = consensus_function,
                                    ccf_thresh = ccf_thresh,
                                    dx_thresh = dx_thresh,
                                    dy_thresh = dy_thresh,
                                    theta_thresh = theta_thresh,
                                    consensus_function_theta = consensus_function_theta),
                    "originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = data.frame(cellNum = integer(0),
                                            cellRange = character(0),
                                            ccf = double(0),
                                            fft.ccf = double(0),
                                            x = integer(0),
                                            y = integer(0),
                                            theta = integer(0),
                                            comparison = character(0))))
      }

      # if(sign(thetaMax_dismissed[1,"theta"]) == sign(thetaMax_dismissed[2,"theta"]) &
      #     abs(thetaMax_dismissed[1,"theta"] - -1*thetaMax_dismissed[2,"theta"]) > theta_thresh){
      # return(list("params" = list(consensus_function = consensus_function,
      #                             ccf_thresh = ccf_thresh,
      #                             dx_thresh = dx_thresh,
      #                             dy_thresh = dy_thresh,
      #                             theta_thresh = theta_thresh,
      #                             consensus_function_theta = consensus_function_theta),
      #             "originalMethodCMCs" = originalMethodCMCs,
      #             "highCMCs" = data.frame(cellNum = integer(0),
      #                                     cellRange = character(0),
      #                                     ccf = double(0),
      #                                     fft.ccf = double(0),
      #                                     x = integer(0),
      #                                     y = integer(0),
      #                                     theta = integer(0),
      #                                     comparison = character(0))))
      # }
      # else if((abs((abs(thetaMax_dismissed[1,"theta"]) - abs(thetaMax_dismissed[2,"theta"]))) > theta_thresh)){
      #   return(list("params" = list(consensus_function = consensus_function,
      #                               ccf_thresh = ccf_thresh,
      #                               dx_thresh = dx_thresh,
      #                               dy_thresh = dy_thresh,
      #                               theta_thresh = theta_thresh,
      #                               consensus_function_theta = consensus_function_theta),
      #               "originalMethodCMCs" = originalMethodCMCs,
      #               "highCMCs" = data.frame(cellNum = integer(0),
      #                                       cellRange = character(0),
      #                                       ccf = double(0),
      #                                       fft.ccf = double(0),
      #                                       x = integer(0),
      #                                       y = integer(0),
      #                                       theta = integer(0),
      #                                       comparison = character(0))))
      # }
      # if we've made it this far, then the theta values should be at least to
      # within theta_thresh of being opposites of each other, so we can return
      # the highCMCs without worrying about a disagreement
      else{
        return(list("params" = list(consensus_function = consensus_function,
                                    ccf_thresh = ccf_thresh,
                                    dx_thresh = dx_thresh,
                                    dy_thresh = dy_thresh,
                                    theta_thresh = theta_thresh,
                                    consensus_function_theta = consensus_function_theta),
                    "originalMethodCMCs" = originalMethodCMCs,
                    "highCMCs" = highCMCs))
      }
    }

    # Most "conservative" decision is to flat-out fail the whole comparison if
    # one direction doesn't pass the high CMC criterion
    else if(missingTheta_decision == "fail"){
      return(list("params" = list(consensus_function = consensus_function,
                                  ccf_thresh = ccf_thresh,
                                  dx_thresh = dx_thresh,
                                  dy_thresh = dy_thresh,
                                  theta_thresh = theta_thresh,
                                  consensus_function_theta = consensus_function_theta),
                  "originalMethodCMCs" = originalMethodCMCs,
                  "highCMCs" = data.frame(cellNum = integer(0),
                                          cellRange = character(0),
                                          ccf = double(0),
                                          fft.ccf = double(0),
                                          x = integer(0),
                                          y = integer(0),
                                          theta = integer(0),
                                          comparison = character(0))))
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

  if(compareInitialAndHighThetas){

    thetaCompareBool <- (((abs((thetaMax$comparison_1to2 - median(originalMethodCMCs$comparison_1to2$theta,na.rm = TRUE)))) > theta_thresh) |
                           (abs((thetaMax$comparison_2to1 - median(originalMethodCMCs$comparison_2to1$theta,na.rm = TRUE))) > theta_thresh))

    if(thetaCompareBool | is.na(thetaCompareBool) | purrr::is_empty(thetaCompareBool)){
      return(list("params" = list(consensus_function = consensus_function,
                                  ccf_thresh = ccf_thresh,
                                  dx_thresh = dx_thresh,
                                  dy_thresh = dy_thresh,
                                  theta_thresh = theta_thresh,
                                  consensus_function_theta = consensus_function_theta),
                  "originalMethodCMCs" = originalMethodCMCs,
                  "highCMCs" = data.frame(cellNum = integer(0),
                                          cellRange = character(0),
                                          ccf = double(0),
                                          fft.ccf = double(0),
                                          x = integer(0),
                                          y = integer(0),
                                          theta = integer(0),
                                          comparison = character(0))))

    }
  }

  #The last contingency is making sure that these theta modes are opposites of
  #each other (or within theta_thresh of being opposites).

  #The following determines whether (1) the theta modes in either direction are
  #the same sign as each other (which they shouldn't be for matches) AND (2)
  #whether they are far away from 0. An older version of this function didn't
  #cover cases in which, for example, the theta_thresh was 6 and both
  #comparisons voted for a consensual theta value of 3. Given the theta_thresh,
  #this should be admissable since they are "close enough" (according to
  #theta_thresh) to each other.
  if(sign(thetaMax$comparison_1to2) == sign(thetaMax$comparison_2to1) & abs(thetaMax$comparison_1to2 - -1*thetaMax$comparison_2to1) > theta_thresh){

    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = data.frame(cellNum = integer(0),
                                        cellRange = character(0),
                                        ccf = double(0),
                                        fft.ccf = double(0),
                                        x = integer(0),
                                        y = integer(0),
                                        theta = integer(0),
                                        comparison = character(0))))
  }
  #Even if the two theta values are of opposite signs, it's possible that they
  #don't agree with each other. For example, one direction might vote for -27
  #degrees as the consensual theta value while the other votes for 12. Such a
  #comparison shouldn't pass the High CMC criterion.
  else if((abs((abs(thetaMax$comparison_1to2) - abs(thetaMax$comparison_2to1))) > theta_thresh)){
    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = data.frame(cellNum = integer(0),
                                        cellRange = character(0),
                                        ccf = double(0),
                                        fft.ccf = double(0),
                                        x = integer(0),
                                        y = integer(0),
                                        theta = integer(0),
                                        comparison = character(0))))
  }
  else{
    highCMCs <- purrr::pmap(.l = list(cmcPerTheta,
                                      names(cmcPerTheta),
                                      thetaMax),
                            function(cmcs,compName,th){
                              purrr::map_dfr(th,~ dplyr::filter(cmcs,theta >= . - theta_thresh & theta <= . + theta_thresh)) %>%
                                dplyr::mutate(comparison = rep(compName,times = nrow(.)))
                            }) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct() %>%
      dplyr::group_by(cellNum) %>% #we don't want a cell being double-counted between the two comparisons
      dplyr::filter(ccf == max(ccf, na.rm = TRUE)) %>%
      dplyr::ungroup()

    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "originalMethodCMCs" = originalMethodCMCs,
                "highCMCs" = highCMCs))
  }
}