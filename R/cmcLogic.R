#' Calculate dx, dy, theta values at which top similarity is achieved in a cell
#' pair.
#'
#' @name topResultsPerCell
#'
#' @description Given a list of data frames like the one returned in the
#'   `ccfResults` element of the list returned from cmcR::cellCCF, returns a
#'   data frame of the theta, dx, dy, ccf, and ccf values at which each cell
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

topResultsPerCell <- function(ccfResults){

  ccfResults  %>%
    purrr::map2_dfr(.x = .,
                    .y = names(.),
                    ~ .x %>%
                      dplyr::mutate(theta = as.numeric(rep(.y,times = nrow(.))))) %>%
    dplyr::group_by(cellID) %>%
    dplyr::filter(ccf == max(ccf)) %>%
    dplyr::arrange(cellID)
}

#' Implements CMC logic on a single data frame as proposed by Song (2013).
#'
#' @name cmcFilter
#'
#' @description Applies initially proposed Congruent Matching Cells method logic to the CCF
#'   results of a comparison between two cartridge case scans.
#'
#' @param ccfDF data frame containing ccf results from a comparison between two
#'   cartridge case scans
#' @param consensus_function function to aggregate the translation (dx and dy)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus dx value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus dy value that a cell
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

cmcFilter <- function(ccfDF,
                      consensus_function = median,
                      ccf_thresh = .6,
                      dx_thresh = 10,
                      dy_thresh = dx_thresh,
                      theta_thresh = 3,
                      consensus_function_theta = consensus_function,...){
  #Required tests:
  # ccfResults needs to be a dataframe containing columns: ccf,theta,dx,dy
  # ccf_thresh should be between 0 and 1
  # dx_thresh should positive
  # consensus_function and consensus_function_theta should be a function name that exists

  consensus_dx <- consensus_function(ccfDF$dx,...)
  consensus_dy <- consensus_function(ccfDF$dy,...)
  consensus_theta <- consensus_function_theta(ccfDF$theta,...)

  ccfDF %>%
    dplyr::filter(ccf >= ccf_thresh &
                    dx >= consensus_dx - dx_thresh & dx <= consensus_dx + dx_thresh &
                    dy >= consensus_dy - dy_thresh & dy <= consensus_dy + dy_thresh &
                    theta >= consensus_theta - theta_thresh & theta <= consensus_theta + theta_thresh)
}

#' @name cmcFilterPerTheta
#'
#' @param ccfResults list of data frames, like the one returned by
#'   cmcR::cellCCF, containing breech face cell comparison results
#' @param consensus_function function to aggregate the translation (dx and dy)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus dx value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus dy value that a cell
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
                    select(-theta) %>%
                    dplyr::mutate(theta = rep(as.numeric(theta),times = nrow(.)))
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
  else{
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
    return(cmcMax$theta)
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
#' @param consensus_function function to aggregate the translation (dx and dy)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param ccf_thresh minimum correlation threshold to call a cell pair
#'   "congruent matching"
#' @param dx_thresh maximum distance from the consensus dx value that a cell
#'   pair can be to be called "congruent matching"
#' @param dy_thresh  maximum distance from the consensus dy value that a cell
#'   pair can be to be called "congruent matching"
#' @param theta_thresh maximum distance from the consensus theta value that a
#'   cell pair can be to be called "congruent matching"
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
#'                            consensus_function_theta = getMode)
#' }
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'
#' @export

cmcFilter_improved <- function(cellCCF_bothDirections_output,
                               consensus_function = median,
                               ccf_thresh = .6,
                               dx_thresh = 10,
                               dy_thresh = dx_thresh,
                               theta_thresh = 3,
                               consensus_function_theta = consensus_function){

  #TODO (potentially): Add an argument "thetaDisagreement" that allows for more
  #conservative decisions if the theta modes are identified in both directions
  #but they don't agree with each other (aren't within theta_thresh of each
  #other). Could only take the min of the Final CMCs identified in either
  #direction or only give that pair its initial CMCs (i.e., it would "fail" the
  #both direction criterion).

  initialCMCs <- cellCCF_bothDirections_output %>%
    purrr::map(~ cmcR::topResultsPerCell(.$ccfResults) %>%
                 cmcR:::cmcFilter(consensus_function = consensus_function,
                                  ccf_thresh = ccf_thresh,
                                  dx_thresh = dx_thresh,
                                  dy_thresh = dy_thresh,
                                  theta_thresh = theta_thresh,
                                  consensus_function_theta = consensus_function_theta) %>%
                 ungroup())

  cmcPerTheta <-  cellCCF_bothDirections_output %>%
    purrr::map(~ cmcR:::cmcFilterPerTheta(ccfResults = .$ccfResults,
                                          consensus_function = consensus_function,
                                          ccf_thresh = ccf_thresh,
                                          dx_thresh = dx_thresh,
                                          dy_thresh = dy_thresh,
                                          theta_thresh = theta_thresh,
                                          consensus_function_theta = consensus_function_theta))

  thetaMax <- purrr::map(cmcPerTheta,~ cmcR:::calcMaxCMCTheta(cmcPerTheta = .,
                                                              highCMC_thresh = 1,
                                                              theta_thresh = theta_thresh,
                                                              distanceToCMCMaxTieBreaker = "minDistance"))

  #This was the least "conservative" option if one direction didn't yield a
  #theta mode but the other direction did. I think this was assigning to many "false positive" CMCs for our liking for KNMs.
  #if(purrr::is_empty(thetaMax$comparison_1to2) &
  #!purrr::is_empty(thetaMax$comparison_2to1)){ thetaMax$comparison_1to2 <-
  #-thetaMax$comparison_2to1 } if(!purrr::is_empty(thetaMax$comparison_1to2) &
  #purrr::is_empty(thetaMax$comparison_2to1)){ thetaMax$comparison_2to1 <-
  #-thetaMax$comparison_1to2 }
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
                "initialCMCs" = initialCMCs,
                "finalCMCs" = NULL))
  }

  #if one direction didn't pass the high CMC criterion...
  if(any(is.na(thetaMax))){
    #the direction that didn't pass gets assigned initial CMCs:
    finalCMCs1 <- initialCMCs[[which(is.na(thetaMax))]] %>%
      ungroup() %>%
      mutate(comparison = rep(names(cmcPerTheta)[which(is.na(thetaMax))],times = nrow(.)))


    finalCMCs2 <- purrr::pmap(.l = list(cmcPerTheta[which(!is.na(thetaMax))],
                                            names(cmcPerTheta)[which(!is.na(thetaMax))],
                                            thetaMax[which(!is.na(thetaMax))]),
                                  function(cmcs,compName,th){
                                    purrr::map_dfr(th,~ dplyr::filter(cmcs,theta >= . - 3 & theta <= . + 3)) %>%
                                      dplyr::mutate(comparison = rep(compName,times = nrow(.)))
                                  })

    finalCMCs <- finalCMCs2 %>%
      dplyr::bind_rows() %>%
      dplyr::bind_rows(finalCMCs1) %>%
      dplyr::distinct() %>%
      dplyr::group_by(cellNum) %>% #we don't want a cell being double-counted between the two comparisons
      dplyr::filter(ccf == max(ccf)) %>%
      dplyr::ungroup()
  } #if both directions pass the high CMC criterion...
  else{
    finalCMCs <- purrr::pmap(.l = list(cmcPerTheta,
                                       names(cmcPerTheta),
                                       thetaMax),
                             function(cmcs,compName,th){
                               purrr::map_dfr(th,~ dplyr::filter(cmcs,theta >= . - 3 & theta <= . + 3)) %>%
                                 dplyr::mutate(comparison = rep(compName,times = nrow(.)))
                             }) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct() %>%
      dplyr::group_by(cellNum) %>% #we don't want a cell being double-counted between the two comparisons
      dplyr::filter(ccf == max(ccf)) %>%
      dplyr::ungroup()
  }

  thetaMed <- finalCMCs %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(theta = median(theta)) %>%
    dplyr::pull(theta)

  #one last check to determine if the median theta in both directions agree with
  #each other up to a sign within theta_thresh
  if(sign(thetaMed[1]) == sign(thetaMed[2]) | (abs((abs(thetaMed[1]) - abs(thetaMed[2]))) > theta_thresh)){
    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "initialCMCs" = initialCMCs,
                "finalCMCs" = NULL))
  }

  return(list("params" = list(consensus_function = consensus_function,
                              ccf_thresh = ccf_thresh,
                              dx_thresh = dx_thresh,
                              dy_thresh = dy_thresh,
                              theta_thresh = theta_thresh,
                              consensus_function_theta = consensus_function_theta),
              "initialCMCs" = initialCMCs,
              "finalCMCs" = finalCMCs))
}