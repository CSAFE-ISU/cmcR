#' Compiles a data frame
#'
#' @name topResultsPerCell
#'
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

#' Applies Congruent Matching Cells method logic to the CCF results of a
#' comparison between two cartridge case scans.
#'
#' @name cmcFilter
#'
#' @param ... arguments to be passed to consensus_function or
#'   consensus_function_theta if necessary
#' @export

cmcFilter <- function(ccfResults,
                      consensus_function = median,
                      ccf_thresh = .4,
                      dx_thresh = 20,
                      dy_thresh = dx_thresh,
                      theta_thresh = 3,
                      consensus_function_theta = consensus_function,...){
  #Required tests:
  # ccfResults needs to be a dataframe containing columns: ccf,theta,dx,dy
  # ccf_thresh should be between 0 and 1
  # dx_thresh should positive
  # consensus_function and consensus_function_theta should be a function name that exists

  consensus_dx <- consensus_function(ccfResults$dx,...)
  consensus_dy <- consensus_function(ccfResults$dy,...)
  consensus_theta <- consensus_function_theta(ccfResults$theta,...)

  ccfResults %>%
    dplyr::filter(ccf >= ccf_thresh &
                    dx >= consensus_dx - dx_thresh & dx <= consensus_dx + dx_thresh &
                    dy >= consensus_dy - dy_thresh & dy <= consensus_dy + dy_thresh &
                    theta >= consensus_theta - theta_thresh & theta <= consensus_theta + theta_thresh)
}

#' @name cmcFilterPerTheta
#' @param ... arguments to be passed to consensus_function or
#'   consensus_function_theta if necessary

cmcFilterPerTheta <- function(ccfResults,
                              consensus_function = median,
                              ccf_thresh = .4,
                              dx_thresh = 20,
                              dy_thresh = dx_thresh,
                              theta_thresh = 3,
                              consensus_function_theta = consensus_function,...){

  ccfResults <- ccfResults %>%
    purrr::map2(.x = .,
                .y = names(.),
                function(thetaSpecificResults,theta){
                  thetaSpecificResults %>%
                    dplyr::mutate(theta = rep(as.numeric(theta),times = nrow(.)))
                }) %>%
    setNames(names(ccfResults))

  ccfResults %>%
    purrr::map(~ cmcFilter(ccfResults = .,
                           consensus_function = consensus_function,
                           ccf_thresh = ccf_thresh,
                           dx_thresh = dx_thresh,
                           dy_thresh = dy_thresh,
                           theta_thresh = theta_thresh,
                           consensus_function_theta = median)) %>%
    dplyr::bind_rows()
}

#' @name getMode
#'
#' @description Calculates the mode of a vector. Can be used as a consensus
#'   function in cmcR::cmcFilter.
#'
#' @export

getMode <- function(x){
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#' @name calcMaxCMCTheta

calcMaxCMCTheta <- function(cmcPerTheta,
                            highCMC_thresh = 1,
                            theta_thresh = 3){
  cmcCountPerTheta <- cmcPerTheta %>%
    dplyr::group_by(theta) %>%
    dplyr::tally()

  cmcMax <- cmcCountPerTheta %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n))

  #there may be more than one theta tied with the maximum CMC - in which case
  #we should determine if such thetas are "close" to each other, where
  #proximity is defined based on the theta_tresh set. If not, then we can rule
  #out that the comparison under scrutiny is a known match.
  if(any(diff(cmcMax$theta) > theta_thresh)){
    return(NA)
  }
  #
  maxDistancetoCMCMax <- cmcCountPerTheta %>%
    dplyr::filter(n >= cmcMax$n - highCMC_thresh) %>%
    dplyr::group_by(theta) %>%
    #if there may be multiple theta values tied for cmcMax that are within
    #theta_thresh of each other, then we can't discount any of them as being
    #the "true" theta value (alternatively, may be we should discount all as
    #NOT being the true theta value?). Instead we'll consider how far the
    #"high CMC" theta values are from the closest tied maxCMC theta value.
    dplyr::summarise(distanceToCMCMax = min(abs(cmcMax$theta - theta))) %>%
    dplyr::pull(distanceToCMCMax) %>%
    max()

  if(maxDistancetoCMCMax > theta_thresh){
    return(NA)
  }
  else{
    return(cmcMax$theta)
  }
}

#'
#'
#' @name cmcFilter_improved
#'
#' @param cellCCF_bothDirections_output list returned by the function
#'   cmcR::cellCCF_bothdirections
#'
#' @export

cmcFilter_improved <- function(cellCCF_bothDirections_output,
                               consensus_function = median,
                               ccf_thresh = .4,
                               dx_thresh = 20,
                               dy_thresh = dx_thresh,
                               theta_thresh = 3,
                               consensus_function_theta = consensus_function,...){

  initialCMCs <- cellCCF_bothDirections_output %>%
    purrr::map(~ cmcR::topResultsPerCell(.$ccfResults) %>%
          cmcR:::cmcFilter(consensus_function = consensus_function,
                           ccf_thresh = ccf_thresh,
                           dx_thresh = dx_thresh,
                           dy_thresh = dy_thresh,
                           theta_thresh = theta_thresh,
                           consensus_function_theta = consensus_function_theta))

  cmcPerTheta <-  cellCCF_bothDirections_output %>%
    purrr::map(~ cmcR:::cmcFilterPerTheta(ccfResults = .$ccfResults,
                                   consensus_function = consensus_function,
                                   ccf_thresh = ccf_thresh,
                                   dx_thresh = dx_thresh,
                                   dy_thresh = dy_thresh,
                                   theta_thresh = theta_thresh,
                                   consensus_function_theta = consensus_function_theta,...))

  thetaMax <- purrr::map(cmcPerTheta,cmcR:::calcMaxCMCTheta)

  if(all(is.na(thetaMax))){
    # print("Note: neither comparison produces a valid max CMC theta value. The initial CMCs based on the top results per cell will be returned.")
    return(list("params" = list(consensus_function = consensus_function,
                                ccf_thresh = ccf_thresh,
                                dx_thresh = dx_thresh,
                                dy_thresh = dy_thresh,
                                theta_thresh = theta_thresh,
                                consensus_function_theta = consensus_function_theta),
                "initialCMCs" = list(initialCMCs)))
  }

  # if(is.na(thetaMax$comparison_1to2 != -thetaMax$comparison_2to1) | thetaMax$comparison_1to2 != -thetaMax$comparison_2to1){
  #   # print(paste0("Note: max CMC thetas disagree. Comparison of x3p1 to x3p2: ",thetaMax$comparison_1to2," degrees vs. Comparison of x3p2 to x3p1: ",thetaMax$comparison_2to1," degrees. If one is NA, then it will be replaced with the opposite of the other for final CMC calculation."))
  # }

  if(is.na(thetaMax$comparison_1to2)){
    thetaMax$comparison_1to2 <- -thetaMax$comparison_2to1
  }
  if(is.na(thetaMax$comparison_2to1)){
    thetaMax$comparison_2to1 <- -thetaMax$comparison_1to2
  }

  finalCMCs <- purrr::pmap(.l = list(cmcPerTheta,
                                     names(cmcPerTheta),
                                     thetaMax),
                           function(cmcs,compName,th){
                             cmcs %>%
                               dplyr::filter(theta >= th - 3 & theta <= th + 3) %>%
                               dplyr::mutate(comparison = rep(compName,times = nrow(.)))
                           }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(cellID) %>%
    dplyr::filter(ccf == max(ccf))

  return(list("params" = list(consensus_function = consensus_function,
                              ccf_thresh = ccf_thresh,
                              dx_thresh = dx_thresh,
                              dy_thresh = dy_thresh,
                              theta_thresh = theta_thresh,
                              consensus_function_theta = consensus_function_theta),
              "initialCMCs" = list(initialCMCs),
              "finalCMCs" = finalCMCs))
}