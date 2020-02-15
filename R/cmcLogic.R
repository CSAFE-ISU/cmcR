#' @name topResultsPerCell
#'
#' @export

topResultsPerCell <- function(cellCCF_results,...){
  cellCCF_results  %>%
    purrr::map2_dfr(.x = .,
                    .y = names(.),
                    ~ .x %>%
                      dplyr::mutate(theta = as.numeric(rep(.y,times = nrow(.))))) %>%
    dplyr::group_by(cellID) %>%
    dplyr::filter(corr == max(corr)) %>%
    dplyr::arrange(cellID)
}

#' @name countInitialCMCs
#'
#' @export

countInitialCMCs <- function(cellCCF_results,
                             consensus_function = median,
                             ccf_thresh = .4,
                             dx_thresh = 20,
                             dy_thresh = dx_thresh,
                             theta_thresh = 3,...){
  topResults <- cmcR::topResultsPerCell(cellCCF_results = cellCCF_results)

  consensus_dx <- consensus_function(topResults$dx)
  consensus_dy <- consensus_function(topResults$dy)
  consensus_theta <- consensus_function(topResults$theta)

  topResults %>%
    dplyr::filter(corr >= ccf_thresh &
                    dx >= consensus_dx - dx_thresh & dx <= consensus_dx + dx_thresh &
                    dy >= consensus_dy - dy_thresh & dy <= consensus_dy + dy_thresh &
                    theta >= consensus_theta - theta_thresh & theta <= consensus_theta + theta_thresh)
}

