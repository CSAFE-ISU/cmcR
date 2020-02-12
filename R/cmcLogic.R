#' @name topResultsPerCell
#'
#' @export

topResultsPerCell <- function(cellCCF_results){
  cellCCF_results %>%
    map2(.x = ,
         .y = names(.),
         function(corrValues,pairName){
           corrValues %>%
             map2_dfr(.x = .,
                      .y = names(.),
                      ~ .x %>%
                        mutate(theta = rep(.y,times = nrow(.)))) %>%
             mutate(pairName = rep(pairName,times = nrow(.))) %>%
             group_by(cellID) %>%
             filter(corr == max(corr))
         }) %>%
    bind_rows()
}