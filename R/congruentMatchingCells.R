#' Identifies the congruent matching cells between two processed cartridge case
#' scans
#'
#' @name congruentMatchingCells
#'
#' @description Wrapper for calling the cellCCF_bothDirections and
#'   cmcFilter_improved functions together
#'
#' @param x3p1 (no default) an x3p object containing the surface matrix of a
#'   breech face impression
#' @param x3p2 (no default) an x3p object containing the surface matrix of a
#'   breech face impression to be compared to that in x3p1
#' @param thetas (default seq(from = -30,to = 30,by = 3)) rotation values (in
#'   degrees) for which x3p2$surface.matrix will be rotated, split into cells,
#'   and compared to x3p1$surface.matrix
#' @param cellNumHoriz (default 7) number of cells along horizontal axis to
#'   divide x3p1$surface.matrix into
#' @param cellNumVert (default equal to cellNumHoriz) number of cells along
#'   vertical axis to divide x3p1$surface.matrix into
#' @param regionToCellProp determines how much larger the x3p2 regions will be
#'   relative to the x3p1 cells. For example, if regionToCellProp = 4 means that
#'   the x3p2 regions will be 4 times times larger (sidelengths multiplied by 2)
#' @param minObservedProp (default 15) the minimum proportion of a cell that
#'   needs to contain observed (i.e., non-NA) values for it to be included in
#'   the CCF calculation procedure
#' @param centerCell **OPTIONAL** (default missing) dictates if cell is to be
#'   shifted prior to CCF calculation. Default is that no shifting is performed.
#'   If set to "wholeMatrix", then each cell is subracted by the mean of the
#'   whole matrix. If set to "individualCell", then each cell is subtracted by
#'   its particular mean.
#' @param scaleCell **OPTIONAL** (default missing) dictates if cell is to be
#'   scaled prior to CCF calculation. Default is that no scaling is performed.
#'   If set to "wholeMatrix", then each cell is divided by the standard
#'   deviation of the whole matrix. If set to "individualCell", then each cell
#'   is divided by its particular standard deviation.
#' @param cellCCF_bothDirections_output list returned by the function
#'   cmcR::cellCCF_bothdirections
#' @param consensus_function function to aggregate the translation (dx and dy)
#'   and rotation (theta) values in the ccfDF data frame to determine
#'   "consensus" values
#' @param corr_thresh minimum correlation threshold to call a cell pair
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
#' #x3p1 and x3p2 are two processed cartridge case scans
#' cmc <- congruentMatchingCells(x3p1,x3p2)
#' }
#'
#' @seealso cmcR::cellCCF_bothDirections and cmcR::cmcFilter_improved
#' @seealso
#' \url{https://pdfs.semanticscholar.org/4bf3/0b3a23c38d8396fa5e0d116cba63a3681494.pdf}
#' @seealso
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf}
#'

congruentMatchingCells <- function(x3p1,
                                   x3p2,
                                   thetas = seq(-30,30,by = 3),
                                   cellNumHoriz = 7,
                                   regionToCellProp = 4,
                                   cellNumVert = cellNumHoriz,
                                   centerCell,
                                   scaleCell,
                                   consensus_function = median,
                                   corr_thresh = .4,
                                   dx_thresh = 20,
                                   dy_thresh = dx_thresh,
                                   theta_thresh = 3,
                                   consensus_function_theta = consensus_function,...){

  ccfResults <- cmcR::cellCCF_bothDirections(x3p1 = x3p1,
                                             x3p2 = x3p2,
                                             thetas = thetas,
                                             cellNumHoriz = cellNumHoriz,
                                             cellNumVert = cellNumVert,
                                             regionToCellProp = regionToCellProp,
                                             centerCell = centerCell,
                                             scaleCell = scaleCell)

  cmcs <- cmcR::cmcFilter_improved(cellCCF_bothDirections_output = ccfResults,
                                   consensus_function = consensus_function,
                                   corr_thresh = corr_thresh,
                                   dx_thresh = dx_thresh,
                                   dy_thresh = dy_thresh,
                                   theta_thresh = theta_thresh,
                                   consensus_function_theta = consensus_function_theta,...)

  return(cmcs)
}