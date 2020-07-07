#' Plot a list of x3ps
#' @name x3pListPlot
#'
#' @description Plots the surface matrices in a list of x3p objects. Either
#'   creates one plot faceted by surface matrix or creates individual plots per
#'   surface matrix and returns them in a list.
#'
#' @param x3pList a list of x3p objects. If the x3p objects are named in the
#'   list, then these names will be included in the title of their respective
#'   plot
#' @param type dictates whether one plot faceted by surface matrix or a list of
#'   plots per surface matrix is returned. The faceted plot will have a
#'   consistent height scale across all surface matrices.
#' @param rotate angle (in degrees) to rotate all surface matrices plotted
#' @param legend.quantiles vector of quantiles to be shown as tick marks on
#'   legend plot
#' @param height.colors vector of colors to be passed to scale_fill_gradientn
#'   that dictates the height value colorscale
#' @param na.value color to be used for NA values (passed to
#'   scale_fill_gradientn)
#' @param guide internal usage
#' @examples
#' \dontrun{
#'  x3pListPlot(list("name1" = x3p1, "name2" = x3p2))
#' }
#' @export
#'
#' @importFrom stats setNames median quantile

utils::globalVariables(c("value","x","y",".","x3p","height"))

x3pListPlot <- function(x3pList,
                        type = "faceted",
                        rotate = 0,
                        legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                        height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                        na.value = "gray80",
                        guide = "colorbar"){
  if(purrr::is_empty(names(x3pList))){
    x3pList <- setNames(x3pList,paste0("x3p",1:length(x3pList)))
  }

  if(type == "faceted"){
    surfaceMat_df <- purrr::pmap_dfr(.l = list(x3pList,
                                               names(x3pList),
                                               rotate),
                                     function(x3p,name,theta){
                                       x3p$surface.matrix <- rotateSurfaceMatrix(x3p$surface.matrix,
                                                                                 theta = theta + 180) #+180 to stay with what rotate_x3p would output

                                       x3p %>%
                                         x3ptools::x3p_to_df() %>%
                                         dplyr::mutate(value = value - median(value,na.rm = TRUE)) %>%
                                         dplyr::mutate(x = x*1e6,
                                                       y = y*1e6,
                                                       height = value*1e6) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = height))  +
      ggplot2::scale_fill_gradientn(colours = height.colors,
                                    values = scales::rescale(quantile(surfaceMat_df$height,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                    breaks = function(lims){
                                      dat <- quantile(surfaceMat_df$height,legend.quantiles,na.rm = TRUE)

                                      dat <- dat %>%
                                        setNames(paste0(names(dat)," [",round(dat,3),"]"))

                                      return(dat)
                                    },
                                    na.value = na.value,
                                    guide = guide) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank()) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(2.5,"in"),
                                                      label.theme = ggplot2::element_text(size = 8),
                                                      title.theme = ggplot2::element_text(size = 10),
                                                      frame.colour = "black",
                                                      ticks.colour = "black"),
                      colour = FALSE) +
      ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
      ggplot2::facet_wrap(~ x3p)
  }
  else if(type == "list"){
    plts <- purrr::pmap(.l = list(x3pList,
                                  names(x3pList),
                                  rotate),
                        function(x3p,name,theta){
                          x3p$surface.matrix <- rotateSurfaceMatrix(x3p$surface.matrix,
                                                                    theta = theta + 180) #+180 to stay with what rotate_x3p would output

                          surfaceMat_df <- x3p %>%
                            x3ptools::x3p_to_df() %>%
                            dplyr::mutate(value = value - median(value,na.rm = TRUE)) %>%
                            dplyr::mutate(x = x*1e6,
                                          y = y*1e6,
                                          height = value*1e6) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = height))  +
                            ggplot2::scale_fill_gradientn(colours = height.colors,
                                                          values = scales::rescale(quantile(surfaceMat_df$height,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                                          breaks = function(lims){
                                                            dat <- quantile(surfaceMat_df$height,legend.quantiles,na.rm = TRUE)

                                                            dat <- dat %>%
                                                              setNames(paste0(names(dat)," [",round(dat,3),"]"))

                                                            return(dat)
                                                          },
                                                          na.value = na.value,
                                                          guide = guide) +
                            ggplot2::theme_minimal() +
                            ggplot2::coord_fixed(expand = FALSE) +
                            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                           axis.text.x = ggplot2::element_blank(),
                                           axis.ticks.x = ggplot2::element_blank(),
                                           axis.title.y = ggplot2::element_blank(),
                                           axis.text.y = ggplot2::element_blank(),
                                           axis.ticks.y = ggplot2::element_blank(),
                                           panel.grid.major = ggplot2::element_blank(),
                                           panel.grid.minor = ggplot2::element_blank(),
                                           panel.background = ggplot2::element_blank()) +
                            ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(3,"in"),
                                                                            label.theme = ggplot2::element_text(size = 8),
                                                                            title.theme = ggplot2::element_text(size = 10),
                                                                            frame.colour = "black",
                                                                            ticks.colour = "black"),
                                            colour = FALSE) +
                            ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
                            ggplot2::ggtitle(name)
                        })
  }
  return(plts)
}

#' @name linear_to_matrix
#' @param index integer vector of indices, must be between 1 and nrow*ncol
#' @param nrow number of rows, integer value defaults to 7
#' @param ncol  number of columns, integer value, defaults to number of rows
#' @param byrow logical value, is linear index folded into matrix by row (default) or by column (`byrow=FALSE`).
#' @examples
#' index <- sample(nrow*ncol, 10, replace = TRUE)
#' linear_to_matrix(index, nrow=4, ncol = 5, byrow=TRUE)
#'
#' @keywords internal

linear_to_matrix <- function(index, nrow = 7, ncol = nrow, byrow = TRUE, sep = ", ") {
  index <- as.integer(index)
  stopifnot(all(index <= nrow*ncol), all(index > 0))
  idx <- 1:(nrow*ncol)
  if (byrow) { # column is the fast index
    idx_out_col <- ((index-1) %% ncol) + 1
    idx_out_row <- ((index-1) %/% ncol) + 1
  } else { # row is the fast index
    idx_out_col <- ((index-1) %/% nrow) + 1
    idx_out_row <- ((index-1) %% nrow) + 1
  }
  paste0(idx_out_row, sep, idx_out_col)
}

#' @name arrangeCMCPlot
#'
#' @keywords internal
#'
#' @importFrom stats median setNames
#' @importFrom ggnewscale new_scale_fill

utils::globalVariables(c("firstRow","lastRow","firstCol","lastCol","cellNum",".","firstColCentered","theta","firstRowCentered","dx","dy","lastColCentered","lastRowCentered","topRightCorner_col","bottomLeftCorner_col","topRightCorner_row","bottomLeftCorner_row","cmc","midCol","midRow","cellInd"))

arrangeCMCPlot <- function(x3p1,
                           x3p2,
                           allCells,
                           x3pNames,
                           pltType = "faceted",
                           legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                           height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           cell.colors = c("#a50026","#313695"),
                           na.value = "gray80"){

  x3p1_cellGrid <- allCells %>%
    dplyr::mutate(firstRow = (x3p1$header.info$incrementY*1e6)*(firstRow),
                  lastRow = (x3p1$header.info$incrementY*1e6)*(lastRow),
                  firstCol = (x3p1$header.info$incrementY*1e6)*(firstCol),
                  lastCol = (x3p1$header.info$incrementY*1e6)*(lastCol)) %>%
    dplyr::mutate(x_1 = firstCol,
                  y_1 = firstRow,
                  x_2 = lastCol,
                  y_2 = firstRow,
                  x_3 = lastCol,
                  y_3 = lastRow,
                  x_4 = firstCol,
                  y_4 = lastRow) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
                        names_to = c(".value","order"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(midCol = (lastCol + firstCol)/2,
                  midRow = (lastRow + firstRow)/2,
                  cellInd = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                                               floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                                               ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                                             nrow = ceiling(sqrt(max(cellNum))),
                                             byrow = TRUE),
                  x3p = rep(x3pNames[1],times = nrow(.)),
                  theta = rep(0,times = nrow(.)))

  x3p2_cellGrid <- allCells %>%
    dplyr::mutate(firstRow = (x3p2$header.info$incrementY*1e6)*(firstRow),
                  lastRow = (x3p2$header.info$incrementY*1e6)*(lastRow),
                  firstCol = (x3p2$header.info$incrementY*1e6)*(firstCol),
                  lastCol = (x3p2$header.info$incrementY*1e6)*(lastCol)) %>%
    dplyr::mutate(firstRowCentered = firstRow - max(lastRow)/2,
                  lastRowCentered = lastRow - max(lastRow)/2,
                  firstColCentered = firstCol - max(lastCol)/2,
                  lastColCentered = lastCol - max(lastCol)/2) %>%
    dplyr::mutate(topLeftCorner_col = firstColCentered*cos((theta - median(theta))*(pi/180)) - lastRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (x3p2$header.info$incrementY*1e6)*dx,
                  topLeftCorner_row = firstColCentered*sin((theta - median(theta))*(pi/180)) + lastRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (x3p2$header.info$incrementY*1e6)*dy,
                  topRightCorner_col = lastColCentered*cos((theta - median(theta))*(pi/180)) - lastRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (x3p2$header.info$incrementY*1e6)*dx,
                  topRightCorner_row = lastColCentered*sin((theta - median(theta))*(pi/180)) + lastRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (x3p2$header.info$incrementY*1e6)*dy,
                  bottomRightCorner_col = lastColCentered*cos((theta - median(theta))*(pi/180)) - firstRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (x3p2$header.info$incrementY*1e6)*dx,
                  bottomRightCorner_row = lastColCentered*sin((theta - median(theta))*(pi/180)) + firstRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (x3p2$header.info$incrementY*1e6)*dy,
                  bottomLeftCorner_col = firstColCentered*cos((theta - median(theta))*(pi/180)) - firstRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (x3p2$header.info$incrementY*1e6)*dx,
                  bottomLeftCorner_row = firstColCentered*sin((theta - median(theta))*(pi/180)) + firstRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (x3p2$header.info$incrementY*1e6)*dy) %>%
    #this is redundant, but is how the x and y columns are set-up down below, so
    #I won't change it
    dplyr::mutate(x_1 = topLeftCorner_col,
                  y_1 = topLeftCorner_row,
                  x_2 = topRightCorner_col,
                  y_2 = topRightCorner_row,
                  x_3 = bottomRightCorner_col,
                  y_3 = bottomRightCorner_row,
                  x_4 = bottomLeftCorner_col,
                  y_4 = bottomLeftCorner_row) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
                        names_to = c(".value","order"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(midCol = (topRightCorner_col + bottomLeftCorner_col)/2,
                  midRow = (topRightCorner_row + bottomLeftCorner_row)/2,
                  cellInd = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                                               floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                                               ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                                             nrow = ceiling(sqrt(max(cellNum))),
                                             byrow = TRUE),
                  x3p = rep(x3pNames[2],times = nrow(.)),
                  theta = theta - median(theta))

  x3p1_rotate <- 90 - median(allCells %>%
                               dplyr::filter(cmc != "non-CMC") %>%
                               dplyr::pull(theta))

  x3pPlt <- x3pListPlot(x3pList = list(x3p1,x3p2) %>%
                          setNames(x3pNames),
                        type = pltType,
                        rotate = c(ifelse(is.na(x3p1_rotate),90,x3p1_rotate),
                                   90),
                        legend.quantiles = legend.quantiles,
                        height.colors = height.colors,
                        na.value = na.value,
                        guide = "none")

  if(pltType == "faceted"){

    x3pPlt <- x3pPlt +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = x3p1_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellNum,
                                                   fill = cmc),
                            alpha = .3,
                            size = 2) +
      ggplot2::geom_polygon(data = x3p2_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellNum,
                                                   fill = cmc),
                            alpha = .3,
                            size = 2) +
      ggplot2::geom_text(data = dplyr::bind_rows(x3p1_cellGrid,
                                                 x3p2_cellGrid),
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellInd,
                                      colour = cmc,
                                      angle = theta),
                         size = 3) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type"))
  }
  else if(pltType == "list"){
    x3pPlt[[1]] <- x3pPlt[[1]] +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = x3p1_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellNum,
                                                   fill = cmc),
                            alpha = .3,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = x3p1_cellGrid,
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellInd,
                                      colour = cmc,
                                      angle = theta),
                         size = 3) +
      ggplot2::guides(colour = "legend")

    x3pPlt[[2]] <- x3pPlt[[2]] +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = x3p2_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellNum,
                                                   fill = cmc),
                            alpha = .3,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = x3p2_cellGrid,
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellInd,
                                      colour = cmc,
                                      angle = theta),
                         size = 3) +
      ggplot2::guides(colour = "legend")
  }

  return(x3pPlt)
}

#'Visualize initial and high CMCs for a cartridge case pair comparison
#'@name cmcPlot
#'
#'@description Constructs either a single faceted plot or a list of plots
#'  depicting the CMCs/non-CMCs under the initially proposed and High CMC
#'  methods for a pair of cartridge case scans
#'
#'@param x3p1 an x3p object
#'@param x3p2 a different x3p object
#'@param cellCCF_bothDirections_output output from the function
#'  cmcR::cellCCF_bothDirections
#'@param cmcFilter_improved_output output from the function
#'  cmcR::cmcFilter_improved
#'@param type argument to be passed to cmcR::x3pListPlot function
#'@param x3pNames (Optional) Names of x3p objects to be included in x3pListPlot
#'  function
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param cell.colors vector of 2 colors for plotting non-matching and matching
#'  (in that order) cells
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
#'
#' @examples
#' \dontrun{
#' comparison <- cellCCF_bothDirections(x3p1,x3p2)
#' cmcs <- cmcFilter_improved(comparison)
#'
#' cmcPlot(x3p1,x3p2,comparison,cmcs,type = "faceted",x3pNames = c("name1","name2"))
#'}
#'
#'@export

utils::globalVariables(c(".","cmc","comparison","dx","dy","theta","cellNum","cellID"))

cmcPlot <- function(x3p1,
                    x3p2,
                    cellCCF_bothDirections_output,
                    cmcFilter_improved_output,
                    type = "faceted",
                    x3pNames = c("x3p1","x3p2"),
                    legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                    height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                    cell.colors = c("#a50026","#313695"),
                    na.value = "gray80"){

  directionIndic <- which.min(c(nrow(cmcFilter_improved_output$initialCMC[[1]]),
                                nrow(cmcFilter_improved_output$initialCMC[[2]])))


  initialCMC <- cmcFilter_improved_output$initialCMC[[directionIndic]] %>%
    dplyr::mutate(cmc = rep("Top Vote CMC",times = nrow(.)))

  nonInitialCMC <- cellCCF_bothDirections_output[[directionIndic]]$ccfResults %>%
    topResultsPerCell() %>%
    dplyr::anti_join(initialCMC,by = "cellNum") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cmc = rep("non-CMC",times = nrow(.)))

  allInitialCells <- dplyr::bind_rows(initialCMC,nonInitialCMC) %>%
    dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","Top Vote CMC"))) %>%
    dplyr::left_join(dplyr::bind_rows(initialCMC,nonInitialCMC) %>%
                       purrr::pmap_dfr(~ {
                         idNum <- ..2 %>%
                           stringr::str_extract_all(string = ..2,
                                                    pattern = "[0-9]{1,}") %>%
                           unlist() %>%
                           as.numeric()

                         data.frame(cellID = ..2,
                                    firstRow = idNum[1],
                                    lastRow = idNum[2],
                                    firstCol = idNum[3],
                                    lastCol = idNum[4],
                                    stringsAsFactors = FALSE)
                       }),
                     by = "cellID")

  initialCMCPlt <- arrangeCMCPlot(x3p1 = list(x3p1,x3p2)[[directionIndic]],
                                  x3p2 = list(x3p2,x3p1)[[directionIndic]],
                                  allCells = allInitialCells,
                                  x3pNames = list(x3pNames,rev(x3pNames))[[directionIndic]],
                                  pltType = type,
                                  legend.quantiles = legend.quantiles,
                                  height.colors = height.colors,
                                  cell.colors = cell.colors,
                                  na.value = na.value)

  highCMCPlt <- NULL #missing by default unless high CMCs exist:

  if(!purrr::is_empty(cmcFilter_improved_output$highCMC)){

    x3p1_allCellIDs <- cellDivision(x3p1$surface.matrix,
                                    cellNumHoriz = ceiling(sqrt(max(cellCCF_bothDirections_output$comparison_1to2$ccfResults[[1]]$cellNum))),
                                    cellNumVert = ceiling(sqrt(max(cellCCF_bothDirections_output$comparison_1to2$ccfResults[[1]]$cellNum)))) %>%
      purrr::map2(names(.),
                  function(mat,horizCell){
                    purrr::map(.x = names(mat),
                               function(vertCell) swapCellIDAxes(paste(horizCell,vertCell,sep = ",")))
                  }) %>%
      unlist() %>%
      data.frame(cellID = .) %>%
      dplyr::mutate(cellNum = 1:nrow(.))

    highCMCs_directionCorrected <- cmcFilter_improved_output$highCMCs %>%
      dplyr::mutate(dx = ifelse(comparison == "comparison_2to1",
                                yes = -dx,
                                no = dx),
                    dy = ifelse(comparison == "comparison_2to1",
                                yes = -dy,
                                no = dy),
                    theta = ifelse(comparison == "comparison_2to1",
                                   yes = -theta,
                                   no = theta)) %>%
      dplyr::arrange(cellNum) %>%
      dplyr::select(-cellID) %>%
      dplyr::left_join(x3p1_allCellIDs,by = "cellNum") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cmc = rep("High CMC",times = nrow(.)),
                    cellID = as.character(cellID))

    nonHighCMCs_directionCorrected <- cellCCF_bothDirections_output$comparison_1to2$ccfResults %>%
      topResultsPerCell() %>%
      dplyr::anti_join(highCMCs_directionCorrected,
                       by = "cellNum") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cmc = rep("non-CMC",times = nrow(.)),
                    cellID = as.character(cellID))

    allHighCMCs <- dplyr::bind_rows(highCMCs_directionCorrected,nonHighCMCs_directionCorrected) %>%
      dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","High CMC"))) %>%
      dplyr::left_join(dplyr::bind_rows(highCMCs_directionCorrected,nonHighCMCs_directionCorrected) %>%
                         purrr::pmap_dfr(~ {
                           idNum <- ..8 %>%
                             stringr::str_extract_all(string = ..8,
                                                      pattern = "[0-9]{1,}") %>%
                             unlist() %>%
                             as.numeric()

                           data.frame(cellID = ..8,
                                      firstRow = idNum[1],
                                      lastRow = idNum[2],
                                      firstCol = idNum[3],
                                      lastCol = idNum[4],
                                      stringsAsFactors = FALSE)
                         }),
                       by = "cellID")

    highCMCPlt <- arrangeCMCPlot(x3p1 = x3p1,
                                 x3p2 = x3p2,
                                 allCells = allHighCMCs,
                                 x3pNames = x3pNames,
                                 pltType = type,
                                 legend.quantiles = legend.quantiles,
                                 height.colors = height.colors,
                                 cell.colors = cell.colors,
                                 na.value = na.value)
  }

  return(list("initialCMC" = initialCMCPlt,
              "highCMC" = highCMCPlt))
}

#' Create a bar plot of congruent matching cells per rotation value
#'
#' @name cmcPerThetaBarPlot
#'
#' @param cellCCF_output list returned by the cellCCF or cellCCF_bothDirections
#'   function. If from the cellCCF_bothDirections, then the ggplot will be
#'   faceted by the "direction" of the comparison (i.e., whether x3p1 or x3p2
#'   played the role as the "questioned" cartridge case scan)
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
#' @param highCMCThresh number of CMCs less than the max CMC count that should
#'   be classified as a "High CMC" count. That is, CMC_high = CMC_max -
#'   highCMCThresh
#' @param x3pNames (Optional) Names of x3p objects to be included in x3pListPlot
#'   function
#' @examples
#' \dontrun{
#'  #x3p1 and x3p2 are two x3p objects containing processed cartridge case scans
#'  comparison1 <- cellCCF(x3p1,x3p2)
#'
#'  #creates a single bar plot of the CMCs per theta
#'  cmcPerThetaBarPlot(comparison1)
#'
#'  comparison2 <- cellCCF_bothDirections(x3p1,x3p2)
#'
#'  #creates a faceted bar plot of the CMCs per theta in both comparison directions
#'  cmcPerThetaBarPlot(comparison2)
#' }
#'
#' @export
#'
#' @importFrom stats median

utils::globalVariables(c("comparison","theta","n","cmcHigh"))

cmcPerThetaBarPlot <- function(cellCCF_output,
                               consensus_function = median,
                               ccf_thresh = .6,
                               dx_thresh = 10,
                               dy_thresh = dx_thresh,
                               theta_thresh = 3,
                               consensus_function_theta = consensus_function,
                               highCMCThresh = 1,
                               x3pNames = c("x3p1","x3p2")){
  #make sure that cellCCF_output is either the output of the cellCCF function or
  #cellCCF_bothDirections
  testthat::expect_true(is.list(cellCCF_output))
  testthat::expect_true(identical(names(cellCCF_output),c("comparison_1to2","comparison_2to1")) | identical(names(cellCCF_output),c("params","ccfResults")),
                        label = "cellCCF_output argument is a list returned by cellCCF or cellCCF_bothDirections")

  #if cellCCF_output is from cellCCF_bothDirections:
  if(identical(names(cellCCF_output),c("comparison_1to2","comparison_2to1"))){
    cmcDat <- dplyr::bind_rows(
      cellCCF_output$comparison_1to2$ccfResults %>%
        cmcFilterPerTheta(consensus_function = consensus_function,
                          ccf_thresh = ccf_thresh,
                          dx_thresh = dx_thresh,
                          dy_thresh = dy_thresh,
                          theta_thresh = theta_thresh) %>%
        dplyr::mutate(comparison = paste0(x3pNames[1]," vs. ",x3pNames[2])),
      cellCCF_output$comparison_2to1$ccfResults %>%
        cmcFilterPerTheta(consensus_function = consensus_function,
                          ccf_thresh = ccf_thresh,
                          dx_thresh = dx_thresh,
                          dy_thresh = dy_thresh,
                          theta_thresh = theta_thresh) %>%
        dplyr::mutate(comparison = paste0(x3pNames[2]," vs. ",x3pNames[1]))
    ) %>%
      dplyr::group_by(comparison,theta) %>%
      dplyr::tally() %>%
      dplyr::ungroup()

    highCMCDat <- cmcDat %>%
      dplyr::group_by(comparison) %>%
      dplyr::summarise(x = min(theta) + 3,
                       y = max(n) + 2,
                       cmcHigh = max(n) - highCMCThresh)

    cmcDat %>%
      ggplot2::ggplot(ggplot2::aes(x = theta,y = n)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_bw()  +
      ggplot2::xlab(expression(theta*" (degree)")) +
      ggplot2::ylab("CMC number") +
      ggplot2::facet_wrap(~ comparison,ncol = 1)  +
      ggplot2::ylim(c(NA,25)) +
      ggplot2::geom_hline(data = highCMCDat,
                          ggplot2::aes(yintercept = cmcHigh),
                          colour = "black",
                          linetype = "dashed") +
      ggplot2::geom_text(data = highCMCDat,
                         ggplot2::aes(x = x,
                                      y = y,
                                      label = paste0("High CMC = ",cmcHigh)),
                         fontface = "plain",
                         family = "sans")
  }
  else{
    cellCCF_output$ccfResults  %>%
      cmcFilterPerTheta(consensus_function = consensus_function,
                        ccf_thresh = ccf_thresh,
                        dx_thresh = dx_thresh,
                        dy_thresh = dy_thresh,
                        theta_thresh = theta_thresh) %>%
      dplyr::group_by(theta) %>%
      dplyr::tally() %>%
      dplyr::mutate(cmcHigh = max(n) - highCMCThresh)  %>%
      ggplot2::ggplot(ggplot2::aes(x = theta,y = n)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_bw()  +
      ggplot2::xlab(expression(theta*" (degree)")) +
      ggplot2::ylab("CMC number") +
      ggplot2::ylim(c(NA,25)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = cmcHigh),
                          colour = "black",
                          linetype = "dashed") +
      ggplot2::geom_text(ggplot2::aes(x = min(theta) + 3,
                                      y = cmcHigh,
                                      label = paste0("High CMC = ",cmcHigh)),
                         nudge_y = 1,
                         fontface = "plain",
                         family = "sans")
  }

}

#'Extract cell/region pair matrices for a comparison of interest
#'@name getCellRegionPairs
#'
#'@description This function is meant as a diagnostic tool to determine what the
#'  cell/region pairs for a particular cartridge case comparison looked like
#'  right before calculation of the CCF.
#'
#'@param x3p1 an x3p object
#'@param x3p2 another x3p object
#'@param ccfDF data frame containing comparison results between x3p1 and x3p2. At a
#'  minimum, this needs to contain "cellNum" and "theta" columns.
#'@param cellCCF_params list of parameters under which the comparison represented in
#'  ccfDF was performed. Such a list is returned by the cellCCF and
#'  cellCCF_bothDirections function
#'
#'@examples
#'\dontrun{
#' comparison <- cellCCF(x3p1,x3p2)
#'
#' topResults <- comparison$ccfResults %>%
#'  cmcR::topResultsPerCell()
#'
#' cellRegionPairs <- getCellRegionPairs(x3p1,x3p2,topResults,comparison$params)
#'}
#'
#'@export
#'
#'@importFrom stats sd setNames

utils::globalVariables(c("theta"))

getCellRegionPairs <- function(x3p1,x3p2,ccfDF,cellCCF_params){
  mat1 <- x3p1$surface.matrix
  mat2 <- x3p2$surface.matrix

  if(is.null(cellCCF_params$centerCell)){
    m1 <- 0
    m2 <- 0
  }
  if(is.null(cellCCF_params$scaleCell)){
    sd1 <- 1
    sd2 <- 1
  }

  if(!is.null(cellCCF_params$centerCell)){
    if(cellCCF_params$centerCell == "wholeMatrix"){
      m1 <- mean(as.vector(mat1),na.rm = TRUE)

      m2 <- mean(as.vector(mat2),na.rm = TRUE)
    }
  }

  if(!is.null(cellCCF_params$scaleCell)){
    if(cellCCF_params$scaleCell == "wholeMatrix"){
      sd1 <- sd(as.vector(mat1),na.rm = TRUE)

      sd2 <- sd(as.vector(mat2),na.rm = TRUE)
    }
  }

  mat1_split <- splitSurfaceMat1(surfaceMat = mat1,
                                 cellNumHoriz = cellCCF_params$cellNumHoriz,
                                 cellNumVert = cellCCF_params$cellNumVert,
                                 minObservedProp = cellCCF_params$minObservedProp)

  sidelengthMultiplier <- floor(sqrt(cellCCF_params$regionToCellProp))

  mat2_splitCorners <- getMat2SplitLocations(cellIDs = mat1_split$cellIDs,
                                             cellSideLengths = mat1_split$cellSideLengths,
                                             mat2Dim = dim(mat2),
                                             sidelengthMultiplier = sidelengthMultiplier)

  thetas <- unique(ccfDF$theta)

  cellRegionPairs <- list()

  ccfDFSplit <- ccfDF %>%
    dplyr::group_by(theta) %>%
    dplyr::group_split()

  for(ind in 1:length(thetas)){
    mat2_rotated <- mat2 %>%
      rotateSurfaceMatrix(thetas[ind])

    mat2_splitRotated <-
      purrr::map(.x = mat2_splitCorners,
                 ~ extractCellbyCornerLocs(cornerLocs = .x,
                                           rotatedSurfaceMat = mat2_rotated,
                                           mat2Dim = dim(mat2)))

    mat1_splitFiltered <- purrr::flatten(mat1_split$surfaceMat_split)[ccfDFSplit[[ind]]$cellNum]
    mat2_splitFiltered <- mat2_splitRotated[ccfDFSplit[[ind]]$cellNum]

    filteredCellID <- mat1_split$cellIDs[ccfDFSplit[[ind]]$cellNum]

    if(!is.null(cellCCF_params$centerCell)){
      if(cellCCF_params$centerCell == "individualCell"){
        m1 <- mat1_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)



        m2 <-  mat2_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    if(!is.null(cellCCF_params$scaleCell)){
      if(cellCCF_params$scaleCell == "individualCell"){
        sd1 <-  mat1_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)

        sd2 <-  mat2_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    mat1_splitShifted <- purrr::pmap(list(mat1_splitFiltered,m1,sd1),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3))

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs).
    mat2_splitShifted <- purrr::pmap(list(mat2_splitFiltered,m2,sd2),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3))

    cellRegionPairs[[ind]] <- purrr::map2(mat1_splitShifted,
                                          mat2_splitShifted,
                                          ~ list(.x,.y)) %>%
      setNames(filteredCellID)
  }

  cellRegionPairs %>%
    setNames(thetas)
}

#' @name ccfMap
#'
#' @keywords internal

ccfMap <- function(mat1,mat2){
  #Test that mat1,mat2 contain no NAs
  if(any(is.na(as.vector(mat1))) | any(is.na(as.vector(mat2)))){
    mat1[is.na(mat1)] <- 0
    mat2[is.na(mat2)] <- 0
  }

  ccfMat <- filterViaFFT(mat1,mat2) / (sqrt(sum(mat1^2)) * sqrt(sum(mat2^2)))

  return(Re(ccfMat))
}

#' Plots the CCF map between two matrices
#'
#' @name ccfMapPlot
#'
#' @description Uses the gridExtra::arrangeGrob function to put 3 plots and a
#'   table on the same grob object.
#'
#' @param mat1 a matrix
#' @param mat2 another matrix
#' @param theta *(OPTIONAL)* if considering multiple cell/region pairs over
#'   various rotation values, it may be useful to include the rotation value in
#'   the fft.ccf summary. If a value of theta is supplied to this argument, it
#'   will be included in the summary table. Otherwise, theta will be NA
#' @param type dictates whether the CCF map is created using geom_raster (type =
#'   "raster") or geom_contour_filled (type = "contour")
#' @param returnGrob if TRUE, then function will return the gridExtra grob
#'   object
#'
#' @examples
#'  mat1 <- imager::imfill(x = 100,y = 100) %>%
#'  imager::draw_rect(x0 = 40,y0 = 40,x1 = 60,y1 = 60,color = 255) %>%
#'  as.matrix()
#'
#'  mat2 <- imager::imfill(x = 100,y = 100) %>%
#'  imager::draw_rect(x0 = 15,y0 = 30,x1 = 35,y1 = 50,color = 255) %>%
#'  as.matrix()
#'
#'  ccfMapPlot(mat1,mat2)
#'
#' @importFrom gridExtra tableGrob arrangeGrob grid.arrange
#' @importFrom colorspace divergingx_hcl
#' @importFrom scales rescale
#' @importFrom stats quantile
#'
#' @export

utils::globalVariables(c("x","y","value","fft.ccf","dx","dy"))

ccfMapPlot <- function(mat1,
                       mat2,
                       theta = NA,
                       returnGrob = FALSE,
                       type = "raster"){

  ccfMat <- ccfMap(mat1,mat2)

  ccfDF <- ccfMat %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(dx = x - max(x)/2,
                  dy = y - max(y)/2) %>%
    dplyr::rename(fft.ccf = value)

  ccfMaxInfo <- ccfDF %>%
    dplyr::filter(fft.ccf == max(fft.ccf)) %>%
    dplyr::mutate(fft.ccf = round(fft.ccf,3))

  # mat1 <- (mat1 - mean(mat1,na.rm = TRUE))/(sd(mat1,na.rm = TRUE))
  # mat2 <- (mat2 - mean(mat2,na.rm = TRUE))/(sd(mat2,na.rm = TRUE))

  mat1Plot <- mat1 %>%
    t() %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    # mutate(x = x + floor(max(nrow(mat1),nrow(mat2))/4),
    #        y = y + floor(max(ncol(mat1),ncol(mat2))/4)) %>%
    ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient2(low = "grey0",
                                  mid = "grey50",
                                  high = "grey100") +
    ggplot2::coord_fixed(xlim = c(0,max(nrow(mat1),nrow(mat2))),
                         ylim = c(0,max(ncol(mat1),ncol(mat2))),
                         expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())

  mat2Plot <- mat2 %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    imager::as.cimg() %>%
    as.data.frame() %>%
    ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient2(low = "grey0",
                                  mid = "grey50",
                                  high = "grey100",
                                  midpoint = 0) +
    ggplot2::coord_fixed(xlim = c(0,max(nrow(mat1),nrow(mat2))),
                         ylim = c(0,max(ncol(mat1),ncol(mat2))),
                         expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::geom_tile(ggplot2::aes(
      x = max(x)/2 - ccfMaxInfo$dx,
      y = max(y)/2 - ccfMaxInfo$dy,
      width = ncol(mat1),
      height = nrow(mat1)),
      alpha = 0,
      colour = "orange")

  if(type == "raster"){
    ccfPlot <- ccfDF %>%
      dplyr::mutate(dx = rev(dx),
                    dy = rev(dy)) %>%
      ggplot2::ggplot(ggplot2::aes(x = dx,y = dy)) +
      ggplot2::geom_raster(ggplot2::aes(fill = fft.ccf)) +
      ggplot2::geom_point(data = ccfDF %>%
                            dplyr::mutate(dx = rev(dx),
                                          dy = rev(dy)) %>%
                            dplyr::filter(fft.ccf == max(fft.ccf)) %>%
                            dplyr::mutate(type = "Maximum"),
                          ggplot2::aes(x = dx,y = dy,shape = type),
                          # shape = 4,
                          colour = "white",
                          fill = "white") +
      # ggplot2::geom_contour(aes(z = fft.ccf),
      #                       breaks = quantile(ccfDF$fft.ccf,seq(0,1,length.out = 5)),
      # colour = "black") +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                                    values = rescale(c(min(ccfDF$fft.ccf),0,max(ccfDF$fft.ccf))),
                                    guide = "colourbar",
                                    limits = c(min(ccfDF$fft.ccf),max(ccfDF$fft.ccf))) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_shape_manual(values = 1,
                                  labels = expression("CCF"[max])) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(title = "CCF",
                                                      barheight = grid::unit(1,"in"),
                                                      label.theme = ggplot2::element_text(size = 8),
                                                      frame.colour = "black",
                                                      ticks.colour = "black"),
                      shape = ggplot2::guide_legend(title = NULL,
                                                    override.aes = list(colour = "black")))

    layoutMat <- matrix(c(1,2,4,3),ncol = 2,byrow = TRUE)
  }
  if(type == "contour"){
    ccfPlot <- ccfDF %>%
      dplyr::mutate(dx = rev(dx),
                    dy = rev(dy)) %>%
      ggplot2::ggplot(ggplot2::aes(x = dx,y = dy)) +
      ggplot2::geom_contour_filled(ggplot2::aes(z = fft.ccf),
                                   breaks = c(quantile(as.vector(ccfDF$fft.ccf)[(as.vector(ccfDF$fft.ccf) <= 0)],
                                                       prob = seq(0,1,length.out = 6)),
                                              0,
                                              quantile(as.vector(ccfDF$fft.ccf)[as.vector((ccfDF$fft.ccf) > 0)],
                                                       prob = seq(0,1,length.out = 6))),
      ) +
      ggplot2::geom_point(data = ccfDF %>%
                            dplyr::mutate(dx = rev(dx),
                                          dy = rev(dy)) %>%
                            dplyr::filter(fft.ccf == min(fft.ccf) | fft.ccf == max(fft.ccf)) %>%
                            dplyr::arrange(fft.ccf) %>%
                            dplyr::mutate(type = factor(c("Minimum","Maximum"))),
                          ggplot2::aes(x = dx,y = dy,shape = type),
                          # shape = 4,
                          colour = "white") +
      ggplot2::scale_shape_manual(values = c(1,4)) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::scale_fill_manual(values = divergingx_hcl(n = 13,
                                                         palette = "PuOr",
                                                         rev = TRUE),
                                 drop = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 7),
                     legend.text = ggplot2::element_text(size = 5)) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE,ncol = 3),
                      shape = ggplot2::guide_legend(direction = "horizontal",override.aes = list(colour = c("#312A56","#492900"))))

    layoutMat <- matrix(c(1,2,3,4,4,4),ncol = 3,byrow = TRUE)
  }

  ccfMaxSummary <- ccfMaxInfo %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::mutate(dx = -dx,
                  dy = -dy,
                  theta = theta) %>%
    t() %>%
    tableGrob(rows = c(expression("CCF"[max]),"dx","dy",expression(theta)),
              cols = "Summary")

  gridPlot <- arrangeGrob(mat1Plot,
                          mat2Plot,
                          ccfMaxSummary,
                          ccfPlot,
                          layout_matrix = layoutMat)

  grid.arrange(gridPlot)

  if(returnGrob){
    return(gridPlot)
  }
}