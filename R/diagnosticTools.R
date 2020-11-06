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
#' data(fadul1.1_processed,fadul1.2_processed,fadul2.1_processed)
#'
#' x3pListPlot(list("Fadul 1-1" = fadul1.1_processed,
#'                  "Fadul 1-2" = fadul1.2_processed,
#'                  "Fadul 2-1" = fadul2.1_processed))
#'
#' @export
#'
#' @importFrom stats setNames median quantile
#' @importFrom rlang .data

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
                                         dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                                         dplyr::mutate(x = .data$x*1e6,
                                                       y = .data$y*1e6,
                                                       height = .data$value*1e6) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(.data$x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$height))  +
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
                      colour = 'none') +
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
                            dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::mutate(x = .data$x*1e6,
                                          y = .data$y*1e6,
                                          height = .data$value*1e6) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = .data$height))  +
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
                                           panel.background = ggplot2::element_blank(),
                                           plot.title = ggplot2::element_text(hjust = .5,
                                                                              size = 11)) +
                            ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(3,"in"),
                                                                            label.theme = ggplot2::element_text(size = 8),
                                                                            title.theme = ggplot2::element_text(size = 10),
                                                                            frame.colour = "black",
                                                                            ticks.colour = "black"),
                                            colour =  'none') +
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
#' @importFrom rlang .data

linear_to_matrix <- function(index, nrow = 7, ncol = nrow, byrow = TRUE, sep = ", ") {
  index <- as.integer(index)
  stopifnot(all(index <= nrow*ncol), all(index > 0))

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
#' @importFrom rlang .data

arrangeCMCPlot <- function(referenceScan,
                           targetScan,
                           allCells,
                           x3pNames,
                           pltType = "faceted",
                           legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                           height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           cell.colors = c("#a50026","#313695"),
                           cell.alpha = .2,
                           na.value = "gray80"){

  targetScan_cellGrid <- allCells %>%
    dplyr::mutate(firstRow = (referenceScan$header.info$incrementY*1e6)*(.data$firstRow),
                  lastRow = (referenceScan$header.info$incrementY*1e6)*(.data$lastRow),
                  firstCol = (referenceScan$header.info$incrementY*1e6)*(.data$firstCol),
                  lastCol = (referenceScan$header.info$incrementY*1e6)*(.data$lastCol)) %>%
    dplyr::mutate(x_1 = .data$firstCol,
                  y_1 = .data$firstRow,
                  x_2 = .data$lastCol,
                  y_2 = .data$firstRow,
                  x_3 = .data$lastCol,
                  y_3 = .data$lastRow,
                  x_4 = .data$firstCol,
                  y_4 = .data$lastRow) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
                        names_to = c(".value","order"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(midCol = (.data$lastCol + .data$firstCol)/2,
                  midRow = (.data$lastRow + .data$firstRow)/2,
                  # cellIndex = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                  #                              floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                  #                              ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                  #                            nrow = ceiling(sqrt(max(cellNum))),
                  # byrow = TRUE),
                  x3p = rep(x3pNames[1],times = nrow(.)),
                  theta = rep(0,times = nrow(.)))

  referenceScan_cellGrid <- allCells %>%
    dplyr::mutate(firstRow = (targetScan$header.info$incrementY*1e6)*(.data$firstRow),
                  lastRow = (targetScan$header.info$incrementY*1e6)*(.data$lastRow),
                  firstCol = (targetScan$header.info$incrementY*1e6)*(.data$firstCol),
                  lastCol = (targetScan$header.info$incrementY*1e6)*(.data$lastCol)) %>%
    dplyr::mutate(firstRowCentered = .data$firstRow - max(.data$lastRow)/2,
                  lastRowCentered = .data$lastRow - max(.data$lastRow)/2,
                  firstColCentered = .data$firstCol - max(.data$lastCol)/2,
                  lastColCentered = .data$lastCol - max(.data$lastCol)/2) %>%
    dplyr::mutate(topLeftCorner_col = .data$firstColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$lastRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (targetScan$header.info$incrementY*1e6)*.data$x/2,
                  topLeftCorner_row = .data$firstColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$lastRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (targetScan$header.info$incrementY*1e6)*.data$y/2,
                  topRightCorner_col = .data$lastColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$lastRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (targetScan$header.info$incrementY*1e6)*.data$x/2,
                  topRightCorner_row = .data$lastColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$lastRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (targetScan$header.info$incrementY*1e6)*.data$y/2,
                  bottomRightCorner_col = .data$lastColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$firstRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (targetScan$header.info$incrementY*1e6)*.data$x/2,
                  bottomRightCorner_row = .data$lastColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$firstRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (targetScan$header.info$incrementY*1e6)*.data$y/2,
                  bottomLeftCorner_col = .data$firstColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$firstRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (targetScan$header.info$incrementY*1e6)*.data$x/2,
                  bottomLeftCorner_row = .data$firstColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$firstRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (targetScan$header.info$incrementY*1e6)*.data$y/2) %>%
    #this is redundant, but are the names attributed to the x and y columns are
    #set-up down below, so I won't change it
    dplyr::mutate(x_1 = .data$topLeftCorner_col,
                  y_1 = .data$topLeftCorner_row,
                  x_2 = .data$topRightCorner_col,
                  y_2 = .data$topRightCorner_row,
                  x_3 = .data$bottomRightCorner_col,
                  y_3 = .data$bottomRightCorner_row,
                  x_4 = .data$bottomLeftCorner_col,
                  y_4 = .data$bottomLeftCorner_row) %>%
    tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
                        names_to = c(".value","order"),
                        names_pattern = "(.+)_(.+)") %>%
    dplyr::mutate(midCol = (.data$topRightCorner_col + .data$bottomLeftCorner_col)/2,
                  midRow = (.data$topRightCorner_row + .data$bottomLeftCorner_row)/2,
                  # cellIndex = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                  #                              floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                  #                              ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                  #                            nrow = ceiling(sqrt(max(cellNum))),
                  # byrow = TRUE),
                  x3p = rep(x3pNames[2],times = nrow(.)),
                  theta = .data$theta - median(.data$theta))

  targetScan_rotate <- 90 #- abs(median(allCells %>%
  # dplyr::filter(cmc != "non-CMC") %>%
  # dplyr::pull(theta)))

  x3pPlt <- x3pListPlot(x3pList = list(referenceScan,targetScan) %>%
                          setNames(x3pNames),
                        type = pltType,
                        rotate = c(90,
                                   ifelse(is.na(targetScan_rotate),90,targetScan_rotate)),
                        legend.quantiles = legend.quantiles,
                        height.colors = height.colors,
                        na.value = na.value,
                        guide = "none")

  if(pltType == "faceted"){

    x3pPlt <- x3pPlt +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = targetScan_cellGrid,
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   group = .data$cellIndex,
                                                   fill = .data$cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::geom_polygon(data = referenceScan_cellGrid,
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   group = .data$cellIndex,
                                                   fill = .data$cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::geom_text(data = dplyr::bind_rows(targetScan_cellGrid,
                                                 referenceScan_cellGrid),
                         ggplot2::aes(x = .data$midCol,
                                      y = .data$midRow,
                                      label = .data$cellIndex,
                                      colour = .data$cmc,
                                      angle = .data$theta),
                         size = 3) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
      ggplot2::theme(legend.position = "bottom")
  }
  else if(pltType == "list"){
    x3pPlt[[1]] <- x3pPlt[[1]] +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = targetScan_cellGrid,
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   group = .data$cellIndex,
                                                   fill = .data$cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = targetScan_cellGrid,
                         ggplot2::aes(x = .data$midCol,
                                      y = .data$midRow,
                                      label = .data$cellIndex,
                                      colour = .data$cmc,
                                      angle = .data$theta),
                         size = 3) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
      ggplot2::theme(legend.position = "bottom")

    x3pPlt[[2]] <- x3pPlt[[2]] +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = referenceScan_cellGrid,
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   group = .data$cellIndex,
                                                   fill = .data$cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = referenceScan_cellGrid,
                         ggplot2::aes(x = .data$midCol,
                                      y = .data$midRow,
                                      label = .data$cellIndex,
                                      colour = .data$cmc,
                                      angle = .data$theta),
                         size = 3) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
      ggplot2::theme(legend.position = "bottom")
  }

  return(x3pPlt)
}

#'Visualize initial and high CMCs for a cartridge case pair comparison
#'@name cmcPlot
#'
#'@description Constructs either a single faceted plot or a list of plots
#'  depicting the CMCs/non-CMCs under the initially proposed and High CMC
#'  methods for a pair of cartridge case scans
#'@param referenceScan an x3p object
#'@param targetScan a different x3p object
#'@param reference_v_target_CMCs CMCs for the comparison between the reference
#'  scan and the target scan.
#'@param target_v_reference_CMCs (optional) CMCs for the comparison between the
#'  target scan and the reference scan. If this is missing, then only the
#'  original method CMCs will be plotted
#' @param corColName name of correlation similarity score column used to
#'   identify the CMCs in the two comparison_*_df data frames (e.g.,
#'   pairwiseCompCor)
#'@param type argument to be passed to cmcR::x3pListPlot function
#'@param x3pNames (Optional) Names of x3p objects to be included in x3pListPlot
#'  function
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param cell.colors vector of 2 colors for plotting non-matching and matching
#'  (in that order) cells
#'@param cell.alpha sets alpha of cells (passed to geom_polygon)
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
#'
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
#'cmcPlot(fadul1.1_processed,
#'        fadul1.2_processed,
#'        comparisonDF_1to2,
#'        comparisonDF_2to1,
#'        corColName = "pairwiseCompCor")
#'
#'@importFrom utils hasName
#'@importFrom rlang .data
#'@export

cmcPlot <- function(referenceScan,
                    targetScan,
                    reference_v_target_CMCs,
                    target_v_reference_CMCs = reference_v_target_CMCs,
                    corColName = "pairwiseCompCor",
                    type = "faceted",
                    x3pNames = c("referenceScan","targetScan"),
                    legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                    height.colors = c("#1B1B1B","#404040","#7B7B7B","#B0B0B0","#DBDBDB","#F7F7F7","#E4E4E4","#C5C5C5","#999999","#717171","#4E4E4E"),
                    cell.colors = c("#a50026","#313695"),
                    cell.alpha = .2,
                    na.value = "gray80"){

  referenceScan_cellCorners <- referenceScan %>%
    comparison_cellDivision() %>%
    purrr::pmap_dfr(~ {
      idNum <- ..2$cmcR.info$cellRange %>%
        stringr::str_extract_all(string = ..2$cmcR.info$cellRange,
                                 pattern = "[0-9]{1,}") %>%
        unlist() %>%
        as.numeric()

      data.frame(cellIndex = ..1,
                 firstRow = idNum[1],
                 lastRow = idNum[2],
                 firstCol = idNum[3],
                 lastCol = idNum[4],
                 stringsAsFactors = FALSE)
    })

  targetScan_cellCorners <- targetScan %>%
    comparison_cellDivision() %>%
    purrr::pmap_dfr(~ {
      idNum <- ..2$cmcR.info$cellRange %>%
        stringr::str_extract_all(string = ..2$cmcR.info$cellRange,
                                 pattern = "[0-9]{1,}") %>%
        unlist() %>%
        as.numeric()

      data.frame(cellIndex = ..1,
                 firstRow = idNum[1],
                 lastRow = idNum[2],
                 firstCol = idNum[3],
                 lastCol = idNum[4],
                 stringsAsFactors = FALSE)
    })

  if(nrow(reference_v_target_CMCs) == 0){

    originalMethodCMCs <- as.data.frame(matrix(nrow = 0,ncol = ncol(reference_v_target_CMCs))) %>%
      setNames(names(reference_v_target_CMCs))

  }
  else{
    originalMethodCMCs <- reference_v_target_CMCs %>%
      dplyr::filter(.data$originalMethodClassif == "CMC")
  }

  nonoriginalMethodCMCs <- reference_v_target_CMCs %>%
    dplyr::filter(!(.data$cellIndex %in% originalMethodCMCs$cellIndex))

  if(nrow(nonoriginalMethodCMCs) > 0){
    nonoriginalMethodCMCs <- nonoriginalMethodCMCs %>%
      dplyr::group_by(.data$cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
  }

  allInitialCells_reference_v_target <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
    dplyr::mutate(cmc = ifelse(.data$originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","Original Method CMC"))) %>%
    dplyr::left_join(referenceScan_cellCorners,
                     by = "cellIndex")

  originalMethodCMCsPlt_reference_v_target <- arrangeCMCPlot(referenceScan = referenceScan,
                                                             targetScan = targetScan,
                                                             allCells = allInitialCells_reference_v_target,
                                                             x3pNames = x3pNames,
                                                             pltType = type,
                                                             legend.quantiles = legend.quantiles,
                                                             height.colors = height.colors,
                                                             cell.colors = cell.colors,
                                                             cell.alpha = cell.alpha,
                                                             na.value = na.value)

  #If only data for one comparison direction were given, only plot the original
  #method CMCs in that direction
  if(assertthat::are_equal(reference_v_target_CMCs, target_v_reference_CMCs)){

    return(list("CMCs" = originalMethodCMCsPlt_reference_v_target))

  }

  #otherwise, create the original CMC plots for the other direction and return:

  if(nrow(target_v_reference_CMCs) == 0){

    originalMethodCMCs <- as.data.frame(matrix(nrow = 0,ncol = ncol(target_v_reference_CMCs))) %>%
      setNames(names(target_v_reference_CMCs))

  }
  else{
    originalMethodCMCs <- target_v_reference_CMCs %>%
      dplyr::filter(.data$originalMethodClassif == "CMC")
  }

  nonoriginalMethodCMCs <- target_v_reference_CMCs %>%
    dplyr::filter(!(.data$cellIndex %in% originalMethodCMCs$cellIndex)) %>%
    dplyr::group_by(.data$cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  allInitialCells_target_v_reference <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
    dplyr::mutate(cmc = ifelse(.data$originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","Original Method CMC"))) %>%
    dplyr::left_join(referenceScan_cellCorners,
                     by = "cellIndex")

  originalMethodCMCsPlt_target_v_reference <- arrangeCMCPlot(referenceScan = targetScan,
                                                             targetScan = referenceScan,
                                                             allCells = allInitialCells_target_v_reference,
                                                             x3pNames = rev(x3pNames),
                                                             pltType = type,
                                                             legend.quantiles = legend.quantiles,
                                                             height.colors = height.colors,
                                                             cell.colors = cell.colors,
                                                             cell.alpha = cell.alpha,
                                                             na.value = na.value)

    if(!hasName(reference_v_target_CMCs,"highCMCClassif") | !hasName(target_v_reference_CMCs,"highCMCClassif")){

      return(list("originalMethodCMCs_reference_v_target" = originalMethodCMCsPlt_reference_v_target,
                  "originalMethodCMCs_target_v_reference" = originalMethodCMCsPlt_target_v_reference))

    }

  #If the necessary data to construct the High CMCs were given, then plot them
  #too.

  highCMCs_reference_v_target <- reference_v_target_CMCs %>%
    dplyr::filter(.data$highCMCClassif == "CMC")

  #Remaining cells not identified as High CMCs
  non_highCMCs_reference_v_target <- reference_v_target_CMCs %>%
    dplyr::filter(!(.data$cellIndex %in% highCMCs_reference_v_target$cellIndex)) %>%
    dplyr::group_by(.data$cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  highCMC_plotData_reference_v_target <- dplyr::bind_rows(highCMCs_reference_v_target,
                                                          non_highCMCs_reference_v_target) %>%
    dplyr::mutate(cmc = ifelse(.data$highCMCClassif == "CMC","High CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","High CMC"))) %>%
    dplyr::left_join(referenceScan_cellCorners,
                     by = "cellIndex")

  highCMCPlt_reference_v_target <- arrangeCMCPlot(referenceScan = referenceScan,
                                                  targetScan = targetScan,
                                                  allCells = highCMC_plotData_reference_v_target,
                                                  x3pNames = x3pNames,
                                                  pltType = type,
                                                  legend.quantiles = legend.quantiles,
                                                  height.colors = height.colors,
                                                  cell.colors = cell.colors,
                                                  cell.alpha = cell.alpha,
                                                  na.value = na.value)

  #Different High CMCs may have been identified in the other direction -- we
  #need to plot those separately

  highCMCs_target_v_reference <- target_v_reference_CMCs %>%
    dplyr::filter(.data$highCMCClassif == "CMC")

  #Remaining cells not identified as High CMCs
  non_highCMCs_target_v_reference <- target_v_reference_CMCs %>%
    dplyr::filter(!(.data$cellIndex %in% highCMCs_target_v_reference$cellIndex)) %>%
    dplyr::group_by(.data$cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  highCMC_plotData_target_v_reference <- dplyr::bind_rows(highCMCs_target_v_reference,
                                                          non_highCMCs_target_v_reference) %>%
    dplyr::mutate(cmc = ifelse(.data$highCMCClassif == "CMC","High CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","High CMC"))) %>%
    dplyr::left_join(targetScan_cellCorners,
                     by = "cellIndex")

  highCMCPlt_target_v_reference <- arrangeCMCPlot(referenceScan = targetScan,
                                                  targetScan = referenceScan,
                                                  allCells = highCMC_plotData_target_v_reference,
                                                  x3pNames = rev(x3pNames),
                                                  pltType = type,
                                                  legend.quantiles = legend.quantiles,
                                                  height.colors = height.colors,
                                                  cell.colors = cell.colors,
                                                  cell.alpha = cell.alpha,
                                                  na.value = na.value)


  return(list("originalMethodCMCs_reference_v_target" = originalMethodCMCsPlt_reference_v_target,
              "originalMethodCMCs_target_v_reference" = originalMethodCMCsPlt_target_v_reference,
              "highCMC_reference_v_target"= highCMCPlt_reference_v_target,
              "highCMC_target_v_reference"= highCMCPlt_target_v_reference))

}
