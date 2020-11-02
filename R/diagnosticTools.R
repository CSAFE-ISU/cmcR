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

utils::globalVariables(c("firstRow","lastRow","firstCol","lastCol","cellNum",".","firstColCentered","theta","firstRowCentered","x","y","lastColCentered","lastRowCentered","topRightCorner_col","bottomLeftCorner_col","topRightCorner_row","bottomLeftCorner_row","cmc","midCol","midRow","cellIndex"))

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
    dplyr::mutate(firstRow = (referenceScan$header.info$incrementY*1e6)*(firstRow),
                  lastRow = (referenceScan$header.info$incrementY*1e6)*(lastRow),
                  firstCol = (referenceScan$header.info$incrementY*1e6)*(firstCol),
                  lastCol = (referenceScan$header.info$incrementY*1e6)*(lastCol)) %>%
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
                  # cellIndex = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                  #                              floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                  #                              ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                  #                            nrow = ceiling(sqrt(max(cellNum))),
                  # byrow = TRUE),
                  x3p = rep(x3pNames[1],times = nrow(.)),
                  theta = rep(0,times = nrow(.)))

  referenceScan_cellGrid <- allCells %>%
    dplyr::mutate(firstRow = (targetScan$header.info$incrementY*1e6)*(firstRow),
                  lastRow = (targetScan$header.info$incrementY*1e6)*(lastRow),
                  firstCol = (targetScan$header.info$incrementY*1e6)*(firstCol),
                  lastCol = (targetScan$header.info$incrementY*1e6)*(lastCol)) %>%
    dplyr::mutate(firstRowCentered = firstRow - max(lastRow)/2,
                  lastRowCentered = lastRow - max(lastRow)/2,
                  firstColCentered = firstCol - max(lastCol)/2,
                  lastColCentered = lastCol - max(lastCol)/2) %>%
    dplyr::mutate(topLeftCorner_col = firstColCentered*cos((theta - median(theta))*(pi/180)) - lastRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (targetScan$header.info$incrementY*1e6)*x/2,
                  topLeftCorner_row = firstColCentered*sin((theta - median(theta))*(pi/180)) + lastRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (targetScan$header.info$incrementY*1e6)*y/2,
                  topRightCorner_col = lastColCentered*cos((theta - median(theta))*(pi/180)) - lastRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (targetScan$header.info$incrementY*1e6)*x/2,
                  topRightCorner_row = lastColCentered*sin((theta - median(theta))*(pi/180)) + lastRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (targetScan$header.info$incrementY*1e6)*y/2,
                  bottomRightCorner_col = lastColCentered*cos((theta - median(theta))*(pi/180)) - firstRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (targetScan$header.info$incrementY*1e6)*x/2,
                  bottomRightCorner_row = lastColCentered*sin((theta - median(theta))*(pi/180)) + firstRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (targetScan$header.info$incrementY*1e6)*y/2,
                  bottomLeftCorner_col = firstColCentered*cos((theta - median(theta))*(pi/180)) - firstRowCentered*sin((theta - median(theta))*(pi/180)) + max(lastCol)/2 - (targetScan$header.info$incrementY*1e6)*x/2,
                  bottomLeftCorner_row = firstColCentered*sin((theta - median(theta))*(pi/180)) + firstRowCentered*cos((theta - median(theta))*(pi/180)) + max(lastRow)/2 - (targetScan$header.info$incrementY*1e6)*y/2) %>%
    #this is redundant, but are the names attributed to the x and y columns are
    #set-up down below, so I won't change it
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
                  # cellIndex = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                  #                              floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                  #                              ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                  #                            nrow = ceiling(sqrt(max(cellNum))),
                  # byrow = TRUE),
                  x3p = rep(x3pNames[2],times = nrow(.)),
                  theta = theta - median(theta))

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
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellIndex,
                                                   fill = cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::geom_polygon(data = referenceScan_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellIndex,
                                                   fill = cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::geom_text(data = dplyr::bind_rows(targetScan_cellGrid,
                                                 referenceScan_cellGrid),
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellIndex,
                                      colour = cmc,
                                      angle = theta),
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
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellIndex,
                                                   fill = cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = targetScan_cellGrid,
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellIndex,
                                      colour = cmc,
                                      angle = theta),
                         size = 3) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
      ggplot2::theme(legend.position = "bottom")

    x3pPlt[[2]] <- x3pPlt[[2]] +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_polygon(data = referenceScan_cellGrid,
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   group = cellIndex,
                                                   fill = cmc),
                            alpha = cell.alpha,
                            size = 2) +
      ggplot2::scale_colour_manual(values = cell.colors,
                                   aesthetics = c("fill","colour")) +
      ggplot2::geom_text(data = referenceScan_cellGrid,
                         ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = cellIndex,
                                      colour = cmc,
                                      angle = theta),
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
#'@export

utils::globalVariables(c(".","cmc","comparison","x","y","theta","cellIndex","cellRange"))

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
      dplyr::filter(originalMethodClassif == "CMC")
  }

  nonoriginalMethodCMCs <- reference_v_target_CMCs %>%
    dplyr::filter(!(cellIndex %in% originalMethodCMCs$cellIndex))

  if(nrow(nonoriginalMethodCMCs) > 0){
    nonoriginalMethodCMCs <- nonoriginalMethodCMCs %>%
      dplyr::group_by(cellIndex) %>%
      dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
  }

  allInitialCells_reference_v_target <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
    dplyr::mutate(cmc = ifelse(originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","Original Method CMC"))) %>%
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
      dplyr::filter(originalMethodClassif == "CMC")
  }

  nonoriginalMethodCMCs <- target_v_reference_CMCs %>%
    dplyr::filter(!(cellIndex %in% originalMethodCMCs$cellIndex)) %>%
    dplyr::group_by(cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  allInitialCells_target_v_reference <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
    dplyr::mutate(cmc = ifelse(originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","Original Method CMC"))) %>%
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

  originalMethod_bothDirections <-

    if(!hasName(reference_v_target_CMCs,"highCMCClassif") | !hasName(target_v_reference_CMCs,"highCMCClassif")){

      return(list("originalMethodCMCs_reference_v_target" = originalMethodCMCsPlt_reference_v_target,
                  "originalMethodCMCs_target_v_reference" = originalMethodCMCsPlt_target_v_reference))

    }

  #If the necessary data to construct the High CMCs were given, then plot them
  #too.

  highCMCs_reference_v_target <- reference_v_target_CMCs %>%
    dplyr::filter(highCMCClassif == "CMC")

  #Remaining cells not identified as High CMCs
  non_highCMCs_reference_v_target <- reference_v_target_CMCs %>%
    dplyr::filter(!(cellIndex %in% highCMCs_reference_v_target$cellIndex)) %>%
    dplyr::group_by(cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  highCMC_plotData_reference_v_target <- dplyr::bind_rows(highCMCs_reference_v_target,
                                                          non_highCMCs_reference_v_target) %>%
    dplyr::mutate(cmc = ifelse(highCMCClassif == "CMC","High CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","High CMC"))) %>%
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
    dplyr::filter(highCMCClassif == "CMC")

  #Remaining cells not identified as High CMCs
  non_highCMCs_target_v_reference <- target_v_reference_CMCs %>%
    dplyr::filter(!(cellIndex %in% highCMCs_target_v_reference$cellIndex)) %>%
    dplyr::group_by(cellIndex) %>%
    dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))

  highCMC_plotData_target_v_reference <- dplyr::bind_rows(highCMCs_target_v_reference,
                                                          non_highCMCs_target_v_reference) %>%
    dplyr::mutate(cmc = ifelse(highCMCClassif == "CMC","High CMC","non-CMC")) %>%
    dplyr::mutate(cmc = factor(cmc,levels = c("non-CMC","High CMC"))) %>%
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
#'   table on the same grob object. The table contains summary statistics
#'   including the rotation value used (if specified in the theta argument), the
#'   estimated maximum CCF value calculated using the Cross-Correlation theorem,
#'   and the estimated optimal translation values based on the location of the
#'   maximum CCF value.
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
#' @importFrom scales rescale
#' @importFrom stats quantile
#'
#' @keywords internal

utils::globalVariables(c("x","y","value","fft.ccf","x","y"))

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
    dplyr::mutate(x = x - max(x)/2,
                  y = y - max(y)/2) %>%
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
      x = max(x)/2 - ccfMaxInfo$x,
      y = max(y)/2 - ccfMaxInfo$y,
      width = ncol(mat1),
      height = nrow(mat1)),
      alpha = 0,
      colour = "orange")

  if(type == "raster"){
    ccfPlot <- ccfDF %>%
      dplyr::mutate(x = rev(x),
                    y = rev(y)) %>%
      ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = fft.ccf)) +
      ggplot2::geom_point(data = ccfDF %>%
                            dplyr::mutate(x = rev(x),
                                          y = rev(y)) %>%
                            dplyr::filter(fft.ccf == max(fft.ccf)) %>%
                            dplyr::mutate(type = "Maximum"),
                          ggplot2::aes(x = x,y = y,shape = type),
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
      dplyr::mutate(x = rev(x),
                    y = rev(y)) %>%
      ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_contour_filled(ggplot2::aes(z = fft.ccf),
                                   breaks = c(quantile(as.vector(ccfDF$fft.ccf)[(as.vector(ccfDF$fft.ccf) <= 0)],
                                                       prob = seq(0,1,length.out = 6)),
                                              0,
                                              quantile(as.vector(ccfDF$fft.ccf)[as.vector((ccfDF$fft.ccf) > 0)],
                                                       prob = seq(0,1,length.out = 6))),
      ) +
      ggplot2::geom_point(data = ccfDF %>%
                            dplyr::mutate(x = rev(x),
                                          y = rev(y)) %>%
                            dplyr::filter(fft.ccf == min(fft.ccf) | fft.ccf == max(fft.ccf)) %>%
                            dplyr::arrange(fft.ccf) %>%
                            dplyr::mutate(type = factor(c("Minimum","Maximum"))),
                          ggplot2::aes(x = x,y = y,shape = type),
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
    dplyr::mutate(x = -x,
                  y = -y,
                  theta = theta) %>%
    t() %>%
    tableGrob(rows = c(expression("FFT CCF"[max]),"x","y",expression(theta)),
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