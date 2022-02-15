#'Plot a list of x3ps
#'@name x3pListPlot
#'
#'@description Plots the surface matrices in a list of x3p objects. Either
#'  creates one plot faceted by surface matrix or creates individual plots per
#'  surface matrix and returns them in a list.
#'
#'@param x3pList a list of x3p objects. If the x3p objects are named in the
#'  list, then these names will be included in the title of their respective
#'  plot
#'@param type dictates whether one plot faceted by surface matrix or a list of
#'  plots per surface matrix is returned. The faceted plot will have a
#'  consistent height scale across all surface matrices.
#'@param rotate angle (in degrees) to rotate all surface matrices plotted
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
#'@param guide internal usage
#'@return A ggplot object or list of ggplot objects showing the surface matrix
#'  height values.
#' @examples
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#' x3pListPlot(list("Fadul 1-1" = fadul1.1_processed,
#'                  "Fadul 1-2" = fadul1.2_processed))
#'@export
#'
#'@importFrom stats setNames median quantile
#'@importFrom rlang .data

x3pListPlot <- function(x3pList,
                        type = "faceted",
                        legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                        height.quantiles = c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                        height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                        na.value = "gray80"){
  if(purrr::is_empty(names(x3pList))){
    x3pList <- setNames(x3pList,paste0("x3p",1:length(x3pList)))
  }

  if(type == "faceted"){
    surfaceMat_df <- purrr::pmap_dfr(.l = list(x3pList,
                                               names(x3pList)),
                                     function(x3p,name){

                                       x3p$header.info$incrementX <- 1
                                       x3p$header.info$incrementY <- 1

                                       x3p %>%
                                         x3ptools::x3p_to_df() %>%
                                         #perform some transformations on the
                                         #x,y values so that the plot is
                                         #representative of the actual surface
                                         #matrix (i.e., element [1,1] of the
                                         #surface matrix is in the top-left
                                         #corner)
                                         dplyr::mutate(xnew = max(y) - y,
                                                       ynew = max(x) - x,
                                                       value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                                         dplyr::select(-c(x,y)) %>%
                                         dplyr::rename(x=xnew,
                                                       y=ynew) %>%
                                         dplyr::mutate(x3p = rep(name,times = nrow(.)))
                                     }) %>%
      dplyr::mutate(x3p = factor(.data$x3p,levels = names(x3pList)))

    plts <- surfaceMat_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
      ggplot2::scale_fill_gradientn(colours = height.colors,
                                    values = scales::rescale(quantile(surfaceMat_df$value,height.quantiles,na.rm = TRUE)),
                                    breaks = function(lims){
                                      dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                      dat <- dat %>%
                                        setNames(paste0(names(dat)," [",round(dat,1),"]"))

                                      return(dat)
                                    },
                                    na.value = na.value) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
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

    return(plts)
  }
  else if(type == "list"){
    plts <- purrr::pmap(.l = list(x3pList,
                                  names(x3pList)),
                        function(x3p,name){

                          surfaceMat_df <- x3p %>%
                            x3ptools::x3p_to_df() %>%
                            #perform some transformations on the
                            #x,y values so that the plot is
                            #representative of the actual surface
                            #matrix (i.e., element [1,1] of the
                            #surface matrix is in the top-left
                            #corner)
                            dplyr::mutate(xnew = max(y) - y,
                                          ynew = max(x) - x,
                                          value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::select(-c(x,y)) %>%
                            dplyr::rename(x=xnew,
                                          y=ynew) %>%
                            dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          plt <- surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
                            ggplot2::scale_fill_gradientn(colours = height.colors,
                                                          values = scales::rescale(quantile(surfaceMat_df$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                                          breaks = function(lims){
                                                            dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                                            dat <- dat %>%
                                                              setNames(paste0(names(dat)," [",round(dat,1),"]"))

                                                            return(dat)
                                                          },
                                                          na.value = na.value) +
                            ggplot2::theme_minimal() +
                            ggplot2::coord_fixed(expand = FALSE) +
                            ggplot2::theme(
                              axis.title.x = ggplot2::element_blank(),
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
                            ggplot2::labs(title = name)

                          return(plt)
                        })

    return(plts)
  }
}

# helper function for x3pListPlot. Rotates a surface matrix, but doesn't crop
# back to the original surface matrix's dimensions.
rotateSurfaceMatrix_noCrop <- function(surfaceMat,
                                       theta = 0,
                                       interpolation = 0){
  surfaceMatFake <- (surfaceMat*10^5) + 1 #scale and shift all non-NA pixels up 1 (meter)
  # imFakeRotated <- :bilinearInterpolation(imFake,theta)
  surfaceMatFakeRotated <- surfaceMatFake %>%
    imager::as.cimg() %>%
    imager::imrotate(angle = theta,
                     interpolation = interpolation, #linear interpolation,
                     boundary = 0) %>% #pad boundary with 0s (dirichlet condition)
    as.matrix()

  surfaceMatFakeRotated[surfaceMatFakeRotated == 0] <- NA
  #shift all of the legitimate pixels back down by 1:
  surfaceMatRotated <- (surfaceMatFakeRotated - 1)/(10^5)

  return(surfaceMatRotated)
}

# @name linear_to_matrix
# @param index integer vector of indices, must be between 1 and nrow*ncol
# @param nrow number of rows, integer value defaults to 7
# @param ncol  number of columns, integer value, defaults to number of rows
# @param byrow logical value, is linear index folded into matrix by row (default) or by column (`byrow=FALSE`).
# @examples
# index <- sample(nrow*ncol, 10, replace = TRUE)
# linear_to_matrix(index, nrow=4, ncol = 5, byrow=TRUE)
#
# @keywords internal
# @importFrom rlang .data

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

#' Plot a scan partitioned into a grid of cells.
#'
#' @name cellGridPlot
#'
#' @export

cellGridPlot <- function(x3p,
                         numCells = 64,
                         legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
                         height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                         na.value = "gray80"){

  surfaceMat_df <- x3p %>%
    #TODO: there's a more efficient way to do the following that doesn't require
    #splitting the scan up only to recombine it immediately.
    comparison_cellDivision(numCells = numCells) %>%
    purrr::pmap_dfr(~ {

      ..2 %>%
        x3p_to_df() %>%
        mutate(cellIndex = ..1)

    }) %>%
    dplyr::mutate(value = value - median(value,na.rm = TRUE)) %>%
    tidyr::separate(col = cellIndex,into = c("row","col"),sep = ", ") %>%
    dplyr::mutate(col = as.numeric(col),
                  row = as.numeric(row),
                  xnew = max(y) - y,
                  ynew = max(x) - x) %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::rename(x=xnew,
                  y=ynew)

  plt <- surfaceMat_df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
    ggplot2::scale_fill_gradientn(colours = height.colors,
                                  values = scales::rescale(quantile(surfaceMat_df$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                  breaks = function(lims){
                                    dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                    dat <- dat %>%
                                      setNames(paste0(names(dat)," [",round(dat,1),"]"))

                                    return(dat)
                                  },
                                  na.value = na.value) +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5,
                                         size = 11)) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(3,"in"),
                                                    label.theme = ggplot2::element_text(size = 8),
                                                    title.theme = ggplot2::element_text(size = 10),
                                                    frame.colour = "black",
                                                    ticks.colour = "black"),
                    colour =  'none') +
    ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
    facet_grid(rows = vars(row),
               cols = vars(col))

  return(plt)

}

targetCellCorners <- function(alignedTargetCell,cellIndex,theta,cmcClassif,target){

  targetScanRows <- alignedTargetCell$cmcR.info$regionIndices[c(3)] + alignedTargetCell$cmcR.info$regionRows - 1
  targetScanCols <- alignedTargetCell$cmcR.info$regionIndices[c(1)] + alignedTargetCell$cmcR.info$regionCols - 1

  rotatedMask <- cmcR:::rotateSurfaceMatrix(target$surface.matrix,theta)

  rowPad <- 0
  colPad <- 0

  if(targetScanRows[1] <= 0){

    rowPad <- abs(targetScanRows[1]) + 1

    rotatedMask <- rbind(matrix(NA,nrow = rowPad,ncol = ncol(rotatedMask)),
                         rotatedMask)

    targetScanRows <- targetScanRows + rowPad
  }

  if(targetScanCols[1] <= 0){

    colPad <- abs(targetScanCols[1]) + 1

    rotatedMask <- cbind(matrix(NA,nrow = nrow(rotatedMask),ncol = colPad),
                         rotatedMask)

    targetScanCols <- targetScanCols + colPad
  }

  if(targetScanRows[2] > nrow(rotatedMask)){

    rowPad <- targetScanRows[2] - nrow(rotatedMask)

    rotatedMask <- rbind(rotatedMask,
                         matrix(NA,nrow = rowPad,ncol = ncol(rotatedMask)))

  }

  if(targetScanCols[2] > ncol(rotatedMask)){

    colPad <- targetScanCols[2] - ncol(rotatedMask)

    rotatedMask <- cbind(rotatedMask,
                         matrix(NA,nrow = nrow(rotatedMask),ncol = colPad))

  }

  rotatedMask[targetScanRows[1]:targetScanRows[2],targetScanCols[1]:targetScanCols[2]] <- 100

  rotatedMask <- cmcR:::rotateSurfaceMatrix_noCrop(rotatedMask,theta = -1*theta)
  #make a copy that isn't going to have the target cell indices added so that we
  #know exactly how many rows/cols we need to translate everything to get back to
  #the original scan indices
  rotatedMaskCopy <- rotatedMask#cmcR:::rotateSurfaceMatrix_noCrop(rotatedMaskCopy,theta = 0)#-1*(-30))

  rotatedMaskCopy[rotatedMaskCopy == 100] <- NA

  newColPad <- 0
  newRowPad <- 0
  if(theta != 0){

    # the cells that are rotated may have been "shifted" due to the cropping
    # performed above relative to the unrotated target scan's indices -- for
    # example, padding the rotatedMask to the left requires a correction after
    # rotating. Unfortunately, because the cells are rotated, the padding done in
    # the rotated domain doesn't come out to be nice (dRow,dCol) translations in
    # the unrotated domain. We need to perform some trig to determine what
    # (dRow,dCol) in the rotated domain translates to in the unrotated domain.

    # In the rotated domain:
    #    ------------* <<- the location of target cell after padding in rotatedMask
    #    |    ^dx^
    #    |
    #    | <- dy
    #    |
    #    * <<- where the target cell *should* be relative to the rotated target's indices
    #

    # no consider rotating this whole space, then the dx, dy will be "tilted" by
    # some theta and we will need to calculate via trig what the correct dx', dy'
    # are in the original, unrotated domain. Draw a diagram with a rotated dx, dy
    # by some theta and draw a straight line between the two * above and you
    # should be able to work out the following formulas again

    #TODO: pay attention to the necessary signs (up/down/left/right) of the
    #corrections below
    psi <- atan2(rowPad,colPad)
    hyp <- sqrt(colPad^2 + rowPad^2)
    phi <- pi/2 - (theta*pi/180 + psi)

    newColPad <- sin(phi)*hyp
    newRowPad <- cos(phi)*hyp

  }

  ret <- rotatedMask %>%
    imager::as.cimg() %>%
    as.data.frame() %>%
    mutate(xnew = y,
           ynew = x) %>%
    select(-c(x,y)) %>%
    rename(x=xnew,y=ynew) %>%
    mutate(x = x - min(which(colSums(rotatedMaskCopy,na.rm = TRUE) > 0)),
           # x = x - newColPad,
           y = y - min(which(rowSums(rotatedMaskCopy,na.rm = TRUE) > 0)),
           # y = y - newRowPad,
           # y = max(y) - y
           y = nrow(target$surface.matrix) - y
    ) %>%
    filter(value == 100) %>%
    select(-value) %>%
    group_by(x,y) %>%
    distinct() %>%
    mutate(cellIndex = cellIndex,
           theta = theta,
           cmcClassif = cmcClassif)

  return(ret)

}

#' Plot CMCs
#' @name cmcPlot
#' @export
#' @importFrom patchwork wrap_plots
cmcPlot <- function(reference,
                    target,
                    cmcClassifs,
                    cmcCol = "originalMethod"){

  #check that the necessary columns are in cmcClassifs

  stopifnot("Make sure that there is a column called 'cellHeightValues' that is the result of the comparison_alignedTargetCell() function." = any(str_detect(names(cmcClassifs),"cellHeightValues")))

  stopifnot("Make sure that there is a column called 'alignedTargetCell' that is the result of the comparison_alignedTargetCell() function." = any(str_detect(names(cmcClassifs),"alignedTargetCell")))

  stopifnot("Make sure that there is a column called 'cellIndex'" = any(str_detect(names(cmcClassifs),"cellIndex")))

  stopifnot("Make sure that there is a column called 'theta'" = any(str_detect(names(cmcClassifs),"theta")))

  stopifnot("Make sure there is a column called 'pairwiseCompCor'" = any(str_detect(names(cmcClassifs),"pairwiseCompCor")))

  # get the indices for the necessary columns
  referenceCellCol <- which(str_detect(names(cmcClassifs),"cellHeightValues"))

  targetCellCol <- which(str_detect(names(cmcClassifs),"alignedTargetCell"))

  cellIndexCol <- which(str_detect(names(cmcClassifs),"cellIndex"))

  thetaCol <- which(str_detect(names(cmcClassifs),"theta"))

  cmcIndexCol <- which(str_detect(names(cmcClassifs),cmcCol))

  cmcClassifs <- cmcClassifs %>%
    group_by(cellIndex) %>%
    filter(pairwiseCompCor == max(pairwiseCompCor))

  targetCellData <- cmcClassifs %>%
    select(c(targetCellCol,cellIndexCol,thetaCol,cmcIndexCol)) %>%
    pmap_dfr(~ targetCellCorners(alignedTargetCell = ..1,
                                 cellIndex = ..2,
                                 theta = ..3,
                                 cmcClassif = ..4,
                                 target = target))

  referenceCells <- cmcClassifs %>%
    pull(referenceCellCol)

  cellData <- cmcClassifs %>%
    select(c(cellIndexCol,referenceCellCol,cmcIndexCol)) %>%
    pmap_dfr(~ {

      cellInds <- ..2$cmcR.info$cellRange %>%
        str_remove("rows: ") %>%
        str_remove("cols: ") %>%
        str_split(pattern = ", ")

      cellInds_rows <- str_split(cellInds[[1]][1]," - ")[[1]]
      cellInds_cols <- str_split(cellInds[[1]][2]," - ")[[1]]

      return(data.frame(rowStart = as.numeric(cellInds_rows[1]),
                        rowEnd = as.numeric(cellInds_rows[2]),
                        colStart = as.numeric(cellInds_cols[1]),
                        colEnd = as.numeric(cellInds_cols[2])) %>%
               mutate(cellIndex = ..1,
                      originalMethod = ..3))

    }) %>%
    mutate(rowStart = max(rowEnd) - rowStart,
           rowEnd = max(rowEnd) - rowEnd,
           colMean = map2_dbl(colStart,colEnd,~ mean(c(.x,.y))),
           rowMean = map2_dbl(rowStart,rowEnd,~ mean(c(.x,.y)))) %>%
    rename(cmcClassif = cmcCol)

  refPlt <- x3pListPlot(list("reference" = reference),
                        height.colors = colorspace::desaturate(rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                                                                     '#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')))) +
    ggnewscale::new_scale_fill() +
    geom_rect(data = cellData,
              aes(xmin = colStart,xmax = colEnd,ymin = rowStart,ymax = rowEnd,fill = cmcClassif),
              alpha = .2,
              inherit.aes = FALSE) +
    scale_fill_manual(values = c("#313695","#a50026")) +
    geom_text(data = cellData,
              aes(x = colMean,y = rowMean,label = cellIndex),inherit.aes = FALSE) +
    guides(fill = ggplot2::guide_legend(order = 1)) +
    theme(
      legend.direction = "horizontal"
      ) +
    labs(fill = "CMC Classif.")

  cmcLegend <- ggplotify::as.ggplot(cowplot::get_legend(refPlt)$grobs[[1]])

  refPlt <- refPlt +
    theme(legend.position = "none")

  # refPlt <- cmcR::x3pListPlot(list(reference) %>% set_names("reference"),
  #                             height.colors = colorspace::desaturate(rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')))) +
  #   theme(legend.position = "none")

  plt <- cmcR::x3pListPlot(list("target" = target),
                           height.colors = colorspace::desaturate(rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')))) +
    theme(legend.position = "none")

  plt <- plt +
    ggnewscale::new_scale_fill() +
    geom_raster(data = targetCellData,
                aes(x = x,y = y,fill = cmcClassif),
                alpha = .2) +
    scale_fill_manual(values = c("#313695","#a50026")) +
    geom_text(data = targetCellData %>%
                group_by(cellIndex) %>%
                summarise(x = mean(x),
                          y = mean(y),
                          theta = unique(theta)),
              aes(x=x,y=y,label = cellIndex,angle = -1*theta))

  # library(patchwork)
  # return((refPlt | plt))
  return(patchwork::wrap_plots(refPlt,plt,cmcLegend,nrow = 2,heights = c(1,.1)))
}

# @name arrangeCMCPlot
#
# @keywords internal
#
# @importFrom stats median setNames
# @importFrom rlang .data
#' @importFrom stringr str_remove_all

# arrangeCMCPlot <- function(reference,
#                            target,
#                            allCells,
#                            x3pNames,
#                            pltType = "faceted",
#                            legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
#                            height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
#                            cell.colors = c("#a50026","#313695"),
#                            cell.alpha = .2,
#                            na.value = "gray80"){
#
#   target_cellGrid <- allCells %>%
#     dplyr::mutate(firstRow = (reference$header.info$incrementY)*(.data$firstRow),
#                   lastRow = (reference$header.info$incrementY)*(.data$lastRow),
#                   firstCol = (reference$header.info$incrementY)*(.data$firstCol),
#                   lastCol = (reference$header.info$incrementY)*(.data$lastCol)) %>%
#     dplyr::mutate(x_1 = .data$firstCol,
#                   y_1 = .data$firstRow,
#                   x_2 = .data$lastCol,
#                   y_2 = .data$firstRow,
#                   x_3 = .data$lastCol,
#                   y_3 = .data$lastRow,
#                   x_4 = .data$firstCol,
#                   y_4 = .data$lastRow) %>%
#     tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
#                         names_to = c(".value","order"),
#                         names_pattern = "(.+)_(.+)") %>%
#     dplyr::mutate(midCol = (.data$lastCol + .data$firstCol)/2,
#                   midRow = (.data$lastRow + .data$firstRow)/2,
#                   x3p = rep(x3pNames[1],times = nrow(.)),
#                   theta = rep(0,times = nrow(.)),
#                   cellIndex = stringr::str_remove_all(string = cellIndex,pattern = " "))
#
#   reference_cellGrid <- allCells %>%
#     dplyr::mutate(firstRow = (target$header.info$incrementY)*(.data$firstRow),
#                   lastRow = (target$header.info$incrementY)*(.data$lastRow),
#                   firstCol = (target$header.info$incrementY)*(.data$firstCol),
#                   lastCol = (target$header.info$incrementY)*(.data$lastCol)) %>%
#     dplyr::mutate(firstRowCentered = .data$firstRow - max(.data$lastRow)/2,
#                   lastRowCentered = .data$lastRow - max(.data$lastRow)/2,
#                   firstColCentered = .data$firstCol - max(.data$lastCol)/2,
#                   lastColCentered = .data$lastCol - max(.data$lastCol)/2) %>%
#     dplyr::mutate(topLeftCorner_col = .data$firstColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$lastRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (target$header.info$incrementY)*.data$x/2,
#                   topLeftCorner_row = .data$firstColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$lastRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (target$header.info$incrementY)*.data$y/2,
#                   topRightCorner_col = .data$lastColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$lastRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (target$header.info$incrementY)*.data$x/2,
#                   topRightCorner_row = .data$lastColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$lastRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (target$header.info$incrementY)*.data$y/2,
#                   bottomRightCorner_col = .data$lastColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$firstRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (target$header.info$incrementY)*.data$x/2,
#                   bottomRightCorner_row = .data$lastColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$firstRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (target$header.info$incrementY)*.data$y/2,
#                   bottomLeftCorner_col = .data$firstColCentered*cos((.data$theta - median(.data$theta))*(pi/180)) - .data$firstRowCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastCol)/2 - (target$header.info$incrementY)*.data$x/2,
#                   bottomLeftCorner_row = .data$firstColCentered*sin((.data$theta - median(.data$theta))*(pi/180)) + .data$firstRowCentered*cos((.data$theta - median(.data$theta))*(pi/180)) + max(.data$lastRow)/2 - (target$header.info$incrementY)*.data$y/2) %>%
#     #this is redundant, but are the names attributed to the x and y columns are
#     #set-up down below, so I won't change it
#     dplyr::mutate(x_1 = .data$topLeftCorner_col,
#                   y_1 = .data$topLeftCorner_row,
#                   x_2 = .data$topRightCorner_col,
#                   y_2 = .data$topRightCorner_row,
#                   x_3 = .data$bottomRightCorner_col,
#                   y_3 = .data$bottomRightCorner_row,
#                   x_4 = .data$bottomLeftCorner_col,
#                   y_4 = .data$bottomLeftCorner_row) %>%
#     tidyr::pivot_longer(cols = tidyr::starts_with(c("x","y")),
#                         names_to = c(".value","order"),
#                         names_pattern = "(.+)_(.+)") %>%
#     dplyr::mutate(midCol = (.data$topRightCorner_col + .data$bottomLeftCorner_col)/2,
#                   midRow = (.data$topRightCorner_row + .data$bottomLeftCorner_row)/2,
#                   x3p = rep(x3pNames[2],times = nrow(.)),
#                   theta = .data$theta - median(.data$theta),
#                   cellIndex = stringr::str_remove_all(string = cellIndex,pattern = " "))
#
#   x3pPlt <- x3pListPlot(x3pList = list(reference,target) %>%
#                           setNames(x3pNames),
#                         type = pltType,
#                         rotate = c(90,90 + median(allCells$theta)),
#                         legend.quantiles = legend.quantiles,
#                         height.colors = height.colors,
#                         na.value = na.value,
#                         guide = "none")
#
#   if(pltType == "faceted"){
#
#     x3pPlt <- x3pPlt +
#       ggnewscale::new_scale_fill() +
#       ggplot2::geom_polygon(data = target_cellGrid,
#                             mapping = ggplot2::aes(x = .data$x,
#                                                    y = .data$y,
#                                                    group = .data$cellIndex,
#                                                    fill = .data$cmc),
#                             alpha = cell.alpha,
#                             size = 2) +
#       ggplot2::geom_polygon(data = reference_cellGrid,
#                             mapping = ggplot2::aes(x = .data$x,
#                                                    y = .data$y,
#                                                    group = .data$cellIndex,
#                                                    fill = .data$cmc),
#                             alpha = cell.alpha,
#                             size = 2) +
#       ggplot2::geom_text(data = dplyr::bind_rows(target_cellGrid,
#                                                  reference_cellGrid),
#                          ggplot2::aes(x = .data$midCol,
#                                       y = .data$midRow,
#                                       label = .data$cellIndex,
#                                       colour = .data$cmc,
#                                       angle = .data$theta),
#                          size = 3) +
#       ggplot2::scale_colour_manual(values = cell.colors,
#                                    aesthetics = c("fill","colour")) +
#       ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
#       ggplot2::theme(legend.position = "bottom")
#   }
#   else if(pltType == "list"){
#     x3pPlt[[1]] <- x3pPlt[[1]] +
#       ggnewscale::new_scale_fill() +
#       ggplot2::geom_polygon(data = target_cellGrid,
#                             mapping = ggplot2::aes(x = .data$x,
#                                                    y = .data$y,
#                                                    group = .data$cellIndex,
#                                                    fill = .data$cmc),
#                             alpha = cell.alpha,
#                             size = 2) +
#       ggplot2::scale_colour_manual(values = cell.colors,
#                                    aesthetics = c("fill","colour")) +
#       ggplot2::geom_text(data = target_cellGrid,
#                          ggplot2::aes(x = .data$midCol,
#                                       y = .data$midRow,
#                                       label = .data$cellIndex,
#                                       colour = .data$cmc,
#                                       angle = .data$theta),
#                          size = 3) +
#       ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
#       ggplot2::theme(legend.position = "bottom")
#
#     x3pPlt[[2]] <- x3pPlt[[2]] +
#       ggnewscale::new_scale_fill() +
#       ggplot2::geom_polygon(data = reference_cellGrid,
#                             mapping = ggplot2::aes(x = .data$x,
#                                                    y = .data$y,
#                                                    group = .data$cellIndex,
#                                                    fill = .data$cmc),
#                             alpha = cell.alpha,
#                             size = 2) +
#       ggplot2::scale_colour_manual(values = cell.colors,
#                                    aesthetics = c("fill","colour")) +
#       ggplot2::geom_text(data = reference_cellGrid,
#                          ggplot2::aes(x = .data$midCol,
#                                       y = .data$midRow,
#                                       label = .data$cellIndex,
#                                       colour = .data$cmc,
#                                       angle = .data$theta),
#                          size = 3) +
#       ggplot2::guides(fill = ggplot2::guide_legend(title = "Cell Type")) +
#       ggplot2::theme(legend.position = "bottom")
#   }
#
#   return(x3pPlt)
# }

# #'Visualize initial and high CMCs for a cartridge case pair comparison
# #'@name cmcPlot
# #'
# #'@description Constructs either a single faceted plot or a list of plots
# #'  depicting the CMCs/non-CMCs under the initially proposed and High CMC
# #'  methods for a pair of cartridge case scans
# #'@param reference an x3p object
# #'@param target a different x3p object
# #'@param reference_v_target_CMCs CMCs for the comparison between the reference
# #'  scan and the target scan.
# #'@param target_v_reference_CMCs (optional) CMCs for the comparison between the
# #'  target scan and the reference scan. If this is missing, then only the
# #'  original method CMCs will be plotted
# #'@param corColName name of correlation similarity score column used to identify
# #'  the CMCs in the two comparison_*_df data frames (e.g., pairwiseCompCor)
# #'@param type argument to be passed to cmcR::x3pListPlot function
# #'@param x3pNames (Optional) Names of x3p objects to be included in x3pListPlot
# #'  function
# #'@param legend.quantiles vector of quantiles to be shown as tick marks on
# #'  legend plot
# #'@param height.colors vector of colors to be passed to scale_fill_gradientn
# #'  that dictates the height value colorscale
# #'@param cell.colors vector of 2 colors for plotting non-matching and matching
# #'  (in that order) cells
# #'@param cell.alpha sets alpha of cells (passed to geom_polygon)
# #'@param numCells the size of the grid used to compare the reference and target
# #'  scans. Must be a perfect square.
# #'@param na.value color to be used for NA values (passed to
# #'  scale_fill_gradientn)
# #'@return A list of 4 ggplot objects showing the CMCs identified under both
# #'  decision rules and in both comparison directions.
# #'@examples
# #'#Takes > 5 seconds to run
# #'\donttest{
# #'data(fadul1.1_processed,fadul1.2_processed)
# #'
# #'comparisonDF_1to2 <- purrr::map_dfr(seq(-30,30,by = 3),
# #'                                    ~ comparison_allTogether(fadul1.1_processed,
# #'                                                        fadul1.2_processed,
# #'                                                        theta = .))
# #'comparisonDF_2to1 <- purrr::map_dfr(seq(-30,30,by = 3),
# #'                                    ~ comparison_allTogether(fadul1.2_processed,
# #'                                                        fadul1.1_processed,
# #'                                                        theta = .))
# #'
# #'comparisonDF_1to2 <- comparisonDF_1to2 %>%
# #' dplyr::mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
# #'                                                    x = x,
# #'                                                    y = y,
# #'                                                    theta = theta,
# #'                                                    corr = pairwiseCompCor),
# #'               highCMCClassif = decision_CMC(cellIndex = cellIndex,
# #'                                            x = x,
# #'                                            y = y,
# #'                                            theta = theta,
# #'                                            corr = pairwiseCompCor,
# #'                                            tau = 1))
# #'
# #'
# #'comparisonDF_2to1 <- comparisonDF_2to1 %>%
# #' dplyr::mutate(originalMethodClassif = decision_CMC(cellIndex = cellIndex,
# #'                                                    x = x,
# #'                                                    y = y,
# #'                                                    theta = theta,
# #'                                                    corr = pairwiseCompCor),
# #'               highCMCClassif = decision_CMC(cellIndex = cellIndex,
# #'                                            x = x,
# #'                                            y = y,
# #'                                            theta = theta,
# #'                                            corr = pairwiseCompCor,
# #'                                            tau = 1))
# #'
# #'cmcPlot(fadul1.1_processed,
# #'        fadul1.2_processed,
# #'        comparisonDF_1to2,
# #'        comparisonDF_2to1,
# #'        corColName = "pairwiseCompCor")
# #'}
# #'@importFrom utils hasName
# #'@importFrom rlang .data
# #'@export

# cmcPlot <- function(reference,
#                     target,
#                     reference_v_target_CMCs,
#                     target_v_reference_CMCs = reference_v_target_CMCs,
#                     corColName = "pairwiseCompCor",
#                     type = "faceted",
#                     x3pNames = c("reference","target"),
#                     legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
#                     height.colors = c("#1B1B1B","#404040","#7B7B7B","#B0B0B0","#DBDBDB","#F7F7F7","#E4E4E4","#C5C5C5","#999999","#717171","#4E4E4E"),
#                     cell.colors = c("#a50026","#313695"),
#                     cell.alpha = .2,
#                     numCells = 64,
#                     na.value = "gray80"){
#
#   reference_cellCorners <- reference %>%
#     comparison_cellDivision(numCells = numCells) %>%
#     purrr::pmap_dfr(~ {
#       idNum <- ..2$cmcR.info$cellRange %>%
#         stringr::str_extract_all(string = ..2$cmcR.info$cellRange,
#                                  pattern = "[0-9]{1,}") %>%
#         unlist() %>%
#         as.numeric()
#
#       data.frame(cellIndex = ..1,
#                  firstRow = idNum[1],
#                  lastRow = idNum[2],
#                  firstCol = idNum[3],
#                  lastCol = idNum[4],
#                  stringsAsFactors = FALSE)
#     })
#
#   target_cellCorners <- target %>%
#     comparison_cellDivision(numCells = numCells) %>%
#     purrr::pmap_dfr(~ {
#       idNum <- ..2$cmcR.info$cellRange %>%
#         stringr::str_extract_all(string = ..2$cmcR.info$cellRange,
#                                  pattern = "[0-9]{1,}") %>%
#         unlist() %>%
#         as.numeric()
#
#       data.frame(cellIndex = ..1,
#                  firstRow = idNum[1],
#                  lastRow = idNum[2],
#                  firstCol = idNum[3],
#                  lastCol = idNum[4],
#                  stringsAsFactors = FALSE)
#     })
#
#   if(nrow(reference_v_target_CMCs) == 0){
#
#     originalMethodCMCs <- as.data.frame(matrix(nrow = 0,ncol = ncol(reference_v_target_CMCs))) %>%
#       setNames(names(reference_v_target_CMCs))
#
#   }
#   else{
#     originalMethodCMCs <- reference_v_target_CMCs %>%
#       dplyr::filter(.data$originalMethodClassif == "CMC")
#   }
#
#   nonoriginalMethodCMCs <- reference_v_target_CMCs %>%
#     dplyr::filter(!(.data$cellIndex %in% originalMethodCMCs$cellIndex))
#
#   if(nrow(nonoriginalMethodCMCs) > 0){
#     nonoriginalMethodCMCs <- nonoriginalMethodCMCs %>%
#       dplyr::group_by(.data$cellIndex) %>%
#       dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
#   }
#
#   allInitialCells_reference_v_target <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
#     dplyr::mutate(cmc = ifelse(.data$originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
#     dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","Original Method CMC"))) %>%
#     dplyr::left_join(reference_cellCorners,
#                      by = "cellIndex")
#
#   originalMethodCMCsPlt_reference_v_target <- arrangeCMCPlot(reference = reference,
#                                                              target = target,
#                                                              allCells = allInitialCells_reference_v_target,
#                                                              x3pNames = x3pNames,
#                                                              pltType = type,
#                                                              legend.quantiles = legend.quantiles,
#                                                              height.colors = height.colors,
#                                                              cell.colors = cell.colors,
#                                                              cell.alpha = cell.alpha,
#                                                              na.value = na.value)
#
#   #If only data for one comparison direction were given, only plot the original
#   #method CMCs in that direction
#   if(assertthat::are_equal(reference_v_target_CMCs, target_v_reference_CMCs)){
#
#     return(list("CMCs" = originalMethodCMCsPlt_reference_v_target))
#
#   }
#
#   #otherwise, create the original CMC plots for the other direction and return:
#
#   if(nrow(target_v_reference_CMCs) == 0){
#
#     originalMethodCMCs <- as.data.frame(matrix(nrow = 0,ncol = ncol(target_v_reference_CMCs))) %>%
#       setNames(names(target_v_reference_CMCs))
#
#   }
#   else{
#     originalMethodCMCs <- target_v_reference_CMCs %>%
#       dplyr::filter(.data$originalMethodClassif == "CMC")
#   }
#
#   nonoriginalMethodCMCs <- target_v_reference_CMCs %>%
#     dplyr::filter(!(.data$cellIndex %in% originalMethodCMCs$cellIndex)) %>%
#     dplyr::group_by(.data$cellIndex) %>%
#     dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
#
#   allInitialCells_target_v_reference <- dplyr::bind_rows(originalMethodCMCs,nonoriginalMethodCMCs) %>%
#     dplyr::mutate(cmc = ifelse(.data$originalMethodClassif == "CMC","Original Method CMC","non-CMC")) %>%
#     dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","Original Method CMC"))) %>%
#     dplyr::left_join(reference_cellCorners,
#                      by = "cellIndex")
#
#   originalMethodCMCsPlt_target_v_reference <- arrangeCMCPlot(reference = target,
#                                                              target = reference,
#                                                              allCells = allInitialCells_target_v_reference,
#                                                              x3pNames = rev(x3pNames),
#                                                              pltType = type,
#                                                              legend.quantiles = legend.quantiles,
#                                                              height.colors = height.colors,
#                                                              cell.colors = cell.colors,
#                                                              cell.alpha = cell.alpha,
#                                                              na.value = na.value)
#
#   if(!hasName(reference_v_target_CMCs,"highCMCClassif") | !hasName(target_v_reference_CMCs,"highCMCClassif")){
#
#     return(list("originalMethodCMCs_reference_v_target" = originalMethodCMCsPlt_reference_v_target,
#                 "originalMethodCMCs_target_v_reference" = originalMethodCMCsPlt_target_v_reference))
#
#   }
#
#   #If the necessary data to construct the High CMCs were given, then plot them
#   #too.
#
#   highCMCs_reference_v_target <- reference_v_target_CMCs %>%
#     dplyr::filter(.data$highCMCClassif == "CMC")
#
#   #Remaining cells not identified as High CMCs
#   non_highCMCs_reference_v_target <- reference_v_target_CMCs %>%
#     dplyr::filter(!(.data$cellIndex %in% highCMCs_reference_v_target$cellIndex)) %>%
#     dplyr::group_by(.data$cellIndex) %>%
#     dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
#
#   highCMC_plotData_reference_v_target <- dplyr::bind_rows(highCMCs_reference_v_target,
#                                                           non_highCMCs_reference_v_target) %>%
#     dplyr::mutate(cmc = ifelse(.data$highCMCClassif == "CMC","High CMC","non-CMC")) %>%
#     dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","High CMC"))) %>%
#     dplyr::left_join(reference_cellCorners,
#                      by = "cellIndex")
#
#   highCMCPlt_reference_v_target <- arrangeCMCPlot(reference = reference,
#                                                   target = target,
#                                                   allCells = highCMC_plotData_reference_v_target,
#                                                   x3pNames = x3pNames,
#                                                   pltType = type,
#                                                   legend.quantiles = legend.quantiles,
#                                                   height.colors = height.colors,
#                                                   cell.colors = cell.colors,
#                                                   cell.alpha = cell.alpha,
#                                                   na.value = na.value)
#
#   #Different High CMCs may have been identified in the other direction -- we
#   #need to plot those separately
#
#   highCMCs_target_v_reference <- target_v_reference_CMCs %>%
#     dplyr::filter(.data$highCMCClassif == "CMC")
#
#   #Remaining cells not identified as High CMCs
#   non_highCMCs_target_v_reference <- target_v_reference_CMCs %>%
#     dplyr::filter(!(.data$cellIndex %in% highCMCs_target_v_reference$cellIndex)) %>%
#     dplyr::group_by(.data$cellIndex) %>%
#     dplyr::filter((!!as.name(corColName)) == max((!!as.name(corColName))))
#
#   highCMC_plotData_target_v_reference <- dplyr::bind_rows(highCMCs_target_v_reference,
#                                                           non_highCMCs_target_v_reference) %>%
#     dplyr::mutate(cmc = ifelse(.data$highCMCClassif == "CMC","High CMC","non-CMC")) %>%
#     dplyr::mutate(cmc = factor(.data$cmc,levels = c("non-CMC","High CMC"))) %>%
#     dplyr::left_join(target_cellCorners,
#                      by = "cellIndex")
#
#   highCMCPlt_target_v_reference <- arrangeCMCPlot(reference = target,
#                                                   target = reference,
#                                                   allCells = highCMC_plotData_target_v_reference,
#                                                   x3pNames = rev(x3pNames),
#                                                   pltType = type,
#                                                   legend.quantiles = legend.quantiles,
#                                                   height.colors = height.colors,
#                                                   cell.colors = cell.colors,
#                                                   cell.alpha = cell.alpha,
#                                                   na.value = na.value)
#
#
#   return(list("originalMethodCMCs_reference_v_target" = originalMethodCMCsPlt_reference_v_target,
#               "originalMethodCMCs_target_v_reference" = originalMethodCMCsPlt_target_v_reference,
#               "highCMC_reference_v_target"= highCMCPlt_reference_v_target,
#               "highCMC_target_v_reference"= highCMCPlt_target_v_reference))
#
# }
