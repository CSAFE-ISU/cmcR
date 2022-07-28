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
#'@param legend.quantiles vector of quantiles to be shown as tick marks on
#'  legend plot
#'@param height.quantiles vector of quantiles associated with each color defined
#'  in the height.colors argument
#'@param height.colors vector of colors to be passed to scale_fill_gradientn
#'  that dictates the height value colorscale
#'@param na.value color to be used for NA values (passed to
#'  scale_fill_gradientn)
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
                        na.value = "gray65"){

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
                                         dplyr::mutate(xnew = max(.data$y) - .data$y,
                                                       ynew = max(.data$x) - .data$x,
                                                       value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                                         dplyr::select(-c(.data$x,.data$y)) %>%
                                         dplyr::rename(x=.data$xnew,
                                                       y=.data$ynew) %>%
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
                                        setNames(paste0(names(dat)," [",round(dat,3),"]"))

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
                            dplyr::mutate(xnew = max(.data$y) - .data$y,
                                          ynew = max(.data$x) - .data$x,
                                          value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::select(-c(.data$x,.data$y)) %>%
                            dplyr::rename(x=.data$xnew,
                                          y=.data$ynew) %>%
                            dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
                            dplyr::mutate(x3p = rep(name,times = nrow(.)))

                          surfaceMat_df %>%
                            ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
                            ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
                            ggplot2::scale_fill_gradientn(colours = height.colors,
                                                          values = scales::rescale(quantile(surfaceMat_df$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
                                                          breaks = function(lims){
                                                            dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)

                                                            dat <- dat %>%
                                                              setNames(paste0(names(dat)," [",round(dat,3),"]"))

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
                        })
  }
  return(plts)
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

# #' Plot a scan partitioned into a grid of cells.
# #'
# #' @name cellGridPlot
# #'
# #' @export

# cellGridPlot <- function(x3p,
#                          numCells = c(8,8),
#                          legend.quantiles = c(0,.01,.25,.5,.75,.99,1),
#                          height.colors = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
#                          na.value = "gray65"){
#
#   surfaceMat_df <- x3p %>%
#     #TODO: there's a more efficient way to do the following that doesn't require
#     #splitting the scan up only to recombine it immediately.
#     comparison_cellDivision(numCells = numCells) %>%
#     purrr::pmap_dfr(~ {
#
#       ..2 %>%
#         x3ptools::x3p_to_df() %>%
#         dplyr::mutate(cellIndex = ..1)
#
#     }) %>%
#     dplyr::mutate(value = .data$value - median(.data$value,na.rm = TRUE)) %>%
#     tidyr::separate(col = cellIndex,into = c("row","col"),sep = ", ") %>%
#     dplyr::mutate(col = as.numeric(col),
#                   row = as.numeric(row),
#                   xnew = max(.data$y) - .data$y,
#                   ynew = max(.data$x) - .data$x) %>%
#     dplyr::select(-c(.data$x,.data$y)) %>%
#     dplyr::rename(x=.data$xnew,
#                   y=.data$ynew)
#
#   plt <- surfaceMat_df %>%
#     ggplot2::ggplot(ggplot2::aes(x = .data$x,y = .data$y)) +
#     ggplot2::geom_raster(ggplot2::aes(fill = .data$value))  +
#     ggplot2::scale_fill_gradientn(colours = height.colors,
#                                   values = scales::rescale(quantile(surfaceMat_df$value,c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),na.rm = TRUE)),
#                                   breaks = function(lims){
#                                     dat <- quantile(surfaceMat_df$value,legend.quantiles,na.rm = TRUE)
#
#                                     dat <- dat %>%
#                                       setNames(paste0(names(dat)," [",round(dat,1),"]"))
#
#                                     return(dat)
#                                   },
#                                   na.value = na.value) +
#     ggplot2::theme_minimal() +
#     ggplot2::coord_fixed(expand = FALSE) +
#     ggplot2::theme(
#       axis.title = ggplot2::element_blank(),
#       axis.text = ggplot2::element_blank(),
#       axis.ticks = ggplot2::element_blank(),
#       panel.grid.major = ggplot2::element_blank(),
#       panel.grid.minor = ggplot2::element_blank(),
#       panel.background = ggplot2::element_blank(),
#       strip.background = ggplot2::element_blank(),
#       strip.text = ggplot2::element_blank(),
#       plot.title = ggplot2::element_text(hjust = .5,
#                                          size = 11)) +
#     ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(3,"in"),
#                                                     label.theme = ggplot2::element_text(size = 8),
#                                                     title.theme = ggplot2::element_text(size = 10),
#                                                     frame.colour = "black",
#                                                     ticks.colour = "black"),
#                     colour =  'none') +
#     ggplot2::labs(fill = expression("Rel. Height ["*mu*"m]")) +
#     ggplot2::facet_grid(rows = ggplot2::vars(row),
#                         cols = ggplot2::vars(col))
#
#   return(plt)
#
# }

targetCellCorners <- function(alignedTargetCell,cellIndex,theta,cmcClassif,target){

  targetScanRows <- alignedTargetCell$cmcR.info$regionIndices[c(3)] + alignedTargetCell$cmcR.info$regionRows - 1
  targetScanCols <- alignedTargetCell$cmcR.info$regionIndices[c(1)] + alignedTargetCell$cmcR.info$regionCols - 1

  rotatedMask <- rotateSurfaceMatrix(target$surface.matrix,theta)

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

  rotatedMask <- rotateSurfaceMatrix_noCrop(rotatedMask,theta = -1*theta)
  #make a copy that isn't going to have the target cell indices added so that we
  #know exactly how many rows/cols we need to translate everything to get back to
  #the original scan indices
  rotatedMaskCopy <- rotatedMask#rotateSurfaceMatrix_noCrop(rotatedMaskCopy,theta = 0)#-1*(-30))

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
    dplyr::mutate(xnew = .data$y,
                  ynew = .data$x) %>%
    dplyr::select(-c(.data$x,.data$y)) %>%
    dplyr::rename(x=.data$xnew,y=.data$ynew) %>%
    dplyr::mutate(x = .data$x - min(which(colSums(abs(rotatedMaskCopy),na.rm = TRUE) > 0)),
                  y = .data$y - min(which(rowSums(abs(rotatedMaskCopy),na.rm = TRUE) > 0)),
                  y = nrow(target$surface.matrix) - .data$y
    ) %>%
    dplyr::filter(.data$value == 100) %>%
    dplyr::select(-.data$value) %>%
    dplyr::group_by(.data$x,.data$y) %>%
    dplyr::distinct() %>%
    dplyr::mutate(cellIndex = cellIndex,
                  theta = theta,
                  cmcClassif = cmcClassif)

  return(ret)

}

#'Plot Congruent Matching Cells results for a pair of cartridge cases.
#'@name cmcPlot
#'
#'@param reference the scan that is partitioned into a grid of cells
#'@param target the scan to which each reference cell is compared during the
#'  cell-based comparison procedure
#'@param cmcClassifs a data frame containing columns cellHeightValues,
#'  alignedTargetCell, cellIndex, theta, and user-defined cmcCol & corrCol
#'@param type the form of the returned plot object(s). Either "faceted," meaning
#'  the reference and target plot will be shown side-by-side or "list" meaning
#'  each element of the plot (referece, target, and legend) will be returned
#'  separately as elements of a list
#'@param cmcCol name of column containing CMC classifications as returned by the
#'  decision_CMC function. Defaults to "originalMethod"
#'@param corrCol name of column containing correlation values for each cell.
#'  Defaults to "pairwiseCompCor," but "fft_ccf" is a common alternative.
#'@export
#'@importFrom patchwork wrap_plots
#'@importFrom ggplotify as.ggplot
cmcPlot <- function(reference,
                    target,
                    cmcClassifs,
                    type = "faceted",
                    cmcCol = "originalMethod",
                    corrCol = "pairwiseCompCor"){

  #check that the necessary columns are in cmcClassifs

  stopifnot("Make sure that there is a column called 'cellHeightValues' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"cellHeightValues")))

  stopifnot("Make sure that there is a column called 'alignedTargetCell' that is the result of the comparison_alignedTargetCell() function." = any(stringr::str_detect(names(cmcClassifs),"alignedTargetCell")))

  stopifnot("Make sure that there is a column called 'cellIndex'" = any(stringr::str_detect(names(cmcClassifs),"cellIndex")))

  stopifnot("Make sure that there is a column called 'theta'" = any(stringr::str_detect(names(cmcClassifs),"theta")))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),cmcCol)))

  stopifnot(any(stringr::str_detect(names(cmcClassifs),corrCol)))

  # get the indices for the necessary columns
  referenceCellCol <- which(stringr::str_detect(names(cmcClassifs),"cellHeightValues"))

  targetCellCol <- which(stringr::str_detect(names(cmcClassifs),"alignedTargetCell"))

  cellIndexCol <- which(stringr::str_detect(names(cmcClassifs),"cellIndex"))

  thetaCol <- which(stringr::str_detect(names(cmcClassifs),"theta"))

  cmcIndexCol <- which(stringr::str_detect(names(cmcClassifs),cmcCol))

  # cmcClassifs <- cmcClassifs %>%
  #   dplyr::group_by(cellIndex) %>%
  #   dplyr::filter(!!as.name(corrCol) == max(!!as.name(corrCol)))

  targetCellData <- cmcClassifs %>%
    dplyr::select(c(targetCellCol,cellIndexCol,thetaCol,cmcIndexCol)) %>%
    purrr::pmap_dfr(~ targetCellCorners(alignedTargetCell = ..1,
                                        cellIndex = ..2,
                                        theta = ..3,
                                        cmcClassif = ..4,
                                        target = target))

  referenceCells <- cmcClassifs %>%
    dplyr::pull(referenceCellCol)

  cellData <- cmcClassifs %>%
    dplyr::select(c(cellIndexCol,referenceCellCol,cmcIndexCol)) %>%
    purrr::pmap_dfr(~ {

      cellInds <- ..2$cmcR.info$cellRange %>%
        stringr::str_remove("rows: ") %>%
        stringr::str_remove("cols: ") %>%
        stringr::str_split(pattern = ", ")

      cellInds_rows <- stringr::str_split(cellInds[[1]][1]," - ")[[1]]
      cellInds_cols <- stringr::str_split(cellInds[[1]][2]," - ")[[1]]

      return(data.frame(rowStart = as.numeric(cellInds_rows[1]),
                        rowEnd = as.numeric(cellInds_rows[2]),
                        colStart = as.numeric(cellInds_cols[1]),
                        colEnd = as.numeric(cellInds_cols[2])) %>%
               dplyr::mutate(cellIndex = ..1,
                             cmcClassif = ..3))

    }) %>%
    dplyr::mutate(rowStart = max(.data$rowEnd) - .data$rowStart,
                  rowEnd = max(.data$rowEnd) - .data$rowEnd,
                  colMean = purrr::map2_dbl(.data$colStart,.data$colEnd,~ mean(c(.x,.y))),
                  rowMean = purrr::map2_dbl(.data$rowStart,.data$rowEnd,~ mean(c(.x,.y))))

  # ggplot2 complains about the guides
  suppressWarnings({

    plt <- x3pListPlot(list("target" = target),
                       height.colors =
                         c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
      ggplot2::theme(legend.position = "none")

    if(all(targetCellData$cmcClassif == "CMC")){

      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#313695")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#313695")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }
    else if(all(targetCellData$cmcClassif == "non-CMC")){

      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#a50026")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#a50026")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }
    else{
      refPlt <- x3pListPlot(list("reference" = reference),
                            height.colors =
                              c('#1B1B1B','#404040','#7B7B7B','#B0B0B0','#DBDBDB','#F7F7F7','#E4E4E4','#C5C5C5','#999999','#717171','#4E4E4E')) +
        ggplot2::guides(fill = "none") +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = cellData,
                           ggplot2::aes(xmin = .data$colStart,xmax = .data$colEnd,ymin = .data$rowStart,ymax = .data$rowEnd,fill = .data$cmcClassif),
                           alpha = .2,
                           inherit.aes = FALSE) +
        ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
        ggplot2::geom_text(data = cellData,
                           ggplot2::aes(x = .data$colMean,y = .data$rowMean,label = .data$cellIndex),inherit.aes = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
        ggplot2::theme(
          legend.direction = "horizontal"
        ) +
        ggplot2::labs(fill = "CMC Classif.")

      plt <- plt +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(data = targetCellData,
                             ggplot2::aes(x = .data$x,y = .data$y,fill = .data$cmcClassif),
                             alpha = .2) +
        ggplot2::scale_fill_manual(values = c("#313695","#a50026")) +
        ggplot2::geom_text(data = targetCellData %>%
                             dplyr::group_by(.data$cellIndex) %>%
                             dplyr::summarise(x = mean(.data$x),
                                              y = mean(.data$y),
                                              theta = unique(.data$theta)),
                           ggplot2::aes(x=.data$x,y=.data$y,label = .data$cellIndex,angle = -1*.data$theta))

    }

    cmcLegend <- ggplotify::as.ggplot(cowplot::get_legend(refPlt)$grobs[[1]])

    refPlt <- refPlt +
      ggplot2::theme(legend.position = "none")

  })

  # library(patchwork)
  # return((refPlt | plt))
  if(type == "list"){
    return(list("reference" = refPlt,
                "target" = plt,
                "legend" = cmcLegend))
  }

  return(patchwork::wrap_plots(refPlt,plt,cmcLegend,nrow = 2,heights = c(1,.1)))
}
