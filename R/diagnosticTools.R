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
#' @export

#TODO: add theta and cellID to the plot table

ccfMapPlot <- function(mat1,
                       mat2,
                       returnGrob = FALSE){

  ccfMat <- cmcR:::ccfMap(mat1,mat2)

  ccfDF <- ccfMat %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    imager::as.cimg() %>%
    as.data.frame() %>%
    dplyr::mutate(dx = x - max(x)/2,
                  dy = y - max(y)/2) %>%
    dplyr::rename(CCF = value)

  ccfMaxInfo <- ccfDF %>%
    dplyr::filter(CCF == max(CCF)) %>%
    dplyr::mutate(CCF = round(CCF,3))

  mat1 <- (mat1 - mean(mat1,na.rm = TRUE))/(sd(mat1,na.rm = TRUE))
  mat2 <- (mat2 - mean(mat2,na.rm = TRUE))/(sd(mat2,na.rm = TRUE))

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
                         ylim = c(0,max(ncol(mat1),ncol(mat2)))) +
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
                         ylim = c(0,max(ncol(mat1),ncol(mat2)))) +
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

  ccfPlot <- ccfDF %>%
    dplyr::mutate(dx = rev(dx),
                  dy = rev(dy)) %>%
    ggplot2::ggplot(ggplot2::aes(x = dx,y = dy,fill = CCF)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(low = "purple",
                                  mid = "white",
                                  high = "orange",
                                  midpoint = 0) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  ccfMaxSummary <- ccfMaxInfo %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::mutate(dx = -dx,
                  dy = -dy) %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    gridExtra::tableGrob(rows = c("CCFmax","dx","dy"),
                         cols = "CCFmax.Summary")

  layoutMat <- matrix(c(1,2,3,4),ncol = 2,byrow = TRUE)

  gridPlot <- gridExtra::arrangeGrob(mat1Plot,
                                     mat2Plot,
                                     ccfPlot,
                                     ccfMaxSummary,
                                     layout_matrix = layoutMat)

  plot(gridPlot)

  if(returnGrob){
    return(gridPlot)
  }
}

#' Plots the Congruent Matching Cells identified in a cartridge case scan
#'
#' @name cmcPlot
#'
#' @description Plots the selected congruent matching cells of a questioned x3p by either
#'   converting the surface matrix contained in the x3p to a cimg object and
#'   using base plot (fast, but doesn't look particularly great) or by plotting
#'   surface matrix using ggplot2::geom_tile (slow, but nice looking).
#'
#' @param x3p The "questioned" cartridge case scan, in the form of an x3p
#'   object, for which CMCs were determined (this would be x3p1 in a call to
#'   cellCCF)
#' @param cmcDF data frame containing the congruent matching cells of the
#'   questioned x3p object
#' @param method dictates whether a cimg plot is returned, or a ggplot2 plot.
#'   The ggplot2 plot looks nicer, but takes much longer to construct than the
#'   cimg plot.
#'
#' @examples
#' \dontrun{
#' comparison1 <- cellCCF_bothDirections(x3p1,x3p2)
#'
#' cmcs <- cmcFilter_improved(comparison1)
#'
#' #initially selected CMCs based on x3p1 vs x3p2 comparison
#' cmcPlot(x3p1,cmcs$initialCMCs[[1]][[1]])
#'
#' #initialy selected CMCs based on x3p2 vs x3p1
#' cmcPlot(x3p2,cmcs$initialCMCs[[1]][[2]])
#'
#' #final selected CMCs (may be empty) based on both x3p1 vs x3p2 and x3p2 vs x3p1
#' cmcPlot(x3p1,cmcs$finalCMCs)
#'
#' #make plot using ggplot2::geom_tile
#' cmcPlot(x3p1,cmcs$finalCMCs,method = "ggplot2")
#' }
#'
#' @export

cmcPlot <- function(x3p,
                    cmcDF,
                    method = "cimg"){
  blankImage <- matrix(NA,
                       nrow = nrow(x3p$surface.matrix),
                       ncol = ncol(x3p$surface.matrix))

  #these will be the corners of the CMCs in the matrix
  cmcLocs <- cmcDF$cellID %>%
    stringr::str_extract_all(pattern = "[0-9]{1,}") %>%
    purrr::map(as.numeric)

  #add CMCs into the blank image at the same location as they are in the original matrix
  for(ind in 1:length(cmcLocs)){
    cellRowInd <- cmcLocs[[ind]][1]:min(nrow(x3p$surface.matrix),cmcLocs[[ind]][2])
    cellColInd <- cmcLocs[[ind]][3]:min(ncol(x3p$surface.matrix),cmcLocs[[ind]][4])

    blankImage[cellRowInd,cellColInd] <-
      as.matrix(x3p$surface.matrix[cellRowInd,cellColInd],
                nrow = max(cellRowInd),ncol = max(cellColInd))
  }


  #cimg plot is much faster, but I tend to prefer the look of ggplot2
  if(method == "cimg"){
    blankImage %>%
      t() %>%
      imager::as.cimg() %>%
      imager::mirror("y") %>%
      plot()
  }
  if(method == "ggplot2"){
    blankImage %>%
      t() %>%
      imager::as.cimg() %>%
      as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(x = x,y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = value)) +
      ggplot2::scale_fill_gradient2(low = "grey0",
                                    mid = "grey50",
                                    high = "grey100",
                                    na.value = "white") +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.position = "none",
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank())
  }
}

#' Create a bar plot of congruent matching cells per rotation value
#'
#' @name cmcPerThetaBarPlot
#'
#' @param cellCCF_output list returned by the cellCCF or cellCCF_bothDirections
#'   function. If from the cellCCF_bothDirections, then the ggplot will be
#'   faceted by the "direction" of the comparison (i.e., whether x3p1 or x3p2
#'   played the role as the "questioned" cartridge case scan)
#'
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

cmcPerThetaBarPlot <- function(cellCCF_output,
                               consensus_function = median,
                               corr_thresh = .5,
                               dx_thresh = 15,
                               dy_thresh = dx_thresh,
                               theta_thresh = 3,
                               consensus_function_theta = consensus_function,
                               highCMCThresh = 1){
  #make sure that cellCCF_output is either the output of the cellCCF function or
  #cellCCF_bothDirections
  testthat::expect_true(is.list(cellCCF_output))
  testthat::expect_true(identical(names(cellCCF_output),c("comparison_1to2","comparison_2to1")) | identical(names(cellCCF_output),c("params","ccfResults")),
                        label = "cellCCF_output argument is a list returned by cellCCF or cellCCF_bothDirections")

  #if cellCCF_output is from cellCCF_bothDirections:
  if(identical(names(cellCCF_output),c("comparison_1to2","comparison_2to1"))){
    dplyr::bind_rows(
      cellCCF_output$comparison_1to2$ccfResults %>%
        cmcR:::cmcFilterPerTheta(consensus_function = consensus_function,
                                 corr_thresh = corr_thresh,
                                 dx_thresh = dx_thresh,
                                 dy_thresh = dy_thresh,
                                 theta_thresh = theta_thresh) %>%
        dplyr::mutate(comparison = "x3p1 vs. x3p2"),
      cellCCF_output$comparison_2to1$ccfResults %>%
        cmcR:::cmcFilterPerTheta(consensus_function = consensus_function,
                                 corr_thresh = corr_thresh,
                                 dx_thresh = dx_thresh,
                                 dy_thresh = dy_thresh,
                                 theta_thresh = theta_thresh) %>%
        dplyr::mutate(comparison = "x3p2 vs. x3p1")
    ) %>%
      dplyr::group_by(comparison,theta) %>%
      dplyr::tally() %>%
      dplyr::mutate(cmcHigh = max(n) - highCMCThresh)  %>%
      ggplot2::ggplot(ggplot2::aes(x = theta,y = n)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_bw()  +
      ggplot2::xlab(expression(theta*" (degree)")) +
      ggplot2::ylab("CMC number") +
      ggplot2::facet_wrap(~ comparison,ncol = 1)  +
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
  else{
    cellCCF_output$ccfResults  %>%
      cmcR:::cmcFilterPerTheta(consensus_function = consensus_function,
                               corr_thresh = corr_thresh,
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