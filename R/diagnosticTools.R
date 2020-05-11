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
#' @export

ccfMapPlot <- function(mat1,
                       mat2,
                       theta = NA,
                       returnGrob = FALSE,
                       type = "raster"){

  ccfMat <- cmcR:::ccfMap(mat1,mat2)

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
      # ggplot2::geom_contour(aes(z = fft.ccf),
      #                       breaks = quantile(ccfDF$fft.ccf,seq(0,1,length.out = 5)),
      # colour = "black") +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_gradient2(low = "purple",
                                    mid = "white",
                                    high = "orange",
                                    midpoint = 0,aesthetics = c("fill"))+
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank())

    layoutMat <- matrix(c(1,2,4,3),ncol = 2,byrow = TRUE)
  }
  if(type == "contour"){
    ccfPlot <- ccfDF %>%
      dplyr::mutate(dx = rev(dx),
                    dy = rev(dy)) %>%
      ggplot2::ggplot(ggplot2::aes(x = dx,y = dy)) +
      # ggplot2::geom_contour_filled(ggplot2::aes(z = fft.ccf),
      #                              breaks = quantile(as.vector(ccfDF$fft.ccf),probs = seq(0,1,by = .1),na.rm = TRUE)) +
      ggplot2::geom_contour_filled(ggplot2::aes(z = fft.ccf),
                                   breaks = c(quantile(as.vector(ccfDF$fft.ccf)[(as.vector(ccfDF$fft.ccf) <= 0)],
                                                       prob = seq(0,1,length.out = 6)),
                                              0,
                                              quantile(as.vector(ccfDF$fft.ccf)[as.vector((ccfDF$fft.ccf) > 0)],
                                                       prob = seq(0,1,length.out = 6))),
                                   # breaks = c(seq(from = -1*max(abs(c(min(ccfDF$fft.ccf),max(ccfDF$fft.ccf)))),
                                   #                to = 0,length.out = 7),
                                   #            0,
                                   #            seq(from = 0,
                                   #                to = max(abs(c(min(ccfDF$fft.ccf),max(ccfDF$fft.ccf)))),
                                   #                length.out = 7)) %>%
                                   #   unique()
      ) +
      ggplot2::geom_point(data = ccfDF %>%
                            mutate(dx = rev(dx),
                                   dy = rev(dy)) %>%
                            filter(fft.ccf == min(fft.ccf) | fft.ccf == max(fft.ccf)) %>%
                            arrange(fft.ccf) %>%
                            mutate(type = c("Minimum","Maximum")),
                          ggplot2::aes(x = dx,y = dy,shape = type),
                          # shape = 4,
                          colour = "white") +
      ggplot2::scale_shape_manual(values = c(1,4)) +
      ggplot2::coord_fixed(expand = FALSE) +
      # scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
      #                      values = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
      #                      # breaks = function(lims) {quantile(as.vector(fadul1.1$x3p$surface.matrix),c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),na.rm = TRUE)},
      #                      na.value = "grey80") +
      # ggplot2::scale_fill_manual(values = rev(c(rev(RColorBrewer::brewer.pal(5,"Oranges")),"white",RColorBrewer::brewer.pal(5,"Blues"),RColorBrewer::brewer.pal(5,"Purples"))),
      # drop = FALSE) +
      ggplot2::scale_fill_manual(values = colorspace::divergingx_hcl(13,"PuOr"),
                                 drop = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.title = element_text(size = 7),
                     legend.text = element_text(size = 5)) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE,ncol = 3),
                      shape = ggplot2::guide_legend(direction = "horizontal",override.aes = list(colour = c("#312A56","#492900"))))

    layoutMat <- matrix(c(1,2,3,4,4,4),ncol = 3,byrow = TRUE)
  }

  ccfMaxSummary <- ccfMaxInfo %>%
    dplyr::select(-c(x,y)) %>%
    dplyr::mutate(dx = -dx,
                  dy = -dy,
                  theta = theta) %>%
    t() %>% #imager treats a matrix as its transpose ("x" axis in imager refers to rows "y" to cols)
    gridExtra::tableGrob(rows = c("fft.ccfMax","dx","dy","theta"),
                         cols = "Summary")

  gridPlot <- gridExtra::arrangeGrob(mat1Plot,
                                     mat2Plot,
                                     ccfMaxSummary,
                                     ccfPlot,
                                     layout_matrix = layoutMat)

  gridExtra::grid.arrange(gridPlot)

  if(returnGrob){
    return(gridPlot)
  }
}

#' @name arrangeCMCPlot
#'
#' @keywords internal

arrangeCMCPlot <- function(x3p1,
                           x3p2,
                           x3p1_cmcs,
                           x3p1_nonCMCs,
                           x3p2_cmcs,
                           x3p2_nonCMCs,
                           type,
                           directionIndic){

  if(type == "Final"){
    x3p1_cmcs <- x3p1_cmcs %>%
      dplyr::mutate(theta = ifelse(comparison == "comparison_1to2",theta,-theta),
                    dx = ifelse(comparison == "comparison_1to2",dx,-dx),
                    dy = ifelse(comparison == "comparison_1to2",dy,-dy),
                    cmc = rep("yes",times = nrow(.))) %>%
      select(-comparison) %>%
      bind_rows(x3p1_nonCMCs %>%
                  mutate(cmc = rep("no",times = nrow(.)))) %>%
      mutate(cmc = factor(cmc, levels = c("yes","no")))

    x3p2_cmcs <- x3p2_cmcs %>%
      dplyr::mutate(theta = ifelse(comparison == "comparison_1to2",theta,-theta),
                    dx = ifelse(comparison == "comparison_1to2",dx,-dx),
                    dy = ifelse(comparison == "comparison_1to2",dy,-dy),
                    cmc = rep("yes",times = nrow(.))) %>%
      select(-comparison) %>%
      bind_rows(x3p2_nonCMCs %>%
                  mutate(cmc = rep("no",times = nrow(.)))) %>%
      mutate(cmc = factor(cmc, levels = c("yes","no")))

    x3p1_cmcPlot <- x3p1_cmcs %>%
      dplyr::left_join(x3p1_cmcs %>%
                         purrr::pmap_dfr(~ {
                           idNum <- ..7 %>%
                             stringr::str_extract_all(string = ..7,
                                                      pattern = "[0-9]{1,}") %>%
                             unlist() %>%
                             as.numeric()

                           data.frame(cellID = ..7,
                                      firstRow = idNum[1],
                                      lastRow = idNum[2],
                                      firstCol = idNum[3],
                                      lastCol = idNum[4],
                                      stringsAsFactors = FALSE)
                         }),
                       by = "cellID") %>%
      dplyr::arrange(cellNum) %>%
      dplyr::mutate(firstRow = 6.25*(firstRow),
                    lastRow = 6.25*(lastRow),
                    firstCol = 6.25*(firstCol),
                    lastCol = 6.25*(lastCol)) %>%
      dplyr::mutate(midCol = (lastCol + firstCol)/2,
                    midRow = (lastRow + firstRow)/2,
                    cellInd = 1:nrow(.)) %>%
      ggplot2::ggplot() +
      ggplot2::geom_raster(data = {
        tmp <- x3p1

        tmp %>%
          x3ptools::x3p_to_df() %>%
          dplyr::rename(height = value) %>%
          dplyr::mutate(x = 1e6*(x),
                        y = 1e6*(y),
                        height = 1e6*height)},
        ggplot2::aes(x = x,y = y,fill = height),interpolate = TRUE) +
      ggplot2::geom_spoke(ggplot2::aes(x = firstCol,
                                       y = firstRow,angle = 0,
                                       radius = lastRow - firstRow,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = firstCol,
                                       y = firstRow,
                                       angle = pi/2,
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = lastCol,
                                       y = lastRow,angle = pi,
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = lastCol,
                                       y = lastRow,
                                       angle = 3*pi/2,
                                       radius = lastRow - firstRow,
                                       colour = cmc))  +
      ggplot2::scale_colour_manual(values = c("black","red")) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::geom_text(ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = paste0("A",cellInd),
                                      colour = cmc),
                         size = 3) +
      scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           values = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                           # breaks = function(lims) {quantile(as.vector(fadul1.1$x3p$surface.matrix),c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),na.rm = TRUE)},
                           na.value = "grey80") +
      # ggplot2::scale_fill_gradient2(low = "#1a6bff",
      #                               mid = "grey75",
      #                               high = "#ffae1a",
      #                               midpoint = median(as.vector(x3p1$surface.matrix),
      #                                                 na.rm = TRUE),
      #                               na.value = "white") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      # guides(fill = guide_colourbar(barheight = 18)) +
      ggplot2::labs(fill = expression("Height ["*mu*"m]")) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::ylab(expression("Y-Position ["*mu*"m]")) +
      ggplot2::xlab(expression("X-Position ["*mu*"m]")) +
      ggplot2::ggtitle(paste0("x3p1 ",type," CMCs"))

    x3p2_cmcPlot <- x3p2_cmcs  %>%
      dplyr::left_join(x3p2_cmcs %>%
                         purrr::pmap_dfr(~ {
                           idNum <- ..7 %>%
                             stringr::str_extract_all(string = ..7,
                                                      pattern = "[0-9]{1,}") %>%
                             unlist() %>%
                             as.numeric()

                           data.frame(cellID = ..7,
                                      firstRow = idNum[1],
                                      lastRow = idNum[2],
                                      firstCol = idNum[3],
                                      lastCol = idNum[4],stringsAsFactors = FALSE)
                         }),
                       by = "cellID") %>%
      dplyr::arrange(cellNum) %>%
      dplyr::mutate(firstRow = 6.25*(firstRow - dy),
                    lastRow = 6.25*(lastRow - dy),
                    firstCol = 6.25*(firstCol - dx),
                    lastCol = 6.25*(lastCol - dx)) %>%
      dplyr::mutate(firstRowCentered = firstRow - max(lastRow)/2,
                    lastRowCentered = lastRow - max(lastRow)/2,
                    firstColCentered = firstCol - max(lastCol)/2,
                    lastColCentered = lastCol - max(lastCol)/2) %>%
      dplyr::mutate(bottomLeftCorner_col = firstColCentered*cos(theta*(pi/180)) - firstRowCentered*sin(theta*(pi/180)) + max(lastCol)/2,
                    bottomLeftCorner_row = firstColCentered*sin(theta*(pi/180)) + firstRowCentered*cos(theta*(pi/180)) + max(lastRow)/2,
                    topRightCorner_col = lastColCentered*cos(theta*(pi/180)) - lastRowCentered*sin(theta*(pi/180)) + max(lastCol)/2,
                    topRightCorner_row = lastColCentered*sin(theta*(pi/180)) + lastRowCentered*cos(theta*(pi/180)) + max(lastRow)/2) %>%
      dplyr::mutate(midCol = (topRightCorner_col + bottomLeftCorner_col)/2,
                    midRow = (topRightCorner_row + bottomLeftCorner_row)/2,
                    cellInd = 1:nrow(.)) %>%
      ggplot2::ggplot() +
      ggplot2::geom_raster(data = {
        tmp <- x3p2
        tmp$surface.matrix <- cmcR:::rotateSurfaceMatrix(x3p2$surface.matrix,theta = median(x3p2_cmcs$theta))
        tmp %>%
          x3ptools::x3p_to_df() %>%
          dplyr::rename(height = value) %>%
          dplyr::mutate(x = 1e6*(x),
                        y = 1e6*(y),
                        height = 1e6*height)},
        ggplot2::aes(x = x,y = y,fill = height),interpolate = TRUE) +
      ggplot2::geom_spoke(ggplot2::aes(x = bottomLeftCorner_col,
                                       y = bottomLeftCorner_row,
                                       angle = theta*(pi/180),
                                       radius = lastRow - firstRow,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = bottomLeftCorner_col,
                                       y = bottomLeftCorner_row,angle = (pi/2 + theta*(pi/180)),
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = topRightCorner_col,
                                       y = topRightCorner_row,angle = (pi + theta*(pi/180)),
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = topRightCorner_col,
                                       y = topRightCorner_row,
                                       angle = (3*pi/2 + theta*(pi/180)),
                                       radius = lastRow - firstRow,
                                       colour = cmc))  +

      ggplot2::guides(colour = FALSE) +
      ggplot2::geom_text(ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = paste0("B",cellInd),
                                      angle = theta,
                                      colour = cmc),
                         size = 3) +
      ggplot2::scale_colour_manual(values = c("black","red")) +
      scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           values = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                           # breaks = function(lims) {quantile(as.vector(fadul1.1$x3p$surface.matrix),c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),na.rm = TRUE)},
                           na.value = "grey80") +
      # ggplot2::scale_fill_gradient2(low = "#1a6bff",
      #                               mid = "grey75",
      #                               high = "#ffae1a",
      #                               midpoint = median(as.vector(fadul1.2$x3p$surface.matrix),
      #                                                 na.rm = TRUE),
      #                               na.value = "white") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      # guides(fill = guide_colourbar(barheight = 18)) +
      ggplot2::labs(fill = expression("Height ["*mu*"m]")) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::ylab(expression("Y-Position ["*mu*"m]")) +
      ggplot2::xlab(expression("X-Position ["*mu*"m]")) +
      ggplot2::ggtitle(paste0("x3p2 ",type," CMCs"))
  }

  else{
    x3p1_cmcs <- x3p1_cmcs %>%
      dplyr::mutate(theta = c(1,-1)[directionIndic]*theta,
                    dx = c(1,-1)[directionIndic]*dx,
                    dy = c(1,-1)[directionIndic]*dy,
                    cmc = rep("yes",times = nrow(.))) %>%
      bind_rows(x3p1_nonCMCs %>%
                  mutate(cmc = rep("no",times = nrow(.)))) %>%
      mutate(cmc = factor(cmc, levels = c("yes","no")))

    x3p2_cmcs <- x3p2_cmcs %>%
      dplyr::mutate(theta = c(1,-1)[directionIndic]*theta,
                    dx = c(1,-1)[directionIndic]*dx,
                    dy = c(1,-1)[directionIndic]*dy,
                    cmc = rep("yes",times = nrow(.))) %>%
      bind_rows(x3p1_nonCMCs %>%
                  mutate(cmc = rep("no",times = nrow(.)))) %>%
      mutate(cmc = factor(cmc, levels = c("yes","no")))

    x3p1_cmcPlot <- x3p1_cmcs %>%
      dplyr::left_join(x3p1_cmcs %>%
                         purrr::pmap_dfr(~ {
                           idNum <- ..7 %>%
                             stringr::str_extract_all(string = ..7,
                                                      pattern = "[0-9]{1,}") %>%
                             unlist() %>%
                             as.numeric()

                           data.frame(cellID = ..7,
                                      firstRow = idNum[1],
                                      lastRow = idNum[2],
                                      firstCol = idNum[3],
                                      lastCol = idNum[4],stringsAsFactors = FALSE)
                         }),
                       by = "cellID") %>%
      dplyr::arrange(cellNum) %>%
      dplyr::mutate(firstRow = 6.25*(firstRow),
                    lastRow = 6.25*(lastRow),
                    firstCol = 6.25*(firstCol),
                    lastCol = 6.25*(lastCol)) %>%
      dplyr::mutate(midCol = (lastCol + firstCol)/2,
                    midRow = (lastRow + firstRow)/2,
                    cellInd = 1:nrow(.)) %>%
      ggplot2::ggplot() +
      ggplot2::geom_raster(data = {
        tmp <- x3p1

        tmp %>%
          x3ptools::x3p_to_df() %>%
          dplyr::rename(height = value) %>%
          dplyr::mutate(x = 1e6*(x),
                        y = 1e6*(y),
                        height = 1e6*height)},
        ggplot2::aes(x = x,y = y,fill = height),interpolate = TRUE) +
      ggplot2::geom_spoke(ggplot2::aes(x = firstCol,
                                       y = firstRow,
                                       angle = 0,
                                       radius = lastRow - firstRow,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = firstCol,
                                       y = firstRow,
                                       angle = pi/2,
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = lastCol,
                                       y = lastRow,
                                       angle = pi,
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = lastCol,
                                       y = lastRow,
                                       angle = 3*pi/2,
                                       radius = lastRow - firstRow,
                                       colour = cmc))  +
      ggplot2::scale_colour_manual(values = c("black","red")) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::geom_text(ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = paste0("A",cellInd),
                                      colour = cmc),
                         size = 3) +
      scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           values = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                           # breaks = function(lims) {quantile(as.vector(fadul1.1$x3p$surface.matrix),c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),na.rm = TRUE)},
                           na.value = "grey80") +
      # ggplot2::scale_fill_gradient2(low = "#1a6bff",
      #                               mid = "grey75",
      #                               high = "#ffae1a",
      #                               midpoint = median(as.vector(x3p1$surface.matrix),
      #                                                 na.rm = TRUE),
      #                               na.value = "white") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
      # guides(fill = guide_colourbar(barheight = 18)) +
      ggplot2::labs(fill = expression("Height ["*mu*"m]")) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::ylab(expression("Y-Position ["*mu*"m]")) +
      ggplot2::xlab(expression("X-Position ["*mu*"m]")) +
      ggplot2::ggtitle(paste0("x3p1 ",type," CMCs"))

    x3p2_cmcPlot <- x3p2_cmcs  %>%
      dplyr::left_join(x3p2_cmcs %>%
                         purrr::pmap_dfr(~ {
                           idNum <- ..7 %>%
                             stringr::str_extract_all(string = ..7,
                                                      pattern = "[0-9]{1,}") %>%
                             unlist() %>%
                             as.numeric()

                           data.frame(cellID = ..7,
                                      firstRow = idNum[1],
                                      lastRow = idNum[2],
                                      firstCol = idNum[3],
                                      lastCol = idNum[4],stringsAsFactors = FALSE)
                         }),
                       by = "cellID") %>%
      dplyr::arrange(cellID) %>%
      dplyr::mutate(firstRow = 6.25*(firstRow - dy),
                    lastRow = 6.25*(lastRow - dy),
                    firstCol = 6.25*(firstCol - dx),
                    lastCol = 6.25*(lastCol - dx)) %>%
      dplyr::mutate(firstRowCentered = firstRow - max(lastRow)/2,
                    lastRowCentered = lastRow - max(lastRow)/2,
                    firstColCentered = firstCol - max(lastCol)/2,
                    lastColCentered = lastCol - max(lastCol)/2) %>%
      dplyr::mutate(bottomLeftCorner_col = firstColCentered*cos(theta*(pi/180)) - firstRowCentered*sin(theta*(pi/180)) + max(lastCol)/2,
                    bottomLeftCorner_row = firstColCentered*sin(theta*(pi/180)) + firstRowCentered*cos(theta*(pi/180)) + max(lastRow)/2,
                    topRightCorner_col = lastColCentered*cos(theta*(pi/180)) - lastRowCentered*sin(theta*(pi/180)) + max(lastCol)/2,
                    topRightCorner_row = lastColCentered*sin(theta*(pi/180)) + lastRowCentered*cos(theta*(pi/180)) + max(lastRow)/2) %>%
      dplyr::mutate(midCol = (topRightCorner_col + bottomLeftCorner_col)/2,
                    midRow = (topRightCorner_row + bottomLeftCorner_row)/2,
                    cellInd = 1:nrow(.)) %>%
      ggplot2::ggplot() +
      ggplot2::geom_raster(data = {
        tmp <- x3p2
        tmp$surface.matrix <- cmcR:::rotateSurfaceMatrix(x3p2$surface.matrix,theta = median(x3p2_cmcs$theta))
        tmp %>%
          x3ptools::x3p_to_df() %>%
          dplyr::rename(height = value) %>%
          dplyr::mutate(x = 1e6*(x),
                        y = 1e6*(y),
                        height = 1e6*height)},
        ggplot2::aes(x = x,y = y,fill = height),interpolate = TRUE) +
      ggplot2::geom_spoke(ggplot2::aes(x = bottomLeftCorner_col,
                                       y = bottomLeftCorner_row,
                                       angle = theta*(pi/180),
                                       radius = lastRow - firstRow,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = bottomLeftCorner_col,
                                       y = bottomLeftCorner_row,
                                       angle = (pi/2 + theta*(pi/180)),
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = topRightCorner_col,
                                       y = topRightCorner_row,
                                       angle = (pi + theta*(pi/180)),
                                       radius = lastCol - firstCol,
                                       colour = cmc)) +
      ggplot2::geom_spoke(ggplot2::aes(x = topRightCorner_col,
                                       y = topRightCorner_row,
                                       angle = (3*pi/2 + theta*(pi/180)),
                                       radius = lastRow - firstRow,colour = cmc))  +
      ggplot2::scale_colour_manual(values = c("black","red")) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::geom_text(ggplot2::aes(x = midCol,
                                      y = midRow,
                                      label = paste0("B",cellInd),
                                      angle = theta,
                                      colour = cmc),
                         size = 3) +
      scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                           values = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                           # breaks = function(lims) {quantile(as.vector(fadul1.1$x3p$surface.matrix),c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),na.rm = TRUE)},
                           na.value = "grey80") +
      # ggplot2::scale_fill_gradient2(low = "#1a6bff",
      #                               mid = "grey75",
      #                               high = "#ffae1a",
      #                               midpoint = median(as.vector(fadul1.2$x3p$surface.matrix),
      #                                                 na.rm = TRUE),
      #                               na.value = "white") +
      ggplot2::theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      # guides(fill = guide_colourbar(barheight = 18)) +
      ggplot2::labs(fill = expression("Height ["*mu*"m]")) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::ylab(expression("Y-Position ["*mu*"m]")) +
      ggplot2::xlab(expression("X-Position ["*mu*"m]")) +
      ggplot2::ggtitle(paste0("x3p2 ",type," CMCs"))
  }

  gridPlot <- gridExtra::arrangeGrob(x3p1_cmcPlot,
                                     x3p2_cmcPlot,
                                     layout_matrix = matrix(c(1,2),ncol = 2))

  return(gridPlot)
}

#' Visualize initial and final CMCs for a cartridge case pair comparison
#' @name cmcPlot
#'
#' @param x3p1 an x3p object
#' @param x3p2 a different x3p object
#' @param cellCCF_bothDirections_output output from the function cmcR::cellCCF_bothDirections
#' @param cmcFilter_improved_output output from the function cmcR::cmcFilter_improved
#'
#' @export

cmcPlot <- function(x3p1,
                    x3p2,
                    cellCCF_bothDirections_output,
                    cmcFilter_improved_output){

  # get all cellIDs considered in the comparison -- including those not
  # identified as CMCs
  x3p1_topResults <- cellCCF_bothDirections_output$comparison_1to2$ccfResults %>%
    cmcR::topResultsPerCell()

  x3p2_topResults <- cellCCF_bothDirections_output$comparison_2to1$ccfResults %>%
    cmcR::topResultsPerCell()

  x3p1_cellID <- x3p1_topResults %>%
    dplyr::select(cellNum,cellID)

  x3p2_cellID <- x3p2_topResults %>%
    dplyr::select(cellNum,cellID)

  # the initial CMCs assigned to particular pair are the minimum CMCs computed
  # in either "direction", so we'll need to keep track of the direction
  directionIndic <- which.min(c(nrow(cmcFilter_improved_output$initialCMC[[1]]),
                                nrow(cmcFilter_improved_output$initialCMC[[2]])))

  # determine which cellIDs are identified as initial CMCs
  x3p1_initialCMC <- cmcFilter_improved_output$initialCMC[[directionIndic]] %>%
    dplyr::ungroup() %>%
    dplyr::select(-cellID) %>%
    dplyr::inner_join(x3p1_cellID,by = "cellNum")

  x3p2_initialCMC <- cmcFilter_improved_output$initialCMC[[directionIndic]] %>%
    dplyr::ungroup() %>%
    dplyr::select(-cellID) %>%
    dplyr::inner_join(x3p2_cellID,by = "cellNum")

  #cells not identified as CMCs will be drawn as red squares
  x3p1_nonInitialCMC <- x3p1_cellID %>%
    ungroup() %>%
    anti_join(cmcFilter_improved_output$initialCMC[[directionIndic]] %>%
                dplyr::ungroup() %>%
                dplyr::select(-cellID),
              by = "cellNum") %>%
    select(-cellID) %>%
    inner_join(cellCCF_bothDirections_output[[1]] %>%
                 .$ccfResults %>%
                 cmcR::topResultsPerCell(),
               by = "cellNum")

  x3p2_nonInitialCMC <- x3p2_cellID %>%
    ungroup() %>%
    anti_join(cmcFilter_improved_output$initialCMC[[directionIndic]] %>%
                dplyr::ungroup() %>%
                dplyr::select(-cellID),
              by = "cellNum") %>%
    select(-cellID) %>%
    inner_join(cellCCF_bothDirections_output[[2]] %>%
                 .$ccfResults %>%
                 cmcR::topResultsPerCell(),
               by = "cellNum")

  # creates plots of the initial CMCs on both x3ps
  initialPlt <- arrangeCMCPlot(x3p1 = x3p1,
                               x3p2 = x3p2,
                               x3p1_cmcs = x3p1_initialCMC,
                               x3p1_nonCMCs = x3p1_nonInitialCMC,
                               x3p2_cmcs = x3p2_initialCMC,
                               x3p2_nonCMCs = x3p2_nonInitialCMC,
                               type = "Initial",
                               directionIndic = directionIndic)

  if(!purrr::is_empty(cmcFilter_improved_output$finalCMCs)){
    x3p1_finalCMC <- cmcFilter_improved_output$finalCMCs %>%
      dplyr::select(-cellID) %>%
      dplyr::inner_join(x3p1_cellID,by = "cellNum")

    x3p2_finalCMC <- cmcFilter_improved_output$finalCMCs %>%
      dplyr::select(-cellID) %>%
      dplyr::inner_join(x3p2_cellID,by = "cellNum")

    x3p1_nonFinalCMC <- x3p1_cellID %>%
      ungroup() %>%
      anti_join(cmcFilter_improved_output$finalCMCs %>%
                  dplyr::select(-cellID),
                by = "cellNum") %>%
      select(-cellID) %>%
      inner_join(cellCCF_bothDirections_output[[1]] %>%
                   .$ccfResults %>%
                   cmcR::topResultsPerCell(),
                 by = "cellNum")

    x3p2_nonFinalCMC <- x3p2_cellID %>%
      ungroup() %>%
      anti_join(cmcFilter_improved_output$finalCMCs %>%
                  dplyr::select(-cellID),
                by = "cellNum") %>%
      select(-cellID) %>%
      inner_join(cellCCF_bothDirections_output[[2]] %>%
                   .$ccfResults %>%
                   cmcR::topResultsPerCell(),
                 by = "cellNum")

    finalPlt <- cmcR:::arrangeCMCPlot(x3p1 = x3p1,
                                      x3p2 = x3p2,
                                      x3p1_cmcs = x3p1_finalCMC,
                                      x3p1_nonCMCs = x3p1_nonFinalCMC,
                                      x3p2_cmcs = x3p2_finalCMC,
                                      x3p2_nonCMCs = x3p2_nonFinalCMC,
                                      type = "Final")

    return(list("initialCMC" = initialPlt,
                "finalCMC" = finalPlt))
  }

  return(list("initialCMC" = initialPlt))
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
                               ccf_thresh = .6,
                               dx_thresh = 10,
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
                                 ccf_thresh = ccf_thresh,
                                 dx_thresh = dx_thresh,
                                 dy_thresh = dy_thresh,
                                 theta_thresh = theta_thresh) %>%
        dplyr::mutate(comparison = "x3p1 vs. x3p2"),
      cellCCF_output$comparison_2to1$ccfResults %>%
        cmcR:::cmcFilterPerTheta(consensus_function = consensus_function,
                                 ccf_thresh = ccf_thresh,
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

#' @name getCellRegionPairs
#'
#' @description This function is meant as a diagnostic tool to determine what
#'   the cell region pairs for a particular cartridge case comparison looked
#'   like right before calculation of the CCF.
#'
#' @export

getCellRegionPairs <- function(x3p1,x3p2,ccfDF,params){
  mat1 <- x3p1$surface.matrix
  mat2 <- x3p2$surface.matrix

  if(is.null(params$centerCell)){
    m1 <- 0
    m2 <- 0
  }
  if(is.null(params$scaleCell)){
    sd1 <- 1
    sd2 <- 1
  }

  if(!is.null(params$centerCell)){
    if(params$centerCell == "wholeMatrix"){
      m1 <- mean(as.vector(mat1),na.rm = TRUE)

      m2 <- mean(as.vector(mat2),na.rm = TRUE)
    }
  }

  if(!is.null(params$scaleCell)){
    if(params$scaleCell == "wholeMatrix"){
      sd1 <- sd(as.vector(mat1),na.rm = TRUE)

      sd2 <- sd(as.vector(mat2),na.rm = TRUE)
    }
  }

  mat1_split <- splitSurfaceMat1(surfaceMat = mat1,
                                 cellNumHoriz = params$cellNumHoriz,
                                 cellNumVert = params$cellNumVert,
                                 minObservedProp = params$minObservedProp)

  sidelengthMultiplier <- floor(sqrt(params$regionToCellProp))

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

    if(!is.null(params$centerCell)){
      if(params$centerCell == "individualCell"){
        m1 <- mat1_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)



        m2 <-  mat2_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    if(!is.null(params$scaleCell)){
      if(params$scaleCell == "individualCell"){
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
