# Cuts out a cell in a surface matrix.
# @name extractCellbyCornerLocs
#
# @param cornerLocs the location (indices) of the desired cell in the image
# @param rotatedSurfaceMat a surface matrix representing a rotated breech face
#   impression image
#
# @description This is a helper for the cellCCF function. This function
#   cuts out a cell in a surface matrix based on the most top, bottom, left,
#   and right indices of the cell's location in the surface matrix.
#
# @keywords internal

extractCellbyCornerLocs <- function(cornerLocs,
                                    rotatedSurfaceMat,
                                    mat2Dim){
  #perform the appropriate subsetting of image A to create a list of larger
  #cells than those in image B. There's a good chance that the
  #splitRotatedSurfaceMat object will need to be padded in various ways, which
  #is done in the rest of the function. Padding occurs when the cornerLocs
  #values extend past the dimension of the surface matrix - either
  #cornerLocs[top/left] <= 0 or cornerLocs[bottom/right] > nrow/ncol(target). In
  #the function getMat2SplitIndices, the important index per target region is
  #the target index in the center of the region - so to make sure that this
  #index stays in the center, we need to perform some padding of the region
  splitRotatedSurfaceMat <- rotatedSurfaceMat[max(cornerLocs[["top"]],1):min(cornerLocs[["bottom"]],mat2Dim[1]),
                                              max(cornerLocs[["left"]],1):min(cornerLocs[["right"]],mat2Dim[2])]


  # pad the rows, as needed
  if(cornerLocs[["top"]] <= 0){
    # e.g., if top == -1, say. Then we need to pad the region 2 rows on the
    # top so that the "new" top index starts at (positive) 1.
    rowPad <- matrix(NA,nrow = abs(cornerLocs[["top"]]) + 1,ncol = ncol(splitRotatedSurfaceMat))
    splitRotatedSurfaceMat <- rbind(rowPad,splitRotatedSurfaceMat)
  }
  # e.g., if bottom == 101 and nrow(target) == 100, then we need to pad the
  # bottom of the region by 1 row
  if(cornerLocs[["bottom"]] > mat2Dim[1]){
    rowPad <- matrix(NA,nrow = cornerLocs[["bottom"]] - mat2Dim[1],ncol = ncol(splitRotatedSurfaceMat))
    splitRotatedSurfaceMat <- rbind(splitRotatedSurfaceMat,rowPad)
  }
  if(cornerLocs[["left"]] <= 0){
    colPad <- matrix(NA,nrow = nrow(splitRotatedSurfaceMat),ncol = abs(cornerLocs[["left"]]) + 1)
    splitRotatedSurfaceMat <- cbind(colPad,splitRotatedSurfaceMat)
  }
  if(cornerLocs[["right"]] > mat2Dim[2]){
    colPad <- matrix(NA,nrow = nrow(splitRotatedSurfaceMat),ncol = cornerLocs[["right"]] - mat2Dim[2])
    splitRotatedSurfaceMat <- cbind(splitRotatedSurfaceMat,colPad)
  }
  #   if(nrow(splitRotatedSurfaceMat) != ncol(splitRotatedSurfaceMat)){ #if the matrix isn't square...
  #
  #     #if the rows need padding...
  #     if(nrow(splitRotatedSurfaceMat) < max(dim(splitRotatedSurfaceMat))){
  #
  #       rowsToPad <- ncol(splitRotatedSurfaceMat) - nrow(splitRotatedSurfaceMat)
  #       rowPadder <- matrix(NA,nrow = rowsToPad,ncol = ncol(splitRotatedSurfaceMat))
  #
  #       #if the split comes from the top of the overall matrix...
  #       if(cornerLocs[["top"]] == 1){
  #         splitRotatedSurfaceMat <- rbind(rowPadder,
  #                                         splitRotatedSurfaceMat)
  #       }
  #
  #       #if the split comes from the bottom of the overall matrix....
  #       if(cornerLocs[["bottom"]] == mat2Dim[1]){
  #         splitRotatedSurfaceMat <- rbind(splitRotatedSurfaceMat,
  #                                         rowPadder)
  #       }
  #     }
  #
  #     #if the cols need padding...
  #     if(ncol(splitRotatedSurfaceMat) < max(dim(splitRotatedSurfaceMat))){
  #
  #       colsToPad <- nrow(splitRotatedSurfaceMat) - ncol(splitRotatedSurfaceMat)
  #       colPadder <- matrix(NA,ncol = colsToPad,nrow = nrow(splitRotatedSurfaceMat))
  #
  #       #if the split comes from the left side of the overall matrix...
  #       if(cornerLocs[["left"]] == 1){
  #         splitRotatedSurfaceMat <- cbind(colPadder,
  #                                         splitRotatedSurfaceMat)
  #       }
  #       #if the split comes from the right side of the overall matrix....
  #       if(cornerLocs[["right"]] == mat2Dim[2]){
  #         splitRotatedSurfaceMat <- cbind(splitRotatedSurfaceMat,
  #                                         colPadder)
  #       }
  #     }
  #   }

  return(splitRotatedSurfaceMat)
}

# @name rotateSurfaceMatrix
#
# @keywords internal
# @importFrom rlang .data

utils::globalVariables(".")

rotateSurfaceMatrix <- function(surfaceMat,
                                theta = 0,
                                interpolation = 0){
  surfaceMatFake <- (surfaceMat*10^5) + 1 #scale and shift all non-NA pixels up 1 (meter)
  # imFakeRotated <- :bilinearInterpolation(imFake,theta)
  surfaceMatFakeRotated <- surfaceMatFake %>%
    imager::as.cimg() %>%
    imager::imrotate(angle = theta,
                     interpolation = interpolation, #linear interpolation,
                     cx = floor(nrow(.)/2), #imager treats the rows as the "x" axis of an image
                     cy = floor(ncol(.)/2),
                     boundary = 0) %>% #pad boundary with 0s (dirichlet condition)
    as.matrix()

  surfaceMatFakeRotated[surfaceMatFakeRotated == 0] <- NA
  #shift all of the legitimate pixels back down by 1:
  surfaceMatRotated <- (surfaceMatFakeRotated - 1)/(10^5)

  return(surfaceMatRotated)
}

# @name getMat2SplitIndices
#
# @keywords internal
#
# @importFrom stats setNames
# @importFrom rlang .data

getMat2SplitIndices <- function(cellRanges,
                                cellSideLengths,
                                mat2Dim,
                                sidelengthMultiplier,
                                ...){
  mat2_splitCorners <- cellRanges %>%
    #pull all numbers from cellRange strings:
    purrr::map(~ stringr::str_extract_all(string = .,pattern = "[0-9]{1,}")) %>%
    purrr::map(~ c(
      #y-position of each cell's center:
      "y" = mean(c(as.numeric(.[[1]][[3]]),as.numeric(.[[1]][[4]]))),
      #x-position of each cell's center:
      "x" = mean(c(as.numeric(.[[1]][[1]]),as.numeric(.[[1]][[2]]))))) %>%
    #determine the indices of a larger cell to search in image B
    purrr::map2(.x = .,
                .y = cellSideLengths,
                function(xyLoc,sideLength){
                  expandedCellCorners <-
                    c(floor(xyLoc["y"] - sidelengthMultiplier*sideLength["col"]/2),
                      ceiling(xyLoc["y"] + sidelengthMultiplier*sideLength["col"]/2),
                      floor(xyLoc["x"] - sidelengthMultiplier*sideLength["row"]/2),
                      ceiling(xyLoc["x"] + sidelengthMultiplier*sideLength["row"]/2)) %>%
                    setNames(c("left","right","top","bottom"))

                  #replace negative indices with 1 (left/upper-most cells):
                  # expandedCellCorners[expandedCellCorners <= 0] <- 1
                  #replace indices greater than the maximum index with the
                  #maximum index (right/bottom-most cells): Note that imager
                  #treats the rows of a matrix as the "x" axis and the columns
                  #as the "y" axis, contrary to intuition - effectively treating
                  #a matrix as its transpose. As such, we need
                  #to swap the dimensions for when we subset the image further
                  #down in the function
                  # if(expandedCellCorners[c("right")] > mat2Dim[2]){
                  #   expandedCellCorners[c("right")] <- mat2Dim[2]
                  # }
                  # if(expandedCellCorners[c("bottom")] > mat2Dim[1]){
                  #   expandedCellCorners[c("bottom")] <- mat2Dim[1]
                  # }
                  # }

                  return(expandedCellCorners)
                }) %>%
    setNames(cellRanges)

  return(mat2_splitCorners)
}

# @name swapcellRangeAxes
#
# @keywords internal

swapcellRangeAxes <- function(cellRange){
  sSplit <- stringr::str_split(string = cellRange,pattern = ",",n = 2)

  paste0(stringr::str_replace(string = sSplit[[1]][1],pattern = "x =",replacement = "rows:"),
         ", ",
         stringr::str_replace(string = sSplit[[1]][2],pattern = "y =",replacement = "cols:"))
}

#' Split a reference scan into a grid of cells
#'
#' @name comparison_cellDivision
#'
#' @param x3p an x3p object containing a breech face scan
#' @param numCells number of cells to partition the breech face scan into. Must be
#'   a perfect square (49, 64, 81, etc.)
#'
#' @return A tibble containing a numCells number of rows. Each row contains a
#'   single cell's index of the form (row #, col #) and an x3p object containing
#'   the breech face scan of that cell.
#'
#'@examples
#'data(fadul1.1_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64)
#'
#'head(cellTibble)
#'
#'@importFrom rlang .data
#' @export

utils::globalVariables(c(".x",".y"))

comparison_cellDivision <- function(x3p,numCells = 64){

  # Future note: put `cellRange` column information in each x3p's metadata?
  # Would then need to change the getTargetRegions function below to look in the
  # metadata rather than expecting the column cellRange

  assertthat::are_equal(sqrt(numCells) %% 1, 0)

  splitSurfaceMat <- x3p$surface.matrix %>%
    imager::as.cimg() %>%
    imager::imsplit(axis = "x",
                    nb = sqrt(numCells)) %>%
    purrr::map(.f = ~ imager::imsplit(.x,
                                      axis = "y",
                                      nb = sqrt(numCells))) %>%
    purrr::map_depth(.depth = 2,
                     .f = ~ purrr::set_names(as.matrix(.),NULL))

  cellRanges <- purrr::map(names(splitSurfaceMat),
                           function(horizCell){
                             purrr::map(.x = names(splitSurfaceMat[[1]]),
                                        function(vertCell) swapcellRangeAxes(paste(horizCell,vertCell,sep = ",")))
                           }) %>%
    unlist()

  splitSurfaceMat <- splitSurfaceMat %>%
    purrr::flatten() %>%
    purrr::map2(.x = .,
                .y = cellRanges,
                function(cellMatrix = .x,cellRange = .y){
                  cell_x3p <- x3ptools::df_to_x3p(data.frame(x = 1,y = 1,value = NA))

                  cell_x3p$surface.matrix <- cellMatrix

                  #update metainformation
                  cell_x3p$header.info <- x3p$header.info
                  cell_x3p$header.info$sizeY <- ncol(cellMatrix)
                  cell_x3p$header.info$sizeX <- nrow(cellMatrix)

                  #include which rows/columns in the original scan each cell was
                  #taken from
                  cell_x3p$cmcR.info$cellRange <- cellRange

                  return(cell_x3p)
                })

  cellTibble <- tibble::tibble(cellNum = 1:numCells,
                               cellHeightValues = splitSurfaceMat) %>%
    dplyr::mutate(cellIndex = linear_to_matrix(.data$cellNum,nrow = ceiling(sqrt(numCells)))) %>%
    # the assignment of cellIndex below will assign cell based on how the old
    # x3pListPlot function plotted scans (i.e., cell 1,1 will actually be in the
    # top-left corner when using the old x3pListPlot). To make things more
    # intuitive, we choose instead to assign cellIndex based on the original
    # surface matrix - meaning cell 1,1 will actually consist of the first
    # rows/cols of the original surface matrix. The new version of the various
    # plotting functions are also changed to reflect this change.
    # dplyr::mutate(cellIndex = linear_to_matrix(index = (.data$cellNum %% ceiling(sqrt(max(.data$cellNum)))) +
    #                                              floor((ceiling(sqrt(max(.data$cellNum)))^2 - .data$cellNum)/ceiling(sqrt(max(.data$cellNum))))*ceiling(sqrt(max(.data$cellNum))) +
    #                                              ifelse(.data$cellNum %% ceiling(sqrt(max(.data$cellNum))) == 0,ceiling(sqrt(max(.data$cellNum))),0),
    #                                            nrow = ceiling(sqrt(max(.data$cellNum))),
    #                                            byrow = TRUE)) %>%
    dplyr::select(.data$cellIndex,.data$cellHeightValues)

  return(cellTibble)
}

#' Calculate the proportion of missing values in a breech face scan
#'
#' @name comparison_calcPropMissing
#' @param heightValues list/tibble column of x3p objects
#' @return a vector of the same length as the input containing the proportion of
#'   missing values in each x3p object's breech face scan.
#'@examples
#'data(fadul1.1_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(cellPropMissing = comparison_calcPropMissing(heightValues = cellHeightValues))
#'
#'head(cellTibble)
#'
#'@importFrom rlang .data
#' @export

comparison_calcPropMissing <- function(heightValues){
  heightValues %>%
    purrr::map_dbl(~ sum(is.na(.$surface.matrix))/length(.$surface.matrix))
}

#'Extract regions from a target scan based on associated cells in reference scan
#'
#'@name comparison_getTargetRegions
#'@param cellHeightValues list/tibble column of x3p objects containing a
#'  reference scan's cells (as returned by comparison_cellDivision)
#'@param target x3p object containing a breech face scan to be compared to the
#'  reference cell.
#'@param theta degrees that the target scan is to be rotated prior extracting
#'  regions.
#'@param regionSizeMultiplier ratio between the area of each target scan regions
#'  and the reference scan cells (e.g., 9 means that the regions' surface
#'  matrices will have thrice the number of rows and columns as the cells'
#'  surface matrices, 4 means twice the number rows and columns, etc.)
#'@return A list of the same length as the input containing x3p objects from the
#'  target scan.
#'@examples
#'
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = cellHeightValues,
#'                                                                target = fadul1.2_processed)) %>%
#' dplyr::mutate(cellPropMissing = comparison_calcPropMissing(heightValues = cellHeightValues),
#'               regionPropMissing = comparison_calcPropMissing(heightValues = regionHeightValues)) %>%
#'dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85)
#'
#'head(cellTibble)
#'
#'@export

comparison_getTargetRegions <- function(cellHeightValues,
                                        target,
                                        theta = 0,
                                        regionSizeMultiplier = 9){

  cellSideLengths <- cellHeightValues %>%
    purrr::map(~ c("row" = nrow(.$surface.matrix),
                   "col" = ncol(.$surface.matrix)))

  cellRange <- cellHeightValues %>%
    purrr::map_chr(~ .$cmcR.info$cellRange)

  target_regionIndices <- getMat2SplitIndices(cellRanges = cellRange,
                                              cellSideLengths = cellSideLengths,
                                              mat2Dim = dim(target$surface.matrix),
                                              sidelengthMultiplier = floor(sqrt(regionSizeMultiplier)))

  target_surfaceMat_rotated <- rotateSurfaceMatrix(target$surface.matrix,
                                                   theta = theta)


  target_splitRotated <-
    purrr::map(.x = target_regionIndices,
               function(cornerIndices){
                 if(cornerIndices["left"] < cornerIndices["right"] & cornerIndices["top"] < cornerIndices["bottom"]){
                   regionMatrix <- extractCellbyCornerLocs(cornerLocs = cornerIndices,
                                                           rotatedSurfaceMat = target_surfaceMat_rotated,
                                                           mat2Dim = dim(target$surface.matrix))
                 }
                 else{
                   regionMatrix <- matrix(NA)
                 }
                 region_x3p <- x3ptools::df_to_x3p(data.frame(x = 1,y = 1,value = NA))
                 region_x3p$surface.matrix <- regionMatrix
                 #update metainformation
                 region_x3p$header.info <- target$header.info
                 region_x3p$header.info$sizeY <- ncol(regionMatrix)
                 region_x3p$header.info$sizeX <- nrow(regionMatrix)
                 region_x3p$cmcR.info$regionIndices <- cornerIndices %>%
                   purrr::set_names(c("colStart","colEnd","rowStart","rowEnd"))
                 return(region_x3p)
               })

  return(target_splitRotated)
}

#'Standardize height values of a scan by centering/scaling by desired statistics
#'and replacing missing values
#'
#'@name comparison_standardizeHeights
#'
#'@param heightValues list/tibble column of x3p objects
#'@param withRespectTo either "individualCell", meaning centering/scaling
#'  statistics are calculated independently per cell, or "entireScan", meaning
#'  centering/scaling statistics are calculated aggregate across all cells
#'@param centerBy statistic by which to center (i.e., subtract from) the height
#'  values
#'@param scaleBy statistic by which to scale (i.e., divide) the height values
#'
#'@return A list of the same length as the input containing x3p objects with
#'  standardized surface matrices
#'
#'@note this function adds information to the metainformation of the x3p scan it
#'  is given that is required for calculating, for example, the
#'  pairwise-complete correlation using the comparison_cor function.
#'@examples
#'
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = cellHeightValues,
#'                                                                target = fadul1.2_processed)) %>%
#' dplyr::mutate(cellPropMissing = comparison_calcPropMissing(heightValues = cellHeightValues),
#'               regionPropMissing = comparison_calcPropMissing(heightValues = regionHeightValues)) %>%
#'dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
#'dplyr::mutate(cellHeightValues = comparison_standardizeHeights(heightValues = cellHeightValues),
#'              regionHeightValues = comparison_standardizeHeights(heightValues = regionHeightValues))
#'
#'head(cellTibble)
#'
#'@importFrom stats sd
#'@export

comparison_standardizeHeights <- function(heightValues,
                                          withRespectTo = "individualCell",
                                          centerBy = mean,
                                          scaleBy = sd){

  if(withRespectTo == "entireScan"){
    centerByVal <- heightValues %>%
      purrr::map(function(x3p){

        as.vector(x3p$surface.matrix[!is.na(x3p$surface.matrix)])

      }) %>%
      unlist() %>%
      centerBy()

    scaleByVal <- heightValues %>%
      purrr::map(function(x3p){

        as.vector(x3p$surface.matrix[!is.na(x3p$surface.matrix)])

      }) %>%
      unlist() %>%
      scaleBy()

    heightValues <- heightValues %>%
      purrr::map(function(x3p){

        x3p$cmcR.info$centerBy <- centerBy
        x3p$cmcR.info$centerByVal <- centerByVal

        x3p$cmcR.info$scaleBy <- scaleBy
        x3p$cmcR.info$scaleByVal <- scaleByVal

        x3p$surface.matrix <- (x3p$surface.matrix - centerByVal)/scaleByVal

        return(x3p)
      })

    return(heightValues)
  }
  else if(withRespectTo == "individualCell"){

    heightValues <- heightValues %>%
      purrr::map(function(x3p){

        x3p$cmcR.info$centerBy <- centerBy
        x3p$cmcR.info$centerByVal <- centerBy(x3p$surface.matrix,na.rm = TRUE)

        x3p$cmcR.info$scaleBy <- scaleBy
        x3p$cmcR.info$scaleByVal <- scaleBy(x3p$surface.matrix,na.rm = TRUE)

        x3p$surface.matrix <- (x3p$surface.matrix - centerBy(x3p$surface.matrix,na.rm = TRUE))/scaleBy(x3p$surface.matrix,na.rm = TRUE)

        return(x3p)
      })

    return(heightValues)
  }
  else{
    return(heightValues)
  }
}

#'Replace missing values in a scan
#'
#'@name comparison_replaceMissing
#'@param heightValues list/tibble column of x3p objects
#'@param replacement value to replace NAs
#'@return A list of the same length as the input containing x3p objects for
#'  which NA values have been replaced.
#'@examples
#'
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(regionHeightValues =
#'              comparison_getTargetRegions(cellHeightValues = cellHeightValues,
#'                                          target = fadul1.2_processed)) %>%
#' dplyr::mutate(cellPropMissing =
#'                  comparison_calcPropMissing(heightValues = cellHeightValues),
#'               regionPropMissing =
#'            comparison_calcPropMissing(heightValues = regionHeightValues)) %>%
#'dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
#'dplyr::mutate(cellHeightValues =
#'               comparison_standardizeHeights(heightValues = cellHeightValues),
#'              regionHeightValues =
#'         comparison_standardizeHeights(heightValues = regionHeightValues)) %>%
#'dplyr::mutate(cellHeightValues =
#'                   comparison_replaceMissing(heightValues = cellHeightValues),
#'              regionHeightValues =
#'                 comparison_replaceMissing(heightValues = regionHeightValues))
#'
#'head(cellTibble)
#'
#'@export

comparison_replaceMissing <- function(heightValues,
                                      replacement = 0){
  replacedHeights <- heightValues %>%
    purrr::map(function(x3p){
      x3p$surface.matrix[is.na(x3p$surface.matrix)] <- 0

      return(x3p)
    })

  return(replacedHeights)
}

#'Estimate translation alignment between a cell/region pair based on the
#'Cross-Correlation Theorem.
#'
#'@name comparison_fft_ccf
#'@param cellHeightValues list/tibble column of x3p objects containing a
#'  reference scan's cells (as returned by comparison_cellDivision)
#'@param regionHeightValues list/tibble column of x3p objects containing a
#'  target scan's regions (as returned by comparison_getTargetRegions)
#'@param ccfMethod algorithm to calculate the cross-correlation. Either "fft",
#'  which uses the Cross-Correlation Theorem (see link below), or "imager",
#'  which uses the imager::correlate function
#'@return A list of the same length as the input containing data frames of the
#'  translation (x,y) values at which each reference cell is estimated to align
#'  in its associated target region and the CCF value at this alignment.
#'@note The FFT is not defined for matrices containing missing values. The
#'  missing values in the cell and region need to be replaced before using this
#'  function. See the \link[cmcR]{comparison_replaceMissing} function to replace
#'  missing values after standardization.
#'@seealso \url{https://mathworld.wolfram.com/Cross-CorrelationTheorem.html}
#'
#'@return a data frame containing the translation (x,y) at which the CCF was
#'  maximized in aligning a target scan region to its associated reference scan
#'  cell.
#'
#'@examples
#'
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(regionHeightValues =
#'              comparison_getTargetRegions(cellHeightValues = cellHeightValues,
#'                                          target = fadul1.2_processed)) %>%
#' dplyr::mutate(cellPropMissing =
#'            comparison_calcPropMissing(heightValues = cellHeightValues),
#'               regionPropMissing =
#'            comparison_calcPropMissing(heightValues = regionHeightValues)) %>%
#'dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
#'dplyr::mutate(cellHeightValues =
#'         comparison_standardizeHeights(heightValues = cellHeightValues),
#'              regionHeightValues =
#'         comparison_standardizeHeights(heightValues = regionHeightValues)) %>%
#'dplyr::mutate(cellHeightValues =
#'                   comparison_replaceMissing(heightValues = cellHeightValues),
#'              regionHeightValues =
#'             comparison_replaceMissing(heightValues = regionHeightValues)) %>%
#'dplyr::mutate(fft_ccf_df = comparison_fft_ccf(cellHeightValues,
#'                                              regionHeightValues))
#'
#'cellTibble %>%
#' tidyr::unnest(cols = fft_ccf_df) %>%
#' head()
#'
#'@export
comparison_fft_ccf <- function(cellHeightValues,regionHeightValues,ccfMethod = "fft"){
  ccfList <- purrr::map2(cellHeightValues,
                         regionHeightValues,
                         ~ ccfComparison(mat1 = .x$surface.matrix,
                                         mat2 = .y$surface.matrix,
                                         ccfMethod = ccfMethod))

  return(ccfList)
}


# This is a helper function for the comparison_alignedTargetCell function. Given
# a reference cell, target region, and estimated translation registration values
# (dx,dy), extracts a matrix from the target region of the same dimension as the
# reference cell that the estimated translation indicates has the highest
# similarity to the reference cell.
extractTargetCell <- function(cell,
                              region,
                              dx = 0,
                              dy = 0,
                              m1 = 0,
                              m2 = 0,
                              sd1 = 1,
                              sd2 = 1){
  # this isn't necessary because what is passed to comparison_alignedTargetCell
  # are the cell/region already standardized:

  #cell <- (cell - m1)/sd1
  # region <- (region - m2)/sd2

  # cell <- cell*sd1 + m1
  # region <- region*sd2 + m2

  #pad the region as it was when the CCF was calculated using the FFT method:
  regionPadded <- matrix(NA,
                         nrow = nrow(cell) + nrow(region) - 1,
                         ncol = ncol(cell) + ncol(region) - 1)
  regionPadded[1:nrow(region),1:ncol(region)] <- region

  regionCenter <- floor(dim(regionPadded)/2)

  alignedCenter <- regionCenter - c(dy + floor(nrow(cell)/2) + 1,dx + floor(ncol(cell)/2) + 1)

  alignedRows <- c((alignedCenter[1] - floor(dim(cell)[1]/2)),(alignedCenter[1] + floor(dim(cell)[1]/2)))

  alignedCols <- c((alignedCenter[2] - floor(dim(cell)[2]/2)),(alignedCenter[2] + floor(dim(cell)[2]/2)))

  regionCroppedInitial <- regionPadded[max(1,alignedRows[1]):min(nrow(region),alignedRows[2]),
                                       max(1,alignedCols[1]):min(ncol(region),alignedCols[2])]

  # Important Note: after looking at how the old version of this function (using
  # regionCroppedInitial as opposed to regionCroppedShifted) behaved on
  # self-comparisons, it seems that the estimated indices were off by like 2
  # rows and columns. I'm going to change it in this function to the
  # shift-corrected version (and continue experimenting with self-comparisons to
  # make sure this works correctly) and keep the original cmcR implementation
  # intact.

  regionCroppedShifted <- regionPadded[max(1,(alignedRows[1] + 2)):min(nrow(region),(alignedRows[2] + 2)),
                                       max(1,(alignedCols[1] + 2)):min(ncol(region),(alignedCols[2] + 2))]

  #build-in some contingencies in case the regionCropped isn't same size as
  #the cell. If the dimensions don't agree, then we need to consider copies
  #of regionCroppedInitial that have had rows/cols (whichever
  #applicable) cropped from the top/bottom or left/right
  regionCroppedList <- list("initial" = regionCroppedInitial,
                            "shifted" = regionCroppedShifted)

  if(all(dim(regionCroppedShifted) == dim(cell))){

    return(list("regionCropped" = regionCroppedShifted,
                "centerRow" = alignedCenter[1] + 2,
                "centerCol" = alignedCenter[2] + 2,
                "regionRows" = c(max(1,(alignedRows[1] + 2)),min(nrow(region),(alignedRows[2] + 2))),
                "regionCols" = c(max(1,(alignedCols[1] + 2)),min(ncol(region),(alignedCols[2] + 2)))
                ))
                # "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                # "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

  }
  if(any(dim(regionCroppedShifted) != dim(cell))){
    #these make sure that the indices that we've cropped the region by aren't
    #less than 1 or larger than the dimension of the region
    if((alignedRows[1] + 2) < 1){
      rowsToPad <- -1*(alignedRows[1] + 2) + 1

      regionCroppedShifted <- rbind(matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedShifted)),
                                    regionCroppedShifted)
    }
    if((alignedRows[2] + 2) > nrow(region)){
      rowsToPad <- (alignedRows[2] + 2) - nrow(region)

      regionCroppedShifted <- rbind(regionCroppedShifted,
                                    matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedShifted)))

    }
    if((alignedCols[1] + 2) < 1){
      colsToPad <- -1*(alignedCols[1] + 2) + 1

      regionCroppedShifted <- cbind(matrix(NA,nrow = nrow(regionCroppedShifted),ncol = colsToPad),
                                    regionCroppedShifted)

    }
    if((alignedCols[2] + 2) > ncol(region)){
      colsToPad <- (alignedCols[2] + 2) - ncol(region)

      regionCroppedShifted <- cbind(regionCroppedShifted,
                                    matrix(NA,nrow = nrow(regionCroppedShifted),ncol = colsToPad))

    }
    # now return correctly-dimensioned matrix after padding
    if(all(dim(regionCroppedShifted) == dim(cell))){

      return(list("regionCropped" = regionCroppedShifted,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

    }

    #two copies if rows are off
    if(nrow(regionCroppedShifted) > nrow(cell) & ncol(regionCroppedShifted) == ncol(cell)){
      rowsToCrop <- abs(nrow(regionCroppedShifted) - nrow(cell))

      regionCroppedRowPre <- regionCroppedShifted[-rowsToCrop,]

      # if we shave off some rows from the top of
      # the aligned target cell, this corresponds to *adding* to the top row
      # index relative to the entire region
      alignedRows <- alignedRows + c(rowsToCrop,0)

      return(list("regionCropped" = regionCroppedRowPre,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

    }


    #(Note for future: I commented-out these < if conditional statements because
    #I wasn't able to replicate a situation in which the regionCroppedShifted
    #rows/cols were/ less than those of the cells after the series of padding if
    #statements above). I'm not sure if the padding completely removes the need
    #to check the < conditions for all possible comparisons in perpetuity, so
    #I'm leaving the code here. However, it drags down the codecov "score" if
    #there isn't a test that hits it.

    #two copies if rows are off
    if(nrow(regionCroppedShifted) < nrow(cell) & ncol(regionCroppedShifted) == ncol(cell)){
      rowsToPad <- abs(nrow(regionCroppedShifted) - nrow(cell))

      regionCroppedRowPre <- rbind(matrix(NA,
                                          nrow = rowsToPad,
                                          ncol = ncol(regionCroppedShifted)),
                                   regionCroppedShifted)

      # need to change the alignedRows object if additional rows are padded
      # above. however, it doesn't seem like this if statement is ever executed
      # anyways
      #alignedRows <- ...something...

      return(list("regionCropped" = regionCroppedRowPre,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

    }

    #2 copies if cols are off
    if(ncol(regionCroppedShifted) > ncol(cell) & nrow(regionCroppedShifted) == nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedShifted) - ncol(cell))

      regionCroppedColPre <- regionCroppedShifted[,-colsToCrop]

      # if we shave off some cols from the left of
      # the aligned target cell, this corresponds to *adding* to the left col
      # index relative to the entire region
      alignedCols <- alignedCols + c(colsToCrop,0)

      return(list("regionCropped" = regionCroppedColPre,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

    }

    #2 copies if cols are off
    if(ncol(regionCroppedShifted) < ncol(cell) & nrow(regionCroppedShifted) == nrow(cell)){
      colsToPad <- abs(ncol(regionCroppedShifted) - ncol(cell))

      regionCroppedColPre <- cbind(matrix(NA,
                                          nrow = nrow(regionCroppedShifted),
                                          ncol = colsToPad),
                                   regionCroppedShifted)

      # need to change the alignedRows object if additional rows are padded
      # above. however, it doesn't seem like this if statement is ever executed
      # anyways
      #alignedCols <- ...something...

      return(list("regionCropped" = regionCroppedColPre,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))
    }

    #4 different copies if both dimensions are too large
    if(ncol(regionCroppedShifted) > ncol(cell) & nrow(regionCroppedShifted) > nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedShifted) - ncol(cell))

      rowsToCrop <- abs(nrow(regionCroppedShifted) - nrow(cell))

      regionCroppedBothPre <- regionCroppedShifted[-rowsToCrop,
                                                   -colsToCrop]

      # same logic to changing alignedRows and alignedCols as above, just
      # applied to both now
      alignedRows <- alignedRows + c(rowsToCrop,0)
      alignedCols <- alignedCols + c(colsToCrop,0)


      return(list("regionCropped" = regionCroppedBothPre,
                  "centerRow" = alignedCenter[1] + 2,
                  "centerCol" = alignedCenter[2] + 2,
                  "regionRows" = c((alignedRows[1] + 2),(alignedRows[2] + 2)),
                  "regionCols" = c((alignedCols[1] + 2),(alignedCols[2] + 2))))

    }
  }
}

#'Extract a matrix from the target region of the same dimension as the reference
#'cell depending on the estimated translation calculated from comparison_fft_ccf
#'
#'@name comparison_alignedTargetCell
#'@param cellHeightValues list/tibble column of x3p objects containing a
#'  reference scan's cells (as returned by comparison_cellDivision)
#'@param regionHeightValues list/tibble column of x3p objects containing a
#'  target scan's regions (as returned by comparison_getTargetRegions)
#'@param target the scan to which each cell in the partitioned scan was
#'  compared.
#'@param theta the theta (rotation) value associated with each cellHeightValues,
#'  regionHeightValues pairing
#'@param fft_ccf_df data frame/tibble column containing the data frame of (x,y)
#'  and CCF values returned by comparison_fft_ccf
#'
#'@return a list of x3p objects containing surface matrices extracted from
#'  regionHeightValues of the same dimension as the x3p objects in
#'  cellHeightValues
#'
# #'@examples
#'
#'@export

comparison_alignedTargetCell <- function(cellHeightValues,
                                         regionHeightValues,
                                         target,
                                         theta,
                                         fft_ccf_df){

  targetCells <- purrr::pmap(.l = list(cellHeightValues,
                                      regionHeightValues,
                                      fft_ccf_df),
                            function(cell,region,translations){
                              targetCell <- extractTargetCell(cell = cell$surface.matrix,
                                                              region = region$surface.matrix,
                                                              dx = translations$x,
                                                              dy = translations$y,
                                                              m1 = cell$cmcR.info$centerByVal,
                                                              sd1 = cell$cmcR.info$scaleByVal,
                                                              m2 = region$cmcR.info$centerByVal,
                                                              sd2 = region$cmcR.info$scaleByVal)

                              # if all we care about is extracting the cell from
                              # the target region, not from the whole scan:

                              # ret <- region
                              # ret$surface.matrix <- targetCell$regionCropped
                              # ret$header.info$sizeX <- nrow(targetCell$regionCropped)
                              # ret$header.info$sizeY <- ncol(targetCell$regionCropped)
                              #
                              # ret$cmcR.info$centerRow <- targetCell$centerRow
                              # ret$cmcR.info$centerCol <- targetCell$centerCol
                              #
                              # ret$cmcR.info$regionRows <- targetCell$regionRows
                              # ret$cmcR.info$regionCols <- targetCell$regionCols
                              #
                              # return(ret)

                              # the information in region$cmcR.info and
                              # targetCell can be put together to determine
                              # exactl which rows/cols of the theta-rotated
                              # target the aligned target cell occupies
                              targetScanRows <- region$cmcR.info$regionIndices[c(3)] + targetCell$regionRows - 1
                              targetScanCols <- region$cmcR.info$regionIndices[c(1)] + targetCell$regionCols - 1

                              # standardize the target scan according to the
                              # region to ensure that the extracted cell
                              # contains the same values as the aligned region
                              targetStandard <- (target$surface.matrix - region$cmcR.info$centerByVal)/region$cmcR.info$scaleByVal
                              # rotate the target by the theta used when
                              # aligning
                              rotatedMask <- rotateSurfaceMatrix(targetStandard,theta)

                              rowPad <- 0
                              colPad <- 0

                              # the target scan rows/cols may be less than the 1
                              # or greater than the dimension of the target. We
                              # will need to pad the rotated mask to get the
                              # correct cell
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

                              # create a new x3p object to contain the aligned
                              # target cell
                              ret <- region
                              ret$surface.matrix <- matrix(rotatedMask[targetScanRows[1]:targetScanRows[2],targetScanCols[1]:targetScanCols[2]],
                                                           nrow = nrow(targetCell$regionCropped),
                                                           ncol = ncol(targetCell$regionCropped))
                              ret$header.info$sizeX <- nrow(targetCell$regionCropped)
                              ret$header.info$sizeY <- ncol(targetCell$regionCropped)

                              ret$cmcR.info$centerRow <- targetCell$centerRow
                              ret$cmcR.info$centerCol <- targetCell$centerCol

                              ret$cmcR.info$regionRows <- targetCell$regionRows
                              ret$cmcR.info$regionCols <- targetCell$regionCols



                              return(ret)
                            })

  return(targetCells)
}


#'Calculates correlation between a cell and a matrix of the same dimensions
#'extracted from the cell's associated region.
#'
#'@name comparison_cor
#'@param cellHeightValues list/tibble column of x3p objects containing a
#'  reference scan's cells (as returned by comparison_cellDivision)
#'@param regionHeightValues list/tibble column of x3p objects containing a
#'  target scan's regions (as returned by comparison_getTargetRegions)
#'@param fft_ccf_df data frame/tibble column containing the data frame of (x,y)
#'  and CCF values returned by comparison_fft_ccf
#'@param use argument for stats::cor
#'@return A vector of the same length as the input containing correlation values
#'  at the estimated alignment between each reference cell and its associated
#'  target region
#'@examples
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- fadul1.1_processed %>%
#' comparison_cellDivision(numCells = 64) %>%
#' dplyr::mutate(regionHeightValues =
#'              comparison_getTargetRegions(cellHeightValues = cellHeightValues,
#'                                          target = fadul1.2_processed)) %>%
#' dplyr::mutate(cellPropMissing =
#'            comparison_calcPropMissing(heightValues = cellHeightValues),
#'               regionPropMissing =
#'            comparison_calcPropMissing(heightValues = regionHeightValues)) %>%
#' dplyr::filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
#' dplyr::mutate(cellHeightValues =
#'         comparison_standardizeHeights(heightValues = cellHeightValues),
#'               regionHeightValues =
#'         comparison_standardizeHeights(heightValues = regionHeightValues)) %>%
#' dplyr::mutate(cellHeightValues =
#'                   comparison_replaceMissing(heightValues = cellHeightValues),
#'               regionHeightValues =
#'             comparison_replaceMissing(heightValues = regionHeightValues)) %>%
#' dplyr::mutate(fft_ccf_df = comparison_fft_ccf(cellHeightValues,
#'                                               regionHeightValues)) %>%
#' dplyr::mutate(pairwiseCompCor = comparison_cor(cellHeightValues,
#'                                                regionHeightValues,
#'                                                fft_ccf_df))
#'
#'head(cellTibble)
#'
#'@export

comparison_cor <- function(cellHeightValues,
                           regionHeightValues,
                           fft_ccf_df,
                           use = "pairwise.complete.obs"){

  rawCors <- purrr::pmap_dbl(.l = list(cellHeightValues,
                                       regionHeightValues,
                                       fft_ccf_df),
                             function(cell,region,translations){
                               rawCor <- calcRawCorr(cell = cell$surface.matrix,
                                                     region = region$surface.matrix,
                                                     dx = translations$x,
                                                     dy = translations$y,
                                                     m1 = cell$cmcR.info$centerByVal,
                                                     sd1 = cell$cmcR.info$scaleByVal,
                                                     m2 = region$cmcR.info$centerByVal,
                                                     sd2 = region$cmcR.info$scaleByVal,
                                                     use = use)

                               return(rawCor)
                             })

  return(rawCors)
}


#'Performs all steps in the cell-based comparison procedure.
#'
#'@name comparison_allTogether
#'
#'@param reference an x3p object containing a breech face scan to be treated as
#'  the "reference scan" partitioned into a grid of cells
#'@param target an x3p object containing a breech face scan to be treated as the
#'  "target scan" that the reference scan's cells are compared to
#'@param theta degrees that the target scan is to be rotated prior extracting
#'  regions.
#'@param numCells number of cells to partition the breech face scan into. Must
#'  be a perfect square (49, 64, 81, etc.)
#'@param maxMissingProp maximum proportion of missing values allowed for each
#'  cell/region.
#'@param regionSizeMultiplier an integer indicating the ratio between the areas
#'  of each target region and reference cell. Must be a perfect square (4, 9,
#'  16, etc.)
#'@param returnX3Ps boolean to return the cellHeightValues and
#'  alignedTargetCells for each cell index. Note that setting this argument to
#'  TRUE significantly increases the size of the returned object.
#'
#'  data(fadul1.1_processed,fadul1.2_processed)
#'
#'  comparisonDF <- comparison_allTogether(reference = fadul1.1_processed,
#'  target = fadul1.2_processed)
#'
#'  head(comparisonDF)
#'
#'@return a tibble object containing cell indices and the x, y, FFT-based CCF,
#'  and pairwise-complete correlation associated with the comparison between
#'  each cell and its associated target scan region (after rotating the target
#'  scan by theta degrees)
#'
#'@examples
#'
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- comparison_allTogether(reference = fadul1.1_processed,target = fadul1.2_processed)
#'
#'head(cellTibble)
#'
#'@importFrom rlang .data
#'@export

comparison_allTogether <- function(reference,
                                   target,
                                   theta = 0,
                                   numCells = 64,
                                   maxMissingProp = .85,
                                   regionSizeMultiplier = 9,
                                   returnX3Ps = FALSE){

  ret <- reference %>%
    comparison_cellDivision(numCells)  %>%
    dplyr::mutate(cellPropMissing = comparison_calcPropMissing(.data$cellHeightValues),
                  refMissingCount = purrr::map_dbl(.data$cellHeightValues,~ sum(is.na(.$surface.matrix)))) %>%
    dplyr::filter(.data$cellPropMissing <= maxMissingProp) %>%
    dplyr::mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = .data$cellHeightValues,
                                                                   target = target,
                                                                   theta = theta,
                                                                   regionSizeMultiplier = regionSizeMultiplier)) %>%
    dplyr::mutate(targMissingProp = comparison_calcPropMissing(.data$regionHeightValues),
           targMissingCount = purrr::map_dbl(.data$regionHeightValues,~ sum(is.na(.$surface.matrix)))) %>%
    dplyr::filter(.data$targMissingProp <= maxMissingProp) %>%
    dplyr::mutate(cellHeightValues = comparison_standardizeHeights(.data$cellHeightValues),
                  regionHeightValues = comparison_standardizeHeights(.data$regionHeightValues)) %>%
    dplyr::mutate(cellHeightValues_replaced = comparison_replaceMissing(.data$cellHeightValues),
                  regionHeightValues_replaced = comparison_replaceMissing(.data$regionHeightValues)) %>%
    dplyr::mutate(fft_ccf_df = comparison_fft_ccf(cellHeightValues = .data$cellHeightValues_replaced,
                                                  regionHeightValues = .data$regionHeightValues_replaced)) %>%
    dplyr::mutate(alignedTargetCell = comparison_alignedTargetCell(cellHeightValues = .data$cellHeightValues,
                                                                   regionHeightValues = .data$regionHeightValues,
                                                                   target = target,
                                                                   theta = theta,
                                                                   fft_ccf_df = .data$fft_ccf_df)) %>%
    dplyr::mutate(jointlyMissing = purrr::map2_dbl(.data$cellHeightValues,.data$alignedTargetCell,~ sum(is.na(.x$surface.matrix) & is.na(.y$surface.matrix))),
                  pairwiseCompCor = purrr::map2_dbl(.data$cellHeightValues,.data$alignedTargetCell,
                                             ~ cor(c(.x$surface.matrix),c(.y$surface.matrix),
                                                   use = "pairwise.complete.obs"))) %>%
    tidyr::unnest(.data$fft_ccf_df) %>%
    dplyr::mutate(theta = theta)

  if(!returnX3Ps){

    ret <- ret %>%
      dplyr::select(.data$cellIndex,.data$x,.data$y,.data$fft_ccf,.data$pairwiseCompCor,.data$theta,.data$refMissingCount,.data$targMissingCount,.data$jointlyMissing)

  }
  else{

    ret <- ret %>%
      dplyr::select(.data$cellIndex,.data$x,.data$y,.data$fft_ccf,.data$pairwiseCompCor,.data$theta,.data$refMissingCount,.data$targMissingCount,.data$jointlyMissing,.data$cellHeightValues,.data$alignedTargetCell)

  }

  return(ret)

}