#' Divides an image into cells
#' @name cellDivision
#'
#' @param surfaceMat surface matrix of a breech face impression
#' @param cellNumHoriz number of splits along horizontal axis
#' @param cellNumVert number of splits along vertical axis
#'
#' @description This is a helper for the cellCCF function.
#'
#' @keywords internal

cellDivision <- function(surfaceMat,
                         cellNumHoriz = 8,
                         cellNumVert = 8){

  splitSurfaceMat <- surfaceMat %>%
    imager::as.cimg() %>%
    imager::imsplit(axis = "x",
                    nb = cellNumHoriz) %>%
    purrr::map(.f = ~ imager::imsplit(.x,
                                      axis = "y",
                                      nb = cellNumVert)) %>%
    purrr::map_depth(.depth = 2,
                     .f = as.matrix)

  return(splitSurfaceMat)
}

#' Cuts out a cell in a surface matrix.
#' @name extractCellbyCornerLocs
#'
#' @param cornerLocs the location (indices) of the desired cell in the image
#' @param rotatedSurfaceMat a surface matrix representing a rotated breech face
#'   impression image
#'
#' @description This is a helper for the cellCCF function. This function
#'   cuts out a cell in a surface matrix based on the most top, bottom, left,
#'   and right indices of the cell's location in the surface matrix.
#'
#' @keywords internal

extractCellbyCornerLocs <- function(cornerLocs,
                                    rotatedSurfaceMat,
                                    mat2Dim){
  #perform the appropriate subsetting of image A to create a list of larger
  #cells than those in image B
  splitRotatedSurfaceMat <- rotatedSurfaceMat[cornerLocs[["top"]]:cornerLocs[["bottom"]],
                                              cornerLocs[["left"]]:cornerLocs[["right"]]]

  if(nrow(splitRotatedSurfaceMat) != ncol(splitRotatedSurfaceMat)){ #if the matrix isn't square...

    #if the rows need padding...
    if(nrow(splitRotatedSurfaceMat) < max(dim(splitRotatedSurfaceMat))){

      rowsToPad <- ncol(splitRotatedSurfaceMat) - nrow(splitRotatedSurfaceMat)
      rowPadder <- matrix(NA,nrow = rowsToPad,ncol = ncol(splitRotatedSurfaceMat))

      #if the split comes from the top of the overall matrix...
      if(cornerLocs[["top"]] == 1){
        splitRotatedSurfaceMat <- rbind(rowPadder,
                                        splitRotatedSurfaceMat)
      }

      #if the split comes from the bottom of the overall matrix....
      if(cornerLocs[["bottom"]] == mat2Dim[1]){
        splitRotatedSurfaceMat <- rbind(splitRotatedSurfaceMat,
                                        rowPadder)
      }
    }

    #if the cols need padding...
    if(ncol(splitRotatedSurfaceMat) < max(dim(splitRotatedSurfaceMat))){

      colsToPad <- nrow(splitRotatedSurfaceMat) - ncol(splitRotatedSurfaceMat)
      colPadder <- matrix(NA,ncol = colsToPad,nrow = nrow(splitRotatedSurfaceMat))

      #if the split comes from the left side of the overall matrix...
      if(cornerLocs[["left"]] == 1){
        splitRotatedSurfaceMat <- cbind(colPadder,
                                        splitRotatedSurfaceMat)
      }
      #if the split comes from the right side of the overall matrix....
      if(cornerLocs[["right"]] == mat2Dim[2]){
        splitRotatedSurfaceMat <- cbind(splitRotatedSurfaceMat,
                                        colPadder)
      }
    }
  }

  return(splitRotatedSurfaceMat)
}

#' @name rotateSurfaceMatrix
#'
#' @keywords internal

utils::globalVariables(c("."))

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

#' Determines whether a matrix contains enough non-NA values.
#' @name checkForBreechface
#'
#' @param cell a matrix representing a cell from a breech face impression image
#' @description This is a helper for the cellCCF function. This function
#'   returns TRUE if a given cell contains more than 225 pixels of observed
#'   values (i.e., non-NA values) and FALSE otherwise.
#'
#' @keywords internal

checkForBreechface <- function(cell,
                               minObservedProp = .1){
  containsBreechfaceBool <- (sum(!is.na(cell)) > minObservedProp*length(as.vector(cell)))
  return(containsBreechfaceBool)
}

#' Subtract the average pixel value from an image and replace NAs with 0
#' @name standardizeSurfaceMat
#'
#' @param surfaceMat surface matrix
#' @param mean average pixel value to shift an image by
#'
#' @description This is a helper for the cellCCF function.
#'
#' @keywords internal

standardizeSurfaceMat <- function(surfaceMat,
                                  m,
                                  s){
  surfaceMat <- (surfaceMat - m)/s
  surfaceMat[is.na(surfaceMat)] <- 0
  return(surfaceMat)
}

#' @name splitSurfaceMat1
#'
#' @keywords internal

splitSurfaceMat1 <- function(surfaceMat,cellNumHoriz,cellNumVert,minObservedProp){

  surfaceMat_split <- cellDivision(surfaceMat,
                                   cellNumHoriz = cellNumHoriz,
                                   cellNumVert = cellNumVert) #split image 1 into cells

  #create a cell ID column based on x,y location in image:
  cellRanges <- purrr::map(names(surfaceMat_split),
                           function(horizCell){
                             purrr::map(.x = names(surfaceMat_split[[1]]),
                                        function(vertCell) paste(horizCell,vertCell,sep = ","))
                           }) %>%
    unlist() #unlist to make it a character vector to be added to a df

  cellSideLengths <- surfaceMat_split %>%
    purrr::map_depth(.depth = 2,
                     ~ c("row" = nrow(.),
                         "col" = ncol(.))) %>%
    purrr::flatten()

  #Determine whether a cell contains more than 10 pixels of breechface
  #information (for some reason checking for more than 0 pixels of breechface
  #was causing the function to fail)
  matPixCounts <- surfaceMat_split %>%
    purrr::map_depth(.depth = 2,
                     ~ checkForBreechface(cell = .,minObservedProp = minObservedProp)) %>%
    unlist()

  return(list(surfaceMat_split = surfaceMat_split,
              cellRanges = cellRanges,
              cellSideLengths = cellSideLengths,
              mat1PixCounts = matPixCounts))
}

#' @name getMat2SplitIndices
#'
#' @keywords internal
#'
#' @importFrom stats setNames

utils::globalVariables(c("."))

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
                  expandedCellCorners[expandedCellCorners <= 0] <- 1
                  #replace indices greater than the maximum index with the
                  #maximum index (right/bottom-most cells): Note that imager
                  #treats the rows of a matrix as the "x" axis and the columns
                  #as the "y" axis, contrary to intuition - effectively treating
                  #a matrix as its transpose. As such, we need
                  #to swap the dimensions for when we subset the image further
                  #down in the function
                  if(expandedCellCorners[c("right")] > mat2Dim[2]){
                    expandedCellCorners[c("right")] <- mat2Dim[2]
                  }
                  if(expandedCellCorners[c("bottom")] > mat2Dim[1]){
                    expandedCellCorners[c("bottom")] <- mat2Dim[1]
                  }

                  return(expandedCellCorners)
                }) %>%
    setNames(cellRanges)

  return(mat2_splitCorners)
}

#' @name swapcellRangeAxes
#'
#' @keywords internal

swapcellRangeAxes <- function(cellRange){
  sSplit <- stringr::str_split(string = cellRange,pattern = ",",n = 2)

  paste0(stringr::str_replace(string = sSplit[[1]][1],pattern = "x =",replacement = "rows:"),
         ", ",
         stringr::str_replace(string = sSplit[[1]][2],pattern = "y =",replacement = "cols:"))
}

#' Calculates the cross-correlation function between cells in one matrix to
#' regions in another.
#'
#' @name cellCCF
#'
#' @export
#'
#' @param x3p1 an x3p object containing the surface matrix of a cartridge case
#'   scan
#' @param x3p2 an x3p object containing the surface matrix of a cartridge case
#'   scan to be compared to that in x3p1
#' @param thetas rotation values (in degrees) for which x3p2$surface.matrix will
#'   be rotated, split into cells, and compared to x3p1$surface.matrix
#' @param cellNumHoriz number of cells along horizontal axis to divide
#'   x3p1$surface.matrix into
#' @param cellNumVert number of cells along vertical axis to divide
#'   x3p1$surface.matrix into
#' @param regionToCellProp determines how much larger the x3p2 regions will be
#'   relative to the x3p1 cells. For example, if regionToCellProp = 4 means that
#'   the x3p2 regions will be 4 times times larger (sidelengths multiplied by
#'   2).
#' @param minObservedProp the minimum proportion of a cell that needs to contain
#'   observed (i.e., non-NA) values for it to be included in the CCF calculation
#'   procedure.
#' @param centerCell  dictates if cell is to be shifted prior to CCF
#'   calculation. Default is that no shifting is performed. If set to
#'   "wholeMatrix", then each cell is subracted by the mean of the whole matrix.
#'   If set to "individualCell", then each cell is subtracted by its particular
#'   mean.
#' @param scaleCell dictates if cell is to be scaled prior to CCF calculation.
#'   Default is that no scaling is performed. If set to "wholeMatrix", then each
#'   cell is divided by the standard deviation of the whole matrix. If set to
#'   "individualCell", then each cell is divided by its particular standard
#'   deviation.
#' @param ccfMethod dictates which of 3 different methods are used to calculate
#'   a final maximum cross-correlation estimate. "fft" estimates the
#'   optimal translation values to align each cell within in its associated
#'   region by using Cross-Correlation theorem + FFT algorithm and then
#'   calculates the pairwise-complete correlation between the cell and a
#'   cell-sized matrix extracted from the paired region. "imager" uses the
#'   normalized imager::correlate function and uses the max CCF value
#'   calculated. "bruteForceReweighted" calculates the pairwise-complete
#'   correlation between a cell and every possible cell-sized matrix that could
#'   be extracted from its paired region (and is thus computationally very
#'   costly). For each cell/region pair, these correlation values are
#'   re-weighted based on the number of non-missing values used to calculate
#'   them by the following: nonMissingCount*cor/max(nonMissingCount).
#' @param rawCorrTieBreaker Only applicable if ccfMethod == "fft".
#'   The way in which the CCF (see ccfMethod) is calculated may require slight
#'   padding/cropping of the mat1-sized matrix extracted from mat2 to make their
#'   dimensions equal (e.g., "center" of mat1-sized matrix in mat2 may be one of
#'   4 pixels). This padding/cropping can occur to the initial or final
#'   rows/cols in the matrix, without a clear way to determine which is
#'   "correct." As such, all possible combinations of pre/post padding/cropping
#'   are considered (only if necessary). To determine a final mat1-sized matrix,
#'   rawCorrTieBreaker can be used to determine which yields the lowest/highest
#'   correlation with mat1 (using rawCorrTieBreaker = which.min or which.max,
#'   respectively).
#' @param use Only applicable if ccfMethod == "fft". argument to be
#'   passed to the cor function. Dictates how NAs are dealt with in computing
#'   the correlation.
#'
#' @return The list allResults contains the CCF values, horizontal, and vertical
#'   translations calculated for each cell, for each rotation value. The data
#'   frame topResults contains the CCF value, associated horizontal/vertical
#'   translation, and rotation value at which each cell in mat1 achieved its
#'   highest CCF with its paired cell in mat2.
#'
#' @description This function performs a cell-based comparison of two cartridge
#'   case scans as described in  'Proposed "Congruent Matching Cells (CMC)
#'   Method for Ballistic Identification and Error Rate Estimation' by John Song
#'   (2015). The method works by first dividing mat1, the surface matrix
#'   representing the height values of a breech face impression microscopy scan,
#'   into pairwise disjoint cells. mat2, the surface matrix of a breech face
#'   impression scan to be compared to mat1, is also broken up into cells
#'   centered at the same location as those in mat1. However, the cells in mat2
#'   are larger than those in mat1 (typically 9 times the size except on the
#'   border of mat2). See ccfMethod argument details for information regarding
#'   how the CCF is calculated.
#'
#'   mat2 is then rotated by some amount (say, 3 degrees) and the process of
#'   dividing mat2 into cells and calculating the max CCF for each cell in mat1
#'   is repeated. The rotations performed on mat2 are dictated by the value(s)
#'   passed to the theta argument.
#'
#' @seealso
#' \url{https://pdfs.semanticscholar.org/4bf3/0b3a23c38d8396fa5e0d116cba63a3681494.pdf}
#' @seealso cmcR::cmcFilter
#'
#' @importFrom stats sd setNames
#'

utils::globalVariables(c("cellRange","."))

cellCCF <- function(x3p1,
                    x3p2,
                    thetas = seq(-30,30,by = 3),
                    cellNumHoriz = 8,
                    cellNumVert = cellNumHoriz,
                    regionToCellProp = 4,
                    minObservedProp = .1,
                    centerCell = "individualCell",
                    scaleCell = "individualCell",
                    ccfMethod = "fft",
                    rawCorrTieBreaker = which.max,
                    use = "pairwise.complete.obs"){
  #Needed tests: mat1 and mat2 must be matrices thetas, cellNumHoriz, and
  #cellNumVert should be integers (cellNumHoriz and cellNumVert would optimally be
  #equal - maybe print a warning if not?)

  mat1 <- x3p1$surface.matrix
  mat2 <- x3p2$surface.matrix

  if(missing(centerCell)){
    m1 <- 0
    m2 <- 0
  }
  if(missing(scaleCell)){
    sd1 <- 1
    sd2 <- 1
  }

  if(!missing(centerCell)){
    if(centerCell == "wholeMatrix"){
      m1 <- mean(as.vector(mat1),na.rm = TRUE)

      m2 <- mean(as.vector(mat2),na.rm = TRUE)
    }
  }

  if(!missing(scaleCell)){
    if(scaleCell == "wholeMatrix"){
      sd1 <- sd(as.vector(mat1),na.rm = TRUE)

      sd2 <- sd(as.vector(mat2),na.rm = TRUE)
    }
  }

  #initialize list to contain all CCF values for each theta value
  allResults <- purrr::map(thetas,
                           function(theta) theta = data.frame(cell_ID = NA,
                                                              ccf = NA,
                                                              dx = NA,
                                                              dy = NA)) %>%
    setNames(thetas)

  mat1_split <- splitSurfaceMat1(surfaceMat = mat1,
                                 cellNumHoriz = cellNumHoriz,
                                 cellNumVert = cellNumVert,
                                 minObservedProp = minObservedProp)

  #creating this df is necessary so that we can compare cells in similar
  #positions between two comparisons (since their "cellRange" may differ since
  #the scans aren't the same size)
  cellRangedf <- data.frame(cellNum = seq_along(mat1_split$cellRanges),
                            cellRange = mat1_split$cellRanges)

  #Now we want to split image B into cells with the same centers as those in
  #image A, but with twice the side length (these wider cells will intersect
  #each other). We will first get the dimensions (the x,y locations of each
  #cells' corners) of where each cell should be in image B.
  sidelengthMultiplier <- floor(sqrt(regionToCellProp))

  mat2_splitCorners <- getMat2SplitIndices(cellRanges = mat1_split$cellRanges,
                                           cellSideLengths = mat1_split$cellSideLengths,
                                           mat2Dim = dim(mat2),
                                           sidelengthMultiplier = sidelengthMultiplier)

  for(theta in thetas){
    #rotate image 2 and split into cells:

    mat2_rotated <- mat2 %>%
      rotateSurfaceMatrix(theta)

    #Now that we've rotated image B, we want to create a list consisting of the
    #cells that we are to compare to image A's cells. These will be larger than
    #Image A's cells (up to image B's boundary) and defined by the
    #expandedCellCorners x,y pairs calculated above. mat2_splitRotated contains
    #regions that lie on the edge of mat2. The getMat2SplitLoctaions function
    #ensures that the way in which we define these regions doesn't extend past
    #the dimensions of mat2. However, the corresponding dx and dy values
    #assocated with the max CCF, while "correct" for a particular matrix pair in
    #the sense that they indeed yield the correct trnaslation values to align
    #the mat2 region to the mat1 cell, are difficult to compare between the
    #nicely square interior regions (the dx and dy values often disagree
    #considerably, which is a problem when we look for consensus among the dx
    #and dy values). Thus, it is necessary to pad these edge images until they
    #are square. Through experimentation, it has been observed that padding the
    #matrices with constant values does not affect the CCF score. We determine
    #which matrices need to be padded based on whether the associated cellRange
    #contains 1 or the maximum value of the cellRanges (containing either part of
    #the first row/col of the matrix or the last)
    mat2_splitRotated <-
      purrr::map(.x = mat2_splitCorners,
                 ~ extractCellbyCornerLocs(cornerLocs = .x,
                                           rotatedSurfaceMat = mat2_rotated,
                                           mat2Dim = dim(mat2)))

    #Determine whether a cell contains more than 10 pixels of breechface
    #information (for some reason checking for more than 0 pixels of
    #breechface was causing the function to fail)
    mat2PixCounts <- mat2_splitRotated %>%
      purrr::map(~ checkForBreechface(cell = .,minObservedProp = minObservedProp)) %>%
      unlist()

    #remove cells that don't include enough breechface:
    mat1_splitFiltered <- purrr::flatten(mat1_split$surfaceMat_split)[mat1_split$mat1PixCounts & mat2PixCounts]
    mat2_splitFiltered <- mat2_splitRotated[mat1_split$mat1PixCounts & mat2PixCounts]

    #grab the cell IDs for each cell not removed above. This is used to update
    #topResults below
    filteredcellRange <- mat1_split$cellRanges[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #creating this df is necessary so that we can compare cells in similar
    #positions between two comparisons (since their "cellRange" may differ since
    #the scans aren't the same size, but their absolute position in the grid
    #(e.g., 4th from the left in the top row) will remain the same)
    cellRangedf_filtered <- cellRangedf %>%
      dplyr::filter(cellRange %in% filteredcellRange) %>%
      #imager swaps x and y axes, so we need to swap them back to be more
      #interpretable:
      dplyr::mutate(cellRange = purrr::map_chr(cellRange,swapcellRangeAxes))

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs)
    if(!missing(centerCell)){
      if(centerCell == "individualCell"){
        m1 <- mat1_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredcellRange)

        m2 <-  mat2_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredcellRange)
      }
    }

    if(!missing(scaleCell)){
      if(scaleCell == "individualCell"){
        sd1 <-  mat1_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredcellRange)

        sd2 <-  mat2_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredcellRange)
      }
    }

    stopifnot(all(sd1 > 0),all(sd2 > 0))

    mat1_splitShifted <- purrr::pmap(list(mat1_splitFiltered,m1,sd1),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3)) %>%
      setNames(filteredcellRange) #Need to set the names to avoid repeated list
    #labels

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs).
    mat2_splitShifted <- purrr::pmap(list(mat2_splitFiltered,m2,sd2),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3)) %>%
      setNames(filteredcellRange) #Need to set the names to avoid repeated list
    #labels

    #calculate the correlation of each cell pair
    ccfValues <- purrr::map2_dfr(.x = mat1_splitShifted,
                                 .y = mat2_splitShifted,
                                 .f = ~ data.frame(purrr::flatten(ccfComparison(.x,.y,ccfMethod = ccfMethod)))) %>% #returns a nested list of ccf,dx,dy values
      dplyr::bind_cols(cellRangedf_filtered,.) %>%
      dplyr::mutate(theta = rep(theta,times = nrow(.)),
                    nonMissingProportion = purrr::map_dbl(mat1_splitFiltered,
                                                          .f = ~ (sum(!is.na(.)))/prod(dim(.))))

    if(ccfMethod == "fft"){
      ccfValues <- ccfValues %>%
        dplyr::mutate(ccf = purrr::pmap_dbl(.l = list(mat1_splitFiltered,
                                                      mat2_splitFiltered,
                                                      .$x,
                                                      .$y,
                                                      m1,
                                                      m2,
                                                      sd1,
                                                      sd2),
                                            .f = ~ calcRawCorr(cell = ..1,
                                                               region = ..2,
                                                               dx = ..3,
                                                               dy = ..4,
                                                               m1 = ..5,
                                                               m2 = ..6,
                                                               sd1 = ..7,
                                                               sd2 = ..8,
                                                               tieBreaker = rawCorrTieBreaker,
                                                               use = use)))
    }

    allResults[paste0(theta)][[1]] <- ccfValues
  }

  #rearrange columns in allResults for readability
  allResults <- allResults %>%
    purrr::map(~ dplyr::select(.,cellNum,cellRange,dplyr::everything()))

  if(missing(centerCell)){
    centerCell <- "none"
  }
  if(missing(scaleCell)){
    scaleCell <- "none"
  }

  return(list(
    "params" = list("theta" = thetas,
                    "cellNumHoriz" = cellNumHoriz,
                    "cellNumVert" = cellNumVert,
                    "regionToCellProp" = regionToCellProp,
                    "minObservedProp" = minObservedProp,
                    "centerCell" = centerCell,
                    "mat1Shift" = m1,
                    "mat2Shift" = m2,
                    "scaleCell" = scaleCell,
                    "rawCorrTieBreaker" = rawCorrTieBreaker,
                    "use" = use,
                    "mat1ScaleFactor" = sd1,
                    "mat2ScaleFactor" = sd2,
                    "ccfMethod" = ccfMethod),
    "ccfResults" = allResults
  ))
}

#' Applies cellCCF function twice for a pair of cartridge case scans.
#'
#' @name cellCCF_bothDirections
#'
#' @description Wrapper for applying the cmcR::cellCCF function to x3p1 vs. x3p2
#'   and again for x3p2 vs. x3p1. See cellCCF function documentation for more
#'   details.
#'
#' @param x3p1 an x3p object containing the surface matrix of a cartridge case
#'   scan
#' @param x3p2 an x3p object containing the surface matrix of a cartridge case
#'   scan to be compared to that in x3p1
#' @param thetas rotation values (in degrees) for which x3p2$surface.matrix will
#'   be rotated, split into cells, and compared to x3p1$surface.matrix
#' @param cellNumHoriz number of cells along horizontal axis to divide
#'   x3p1$surface.matrix into
#' @param cellNumVert number of cells along vertical axis to divide
#'   x3p1$surface.matrix into
#' @param regionToCellProp determines how much larger the x3p2 regions will be
#'   relative to the x3p1 cells. For example, if regionToCellProp = 4 means that
#'   the x3p2 regions will be 4 times times larger (sidelengths multiplied by
#'   2).
#' @param minObservedProp the minimum proportion of a cell that needs to contain
#'   observed (i.e., non-NA) values for it to be included in the CCF calculation
#'   procedure.
#' @param centerCell  dictates if cell is to be shifted prior to CCF
#'   calculation. Default is that no shifting is performed. If set to
#'   "wholeMatrix", then each cell is subracted by the mean of the whole matrix.
#'   If set to "individualCell", then each cell is subtracted by its particular
#'   mean.
#' @param scaleCell dictates if cell is to be scaled prior to CCF calculation.
#'   Default is that no scaling is performed. If set to "wholeMatrix", then each
#'   cell is divided by the standard deviation of the whole matrix. If set to
#'   "individualCell", then each cell is divided by its particular standard
#'   deviation.
#' @param ccfMethod dictates which of 3 different methods are used to calculate
#'   a final maximum cross-correlation estimate. "fft" estimates the
#'   optimal translation values to align each cell within in its associated
#'   region by using Cross-Correlation theorem + FFT algorithm and then
#'   calculates the pairwise-complete correlation between the cell and a
#'   cell-sized matrix extracted from the paired region. "imager" uses the
#'   normalized imager::correlate function and uses the max CCF value
#'   calculated. "bruteForceReweighted" calculates the pairwise-complete
#'   correlation between a cell and every possible cell-sized matrix that could
#'   be extracted from its paired region (and is thus computationally very
#'   costly). For each cell/region pair, these correlation values are
#'   re-weighted based on the number of non-missing values used to calculate
#'   them by the following: nonMissingCount*cor/max(nonMissingCount).
#' @param rawCorrTieBreaker Only applicable if ccfMethod == "fft".
#'   The way in which the CCF (see ccfMethod) is calculated may require slight
#'   padding/cropping of the mat1-sized matrix extracted from mat2 to make their
#'   dimensions equal (e.g., "center" of mat1-sized matrix in mat2 may be one of
#'   4 pixels). This padding/cropping can occur to the initial or final
#'   rows/cols in the matrix, without a clear way to determine which is
#'   "correct." As such, all possible combinations of pre/post padding/cropping
#'   are considered (only if necessary). To determine a final mat1-sized matrix,
#'   rawCorrTieBreaker can be used to determine which yields the lowest/highest
#'   correlation with mat1 (using rawCorrTieBreaker = which.min or which.max,
#'   respectively).
#' @param use Only applicable if ccfMethod == "fft". argument to be
#'   passed to the cor function. Dictates how NAs are dealt with in computing
#'   the correlation.
#'
#' @examples
#' \dontrun{
#'  #x3p1 and x3p2 are assumed to be 2 processed scans
#' comparison1 <- cellCCF_bothDirections(x3p1 = x3p1,
#'                                       x3p2 = x3p2,
#'                                       thetas = seq(-30,30,by = 30),
#'                                       cellNumHoriz = 8,
#'                                       regionToCellProp = 4,
#'                                       centerCell = "individualCell",
#'                                       scaleCell = "individualCell")
#'
#' comparison1$comparison_1to2 #comparison of x3p1 vs x3p2
#' comparison1$comparison_2to1 #comparison of x3p2 vs x3p1
#' }
#'
#' @seealso cmcR::cellCCF_improved
#' @export

cellCCF_bothDirections <- function(x3p1,
                                   x3p2,
                                   thetas = seq(-30,30,by = 3),
                                   cellNumHoriz = 8,
                                   cellNumVert = cellNumHoriz,
                                   regionToCellProp = 4,
                                   minObservedProp = .1,
                                   centerCell = "individualCell",
                                   scaleCell = "individualCell",
                                   ccfMethod = "fft",
                                   rawCorrTieBreaker = which.max,
                                   use = "pairwise.complete.obs"){



  comparison_1to2 <- cellCCF(x3p1 = x3p1,
                             x3p2 = x3p2,
                             thetas = thetas,
                             cellNumHoriz = cellNumHoriz,
                             cellNumVert = cellNumVert,
                             regionToCellProp = regionToCellProp,
                             minObservedProp = minObservedProp,
                             rawCorrTieBreaker = rawCorrTieBreaker,
                             use = use,
                             centerCell = centerCell,
                             scaleCell = scaleCell,
                             ccfMethod = ccfMethod)

  comparison_2to1 <- cellCCF(x3p1 = x3p2,
                             x3p2 = x3p1,
                             thetas = thetas,
                             cellNumHoriz = cellNumHoriz,
                             cellNumVert = cellNumVert,
                             regionToCellProp = regionToCellProp,
                             minObservedProp = minObservedProp,
                             rawCorrTieBreaker = rawCorrTieBreaker,
                             use = use,
                             centerCell = centerCell,
                             scaleCell = scaleCell,
                             ccfMethod = ccfMethod)

  return(list("comparison_1to2" = comparison_1to2,
              "comparison_2to1" = comparison_2to1))
}
#' Split a reference scan into a grid of cells
#'
#' @name comparison_cellDivision
#'
#' @export

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
                     .f = ~ set_names(as.matrix(.),NULL))

  cellRanges <- purrr::map(names(splitSurfaceMat),
                           function(horizCell){
                             purrr::map(.x = names(splitSurfaceMat[[1]]),
                                        function(vertCell) cmcR:::swapcellRangeAxes(paste(horizCell,vertCell,sep = ",")))
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

                  #include which rows/column in the original scan each cell was taken from
                  cell_x3p$cmcR.info$cellRange <- cellRange

                  return(cell_x3p)
                })

  cellTibble <- tibble::tibble(cellNum = 1:numCells,
                               cellHeightValues = splitSurfaceMat) %>%
    dplyr::mutate(cellIndex = cmcR:::linear_to_matrix(index = cellNum,
                                                      nrow = ceiling(sqrt(max(cellNum))),
                                                      byrow = TRUE)) %>%
    dplyr::select(cellIndex,cellHeightValues)

  return(cellTibble)
}

#' Calculate the proportion of missing values in a surface matrix
#'
#' @name comparison_calcPropMissing
#'
#' @export

comparison_calcPropMissing <- function(heightValues){
  heightValues %>%
    purrr::map_dbl(~ sum(is.na(.$surface.matrix))/length(.$surface.matrix))
}

#' Extract regions from a target scan based on associated cells in reference
#' scan
#'
#' @name comparison_getTargetRegions
#'
#' @export

comparison_getTargetRegions <- function(cellHeightValues,
                                        target_x3p,
                                        rotation = 0,
                                        regionSizeMultiplier = 9){

  cellSideLengths <- cellHeightValues %>%
    purrr::map(~ c("row" = nrow(.$surface.matrix),
                   "col" = ncol(.$surface.matrix)))

  cellRange <- cellHeightValues %>%
    purrr::map_chr(~ .$cmcR.info$cellRange)

  target_x3p_regionIndices <- cmcR:::getMat2SplitIndices(cellRanges = cellRange,
                                                         cellSideLengths = cellSideLengths,
                                                         mat2Dim = dim(target_x3p$surface.matrix),
                                                         sidelengthMultiplier = floor(sqrt(regionSizeMultiplier)))

  target_surfaceMat_rotated <- cmcR:::rotateSurfaceMatrix(target_x3p$surface.matrix,
                                                          theta = rotation)


  target_x3p_splitRotated <-
    purrr::map(.x = target_x3p_regionIndices,
               function(cornerIndices){
                 regionMatrix <- cmcR:::extractCellbyCornerLocs(cornerLocs = cornerIndices,
                                                                rotatedSurfaceMat = target_surfaceMat_rotated,
                                                                mat2Dim = dim(target_x3p$surface.matrix))

                 region_x3p <- x3ptools::df_to_x3p(data.frame(x = 1,y = 1,value = NA))

                 region_x3p$surface.matrix <- regionMatrix

                 #update metainformation
                 region_x3p$header.info <- target_x3p$header.info
                 region_x3p$header.info$sizeY <- ncol(regionMatrix)
                 region_x3p$header.info$sizeX <- nrow(regionMatrix)

                 return(region_x3p)
               } )
}

#' Standardize height values of a scan by centering/scaling by desired
#' statistics and replacing missing values
#'
#' @name comparison_standardizeHeightValues
#'
#' @note this function adds information to the metainformation of the x3p scan
#'   it is given that is required for calculating, for example, the
#'   pairwise-complete correlation using the comparison_cor function.
#'
#' @export

comparison_standardizeHeightValues <- function(heightValues,
                                               withRespectTo = "individualCell",
                                               centerBy = mean,
                                               scaleBy = sd){

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

#' Replace missing values in a scan
#'
#' @name comparison_replacingMissingValues
#'
#' @export

comparison_replacingMissingValues <- function(heightValues,
                                              replacement = 0){
  replacedHeightValues <- heightValues %>%
    purrr::map(function(x3p){
      x3p$surface.matrix[is.na(x3p$surface.matrix)] <- 0

      return(x3p)
    })

  return(replacedHeightValues)
}

#' Estimate translation alignment between a cell/region pair based on the
#' [Cross-Correlation
#' Theorem](https://mathworld.wolfram.com/Cross-CorrelationTheorem.html)
#'
#' @name comparison_fft.ccf
#'
#' @note The FFT is not defined for matrices containing missing values. The
#'   missing values in the cell and region need to be replaced before using this
#'   function. See the \link[cmcR](comparison_standardizeHeightValues) function
#'   to replace missing values after standardization.
#'
#' @export
comparison_fft.ccf <- function(cellHeightValues,regionHeightValues){
  ccfList <- purrr::map2(cellHeightValues,
                         regionHeightValues,
                         ~ cmcR:::ccfComparison(mat1 = .x$surface.matrix,mat2 = .y$surface.matrix,ccfMethod = "fft"))

  return(ccfList)
}

#' Calculates correlation between a cell and a matrix of the same dimensions
#' extracted from the cell's associated region.
#'
#' @name comparison_cor
#'
#'
#' @export

comparison_cor <- function(cellHeightValues,
                           regionHeightValues,
                           fft.ccf_df,
                           use = "pairwise.complete.obs"){

  rawCors <- purrr::pmap_dbl(.l = list(cellHeightValues,
                                       regionHeightValues,
                                       fft.ccf_df),
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