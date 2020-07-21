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
                               minObservedProp = .15){
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
  cellIDs <- purrr::map(names(surfaceMat_split),
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
              cellIDs = cellIDs,
              cellSideLengths = cellSideLengths,
              mat1PixCounts = matPixCounts))
}

#' @name getMat2SplitLocations
#'
#' @keywords internal
#'
#' @importFrom stats setNames

utils::globalVariables(c("."))

getMat2SplitLocations <- function(cellIDs,
                                  cellSideLengths,
                                  mat2Dim,
                                  sidelengthMultiplier,
                                  ...){
  mat2_splitCorners <- cellIDs %>%
    #pull all numbers from cellID strings:
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
    setNames(cellIDs)

  return(mat2_splitCorners)
}

#' @name swapCellIDAxes
#'
#' @keywords internal

swapCellIDAxes <- function(cellID){
  sSplit <- stringr::str_split(string = cellID,pattern = ",",n = 2)

  paste0(stringr::str_replace(string = sSplit[[1]][1],pattern = "x",replacement = "y"),
         ", ",
         stringr::str_replace(string = sSplit[[1]][2],pattern = "y",replacement = "x"))
}

#' @name calcRawCorr
#'
#' @description Given the dy,dx values at which CCF_max occurs, this function
#'   extract from the larger "region" matrix a matrix of the same dimension as
#'   "cell" with appropriate shifting/scaling relative to the center based on
#'   the dx,dy arguments. It is possible that the "cell" shaped matrix, after
#'   shifting by dx,dy, lies outside of the bounds of the region matrix (this is
#'   because both "cell" and "region" are padded with 0s when calculating the
#'   CCF). If this occurs, the "cell"-sized matrix that would be extracted from
#'   "region" is padded with NAs to appropriately reflect the conditions under
#'   which the CCF_max value was determined. Even with this padding, however,
#'   there may be a slight mis-match (like one row/col) between the dimensions
#'   of "cell" and of the matrix extracted from "region," typically because
#'   "region" has an odd dimension while "cell" has an even dimension or vice
#'   versa. Thus, additional padding/cropping needs to be performed in the
#'   "cell"-sized matrix extracted from "region." However, we don't know whether
#'   these rows/cols should pre/post padded/cropped, so all possible
#'   combinations of padding/cropping are considered. To determine which of
#'   these these combinations should be ultimately chosen as the final
#'   "cell"-sized matrix, the tieBreaker argument can be used to determine, for
#'   example, which "cell"-sized matrix has the highest correlation with "cell"
#'   (tieBreaker = which.max).
#'
#'
#' @keywords internal
#'
#' @importFrom stats cor

calcRawCorr <- function(cell,
                        region,
                        dx,
                        dy,
                        m1,
                        m2,
                        sd1,
                        sd2,
                        tieBreaker = which.max,
                        use = "pairwise.complete.obs"){
  cell <- (cell - m1)/sd1
  region <- (region - m2)/sd2

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

  #build-in some contingencies in case the regionCropped isn't same size as
  #the cell. If the dimensions don't agree, then we need to consider copies
  #of regionCroppedInitial that have been had rows/cols (whichever
  #applicable) cropped from the top/bottom or left/right
  regionCroppedList <- list("initial" = regionCroppedInitial)

  if(any(dim(regionCroppedInitial) != dim(cell))){
    #these make sure that the indices that we've cropped the region by aren't
    #less than 1 or larger than the dimension of the region
    if(alignedRows[1] < 1){
      rowsToPad <- -1*alignedRows[1] + 1

      regionCroppedInitial <- rbind(matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedInitial)),
                                         regionCroppedInitial)

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedRows[2] > nrow(region)){
      rowsToPad <- alignedRows[2] - nrow(region)

      regionCroppedInitial <- rbind(regionCroppedInitial,
                                         matrix(NA,nrow = rowsToPad,ncol = ncol(regionCroppedInitial)))

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedCols[1] < 1){
      colsToPad <- -1*alignedCols[1] + 1

      regionCroppedInitial <- cbind(matrix(NA,nrow = nrow(regionCroppedInitial),ncol = colsToPad),
                                         regionCroppedInitial)

      regionCroppedList$initial <- regionCroppedInitial
    }
    if(alignedCols[2] > ncol(region)){
      colsToPad <- alignedCols[2] - ncol(region)

      regionCroppedInitial <- cbind(regionCroppedInitial,
                                         matrix(NA,nrow = nrow(regionCroppedInitial),ncol = colsToPad))

      regionCroppedList$initial <- regionCroppedInitial
    }

    #two copies if rows are off
    if(nrow(regionCroppedInitial) > nrow(cell) & ncol(regionCroppedInitial) == ncol(cell)){
      rowsToCrop <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedRowPre <- regionCroppedInitial[-rowsToCrop,]
      regionCroppedRowPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),]

      regionCroppedList$rowPre <- regionCroppedRowPre
      regionCroppedList$rowPost <- regionCroppedRowPost
    }

    #two copies if rows are off
    if(nrow(regionCroppedInitial) < nrow(cell) & ncol(regionCroppedInitial) == ncol(cell)){
      rowsToPad <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedRowPre <- rbind(matrix(NA,
                                          nrow = rowsToPad,
                                          ncol = ncol(regionCroppedInitial)),
                                   regionCroppedInitial)

      regionCroppedRowPost <- rbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedInitial)))

      regionCroppedList$rowPre <- regionCroppedRowPre
      regionCroppedList$rowPost <- regionCroppedRowPost
    }

    #2 copies if cols are off
    if(ncol(regionCroppedInitial) > ncol(cell) & nrow(regionCroppedInitial) == nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedInitial) - ncol(cell))

      regionCroppedColPre <- regionCroppedInitial[,-colsToCrop]
      regionCroppedColPost <- regionCroppedInitial[,-(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedList$colPre <- regionCroppedColPre
      regionCroppedList$colPost <- regionCroppedColPost
    }

    #2 copies if cols are off
    if(ncol(regionCroppedInitial) < ncol(cell) & nrow(regionCroppedInitial) == nrow(cell)){
      colsToPad <- abs(ncol(regionCroppedInitial) - ncol(cell))

      regionCroppedColPre <- cbind(matrix(NA,
                                          nrow = nrow(regionCroppedInitial),
                                          ncol = colsToPad),
                                   regionCroppedInitial)

      regionCroppedColPost <- cbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad))

      regionCroppedList$colPre <- regionCroppedColPre
      regionCroppedList$colPost <- regionCroppedColPost
    }

    #4 different copies if both dimensions are too large
    if(ncol(regionCroppedInitial) > ncol(cell) & nrow(regionCroppedInitial) > nrow(cell)){
      colsToCrop <- abs(ncol(regionCroppedInitial) - ncol(cell))

      rowsToCrop <- abs(nrow(regionCroppedInitial) - nrow(cell))

      regionCroppedBothPre <- regionCroppedInitial[-rowsToCrop,
                                                   -colsToCrop]

      regionCroppedRowPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),
                                                   -colsToCrop]

      regionCroppedColPost <- regionCroppedInitial[-rowsToCrop,
                                                   -(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedBothPost <- regionCroppedInitial[-(nrow(regionCroppedInitial) - rowsToCrop),
                                                    -(ncol(regionCroppedInitial) - colsToCrop)]

      regionCroppedList$bothPre <- regionCroppedBothPre
      regionCroppedList$rowPost <- regionCroppedRowPost
      regionCroppedList$colPost <- regionCroppedColPost
      regionCroppedList$bothPost <- regionCroppedBothPost
    }

    #4 different copies if both dimensions are too small
    if(ncol(regionCroppedInitial) < ncol(cell) & nrow(regionCroppedInitial) < nrow(cell)){
      colsToPad <- abs(ncol(regionCroppedInitial) - ncol(cell))

      rowsToPad <- abs(nrow(regionCroppedInitial) - nrow(cell))

      #both rows and cols pre-padded
      regionCroppedBothPre <- cbind(matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad),
                                    regionCroppedInitial)
      regionCroppedBothPre <- rbind(matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedBothPre)),
                                    regionCroppedBothPre)

      #rows post-padded and cols pre-padded
      regionCroppedRowPost <- cbind(matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad),
                                    regionCroppedInitial)
      regionCroppedRowPost <- rbind(regionCroppedRowPost,
                                    matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedRowPost)))

      #rows pre-padded and cols post-padded
      regionCroppedColPost <- cbind(regionCroppedInitial,
                                    matrix(NA,
                                           nrow = nrow(regionCroppedInitial),
                                           ncol = colsToPad))
      regionCroppedColPost <- rbind(matrix(NA,
                                           nrow = rowsToPad,
                                           ncol = ncol(regionCroppedColPost)),
                                    regionCroppedColPost)

      #rows and cols both post-padded
      regionCroppedBothPost <- cbind(regionCroppedInitial,
                                     matrix(NA,
                                            nrow = nrow(regionCroppedInitial),
                                            ncol = colsToPad))
      regionCroppedBothPost <- rbind(regionCroppedBothPost,
                                     matrix(NA,
                                            nrow = rowsToPad,
                                            ncol = ncol(regionCroppedBothPost)))

      regionCroppedList$bothPre <- regionCroppedBothPre
      regionCroppedList$rowPost <- regionCroppedRowPost
      regionCroppedList$colPost <- regionCroppedColPost
      regionCroppedList$bothPost <- regionCroppedBothPost
    }
  }

  #return numeric(0) if cor fails.
  corSafe <- purrr::safely(cor,otherwise = numeric(0),quiet = TRUE)

  corrVals <- purrr::map_dbl(regionCroppedList,function(croppedRegion){

    corVal <- corSafe(as.vector(cell),
                      as.vector(croppedRegion),
                      use = use)

    as.numeric(corVal[1])
  })

  maxCorr <- as.numeric(corrVals[tieBreaker(corrVals)])
  if(length(maxCorr) == 0){
    return(NA)
  }
  else(return(maxCorr))
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
#' @param rawCorrTieBreaker the way in which the "raw" correlation (see
#'   description) is calculated may require slight padding/cropping of the
#'   mat1-sized matrix extracted from mat2 to make their dimensions equal. This
#'   padding/cropping can occur to the initial or final rows/cols in the matrix,
#'   without a clear way to determine which is "correct." As such, all possible
#'   combinations of pre/post padding/cropping are considered (only if
#'   necessary). To determine a final mat1-sized matrix, rawCorrTieBreaker can
#'   be used to determine which yields the lowest/highest correlation with matt1
#'   (using rawCorrTieBreaker = which.min or which.max, respectively).
#' @param use argument to be passed to the cor function. Dictates how NAs are
#'   dealt with in computing the correlation.
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
#'   are larger than those in mat1 (typically 4 times the size except on the
#'   border of mat2).
#'
#'   The cross-correlation function (CCF) is then calculated between each cell
#'   in mat1 and its larger, paired cell in mat2 using a Fast Fourier Transform.
#'   For each mat1, mat2 cell pair, the maximum CCF value as well as the
#'   necessary horizontal and vertical translations to align the mat1 cell with
#'   the mat2 cell to achieve this max CCF are recorded. Note that the FFT-based
#'   max CCF values aren't necessarily reliable estimates of the true similarity
#'   between mat1 and mat2 (since NAs need to be replaced prior to calculation).
#'   As such, a "raw" correlation is calculated by extracting from mat2 a
#'   mat1-sized matrix based on the dx, dy translation values obtained from the
#'   FFT-based CCF method. Rather than replacing NA values, these two matrices
#'   of equivalent dimension are turned into vectors and the pairwise-complete
#'   observation correlation (in which only non-NA pairs are considered) is
#'   calculated.
#'
#'   mat2 is then rotated by some amount (say, 2.5 degrees) and the process of
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

utils::globalVariables(c("cellID","."))

cellCCF <- function(x3p1,
                    x3p2,
                    thetas = seq(-30,30,by = 3),
                    cellNumHoriz = 8,
                    cellNumVert = cellNumHoriz,
                    regionToCellProp = 4,
                    minObservedProp = .15,
                    rawCorrTieBreaker = which.max,
                    use = "pairwise.complete.obs",
                    centerCell = "individualCell",
                    scaleCell = "individualCell"){
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
  #positions between two comparisons (since their "cellID" may differ since
  #the scans aren't the same size)
  cellIDdf <- data.frame(cellNum = seq_along(mat1_split$cellIDs),
                         cellID = mat1_split$cellIDs)

  #Now we want to split image B into cells with the same centers as those in
  #image A, but with twice the side length (these wider cells will intersect
  #each other). We will first get the dimensions (the x,y locations of each
  #cells' corners) of where each cell should be in image B.
  sidelengthMultiplier <- floor(sqrt(regionToCellProp))

  mat2_splitCorners <- getMat2SplitLocations(cellIDs = mat1_split$cellIDs,
                                             cellSideLengths = mat1_split$cellSideLengths,
                                             mat2Dim = dim(mat2),
                                             sidelengthMultiplier = sidelengthMultiplier)

  for(theta in thetas){
    #rotate image 2 and split into cells:

    mat2_rotated <- mat2 %>%
      rotateSurfaceMatrix(theta)

    #Now that we've rotated image B, we want to create a list consisting of the
    #cells that we are to compare to image A's cells. These will be 4 times
    #larger than Image A's cells (up to image B's boundary) and defined by the
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
    #which matrices need to be padded based on whether the associated cellID
    #contains 1 or the maximum value of the cellIDs (containing either part of
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
    filteredCellID <- mat1_split$cellIDs[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #creating this df is necessary so that we can compare cells in similar
    #positions between two comparisons (since their "cellID" may differ since
    #the scans aren't the same size, but their absolute position in the grid
    #(e.g., 4th from the left in the top row) will remain the same)
    cellIDdf_filtered <- cellIDdf %>%
      dplyr::filter(cellID %in% filteredCellID) %>%
      #imager swaps x and y axes, so we need to swap them back to be more
      #interpretable:
      dplyr::mutate(cellID = purrr::map_chr(cellID,swapCellIDAxes))

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs)
    if(!missing(centerCell)){
      if(centerCell == "individualCell"){
        m1 <- mat1_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)

        m2 <-  mat2_splitFiltered %>%
          purrr::map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    if(!missing(scaleCell)){
      if(scaleCell == "individualCell"){
        sd1 <-  mat1_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)

        sd2 <-  mat2_splitFiltered %>%
          purrr::map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    stopifnot(all(sd1 > 0),all(sd2 > 0))

    mat1_splitShifted <- purrr::pmap(list(mat1_splitFiltered,m1,sd1),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3)) %>%
      setNames(filteredCellID) #Need to set the names to avoid repeated list
    #labels

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs).
    mat2_splitShifted <- purrr::pmap(list(mat2_splitFiltered,m2,sd2),
                                     ~ standardizeSurfaceMat(surfaceMat = ..1,
                                                             m = ..2,
                                                             s = ..3)) %>%
      setNames(filteredCellID) #Need to set the names to avoid repeated list
    #labels

    #calculate the correlation of each cell pair
    ccfValues <- purrr::map2_dfr(.x = mat1_splitShifted,
                                 .y = mat2_splitShifted,
                                 .f = ~ data.frame(purrr::flatten(ccfComparison(.x,.y)))) %>% #returns a nested list of ccf,dx,dy values
      dplyr::bind_cols(cellIDdf_filtered,.) %>%
      dplyr::mutate(ccf = purrr::pmap_dbl(.l = list(mat1_splitFiltered,
                                                    mat2_splitFiltered,
                                                    .$dx,
                                                    .$dy,
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
                                                             use = use)),
                    theta = rep(theta,times = nrow(.)),
                    nonMissingProportion = purrr::map_dbl(mat1_splitFiltered,
                                                          .f = ~ (sum(!is.na(.)))/prod(dim(.))))

    allResults[paste0(theta)][[1]] <- ccfValues
  }

  allResults <- allResults %>%
    purrr::map(~ dplyr::select(.,cellNum,cellID,ccf,fft.ccf,dx,dy,theta,nonMissingProportion)) #rearrange columns in allResults

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
                    "mat2ScaleFactor" = sd2),
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
#' @param rawCorrTieBreaker the way in which the "raw" correlation (see
#'   description) is calculated may require slight padding/cropping of the
#'   mat1-sized matrix extracted from mat2 to make their dimensions equal. This
#'   padding/cropping can occur to the initial or final rows/cols in the matrix,
#'   without a clear way to determine which is "correct." As such, all possible
#'   combinations of pre/post padding/cropping are considered (only if
#'   necessary). To determine a final mat1-sized matrix, rawCorrTieBreaker can
#'   be used to determine which yields the lowest/highest correlation with matt1
#'   (using rawCorrTieBreaker = which.min or which.max, respectively).
#' @param use argument to be passed to the cor function. Dictates how NAs are
#'   dealt with in computing the correlation.
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
                                   minObservedProp = .15,
                                   rawCorrTieBreaker = which.max,
                                   use = "pairwise.complete.obs",
                                   centerCell = "individualCell",
                                   scaleCell = "individualCell"){



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
                             scaleCell = scaleCell)

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
                             scaleCell = scaleCell)

  return(list("comparison_1to2" = comparison_1to2,
              "comparison_2to1" = comparison_2to1))
}