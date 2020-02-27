#' Divides an image into cells
#' @name cellDivision
#'
#' @param surfaceMat surface matrix of a breech face impression
#' @param cellNumHoriz number of splits along horizontal axis
#' @param cellNumVert number of splits along vertical axis
#'
#' @description This is a helper for the cellCCF function.

cellDivision <- function(surfaceMat,
                         cellNumHoriz = 8,
                         cellNumVert = 8,...){

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

rotateSurfaceMatrix <- function(surfaceMat,
                                theta = 0,
                                interpolation = 0,...){
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

checkForBreechface <- function(cell,
                               minObservedProp = .15,...){
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
standardizeSurfaceMat <- function(surfaceMat,
                                  m,
                                  s,...){
  surfaceMat <- (surfaceMat - m)/s
  surfaceMat[is.na(surfaceMat)] <- 0
  return(surfaceMat)
}

splitSurfaceMat1 <- function(surfaceMat,cellNumHoriz,cellNumVert,minObservedProp,...){

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

getMat2SplitLocations <- function(cellIDs,
                                  cellSideLengths,
                                  mat2Dim,...){
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
                    c(floor(xyLoc["y"] - sideLength["col"]),
                      ceiling(xyLoc["y"] + sideLength["col"]),
                      floor(xyLoc["x"] - sideLength["row"]),
                      ceiling(xyLoc["x"] + sideLength["row"])) %>%
                    setNames(c("left","right","top","bottom"))

                  #replace negative indices with 0 (left/upper-most cells):
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

swapCellIDAxes <- function(cellID){
  sSplit <- stringr::str_split(string = cellID,pattern = ",",n = 2)

  paste0(stringr::str_replace(string = sSplit[[1]][1],pattern = "x",replacement = "y"),
         ",",
         stringr::str_replace(string = sSplit[[1]][2],pattern = "y",replacement = "x"))
}

#' Calculate the maximum correlation between two breech face impressions split
#' into cells for range of rotation values
#' @name cellCCF
#' @export
#'
#' @param x3p1 (no default) an x3p object containing the surface matrix of a
#'   breech face impression
#' @param x3p2 (no default) an x3p object containing the surface matrix of a
#'   breech face impression to be compared to that in x3p1
#' @param thetas (default seq(from = -30,to = 30,by = 3)) rotation values (in
#'   degrees) for which x3p2$surface.matrix will be rotated, split into cells,
#'   and compared to x3p1$surface.matrix
#' @param cellNumHoriz (default 7) number of cells along horizontal axis to
#'   divide x3p1$surface.matrix into
#' @param cellNumVert (default equal to cellNumHoriz) number of cells along
#'   vertical axis to divide x3p1$surface.matrix into
#' @param minObservedProp (default 15) the minimum proportion of a cell that
#'   needs to contain observed (i.e., non-NA) values for it to be included in
#'   the CCF calculation procedure
#' @param centerCell **OPTIONAL** (default missing) dictates if cell is to be
#'   shifted prior to CCF calculation. Default is that no shifting is performed.
#'   If set to "wholeMatrix", then each cell is subracted by the mean of the
#'   whole matrix. If set to "individualCell", then each cell is subtracted by
#'   its particular mean.
#' @param scaleCell **OPTIONAL** (default missing) dictates if cell is to be
#'   scaled prior to CCF calculation. Default is that no scaling is performed.
#'   If set to "wholeMatrix", then each cell is divided by the standard
#'   deviation of the whole matrix. If set to "individualCell", then each cell
#'   is divided by its particular standard deviation.
#'
#' @return The list allResults contains the CCF values, horizontal, and vertical
#'   translations calculated for each cell, for each rotation value. The data
#'   frame topResults contains the CCF value, associated horizontal/vertical
#'   translation, and rotation value at which each cell in mat1 achieved its
#'   highest CCF with its paired cell in mat2.
#'
#' @description (Move this description elsewhere) This is an R implementation of
#'   the Congruent Matching Cells algorithm as described in 'Proposed "Congruent
#'   Matching Cells (CMC) Method for Ballistic Identification and Error Rate
#'   Estimation' by John Song (2015). The method works by first dividing mat1,
#'   the surface matrix representing the height values of a breech face
#'   impression microscopy scan, into pairwise disjoint cells. mat2, the surface
#'   matrix of a breech face impression scan to be compared to mat1, is also
#'   broken up into cells centered at the same location as those in mat1.
#'   However, the cells in mat2 are larger than those in mat1 (typically 4 times
#'   the size except on the border of mat2). The cross-correlation function
#'   (CCF) is then calculated between each cell in mat1 and its larger, paired
#'   cell in mat2 using a Fast Fourier Transform. For each mat1, mat2 cell pair,
#'   the maximum CCF value as well as the necessary horizontal and vertical
#'   translations to align the mat1 cell with the mat2 cell to achieve this max
#'   CCF are recorded. mat2 is then rotated by some amount (say, 2.5 degrees)
#'   and the process of dividing mat2 into cells and calculating the max CCF for
#'   each cell in mat1 is repeated. The rotations performed on mat2 are dictated
#'   by the value(s) passed to the theta argument.
#'
#' @seealso
#' \url{https://pdfs.semanticscholar.org/4bf3/0b3a23c38d8396fa5e0d116cba63a3681494.pdf}
#'

cellCCF <- function(x3p1,
                    x3p2,
                    thetas = seq(-30,30,by = 3),
                    cellNumHoriz = 7,
                    cellNumVert = cellNumHoriz,
                    minObservedProp = .15,
                    centerCell,
                    scaleCell,...){
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

  #Now we want to split image B into cells with the same centers as those in
  #image A, but with twice the side length (these wider cells will intersect
  #each other). We will first get the dimensions (the x,y locations of each
  #cells' corners) of where each cell should be in image B.
  mat2_splitCorners <- getMat2SplitLocations(cellIDs = mat1_split$cellIDs,
                                             cellSideLengths = mat1_split$cellSideLengths,
                                             mat2Dim = dim(mat2))

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
    mat1_splitFiltered <- purrr::flatten(mat1_split$surfaceMat_split)[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]
    mat2_splitFiltered <- mat2_splitRotated[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #grab the cell IDs for each cell not removed above. This is used to update
    #topResults below
    filteredCellID <- mat1_split$cellIDs[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs)
    if(!missing(centerCell)){
      if(centerCell == "individualCell"){
        m1 <- mat1_splitFiltered %>%
          map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)

        m2 <-  mat2_splitFiltered %>%
          map(~ mean(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

    if(!missing(scaleCell)){
      if(scaleCell == "individualCell"){
        sd1 <-  mat1_splitFiltered %>%
          map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)

        sd2 <-  mat2_splitFiltered %>%
          map(~ sd(.,na.rm = TRUE)) %>%
          setNames(filteredCellID)
      }
    }

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
                                 .f = ~ data.frame(purrr::flatten(cmcR:::comparison(.x,.y)))) %>% #returns a nested list of ccf,dx,dy values
      dplyr::mutate(cellID = filteredCellID %>%
                      purrr::map_chr(swapCellIDAxes)) #imager swaps x and y axes, so we need to swap them back to be more interpretable

    allResults[paste0(theta)][[1]] <- ccfValues
  }

  allResults <- allResults %>%
    purrr::map(~ dplyr::select(.,cellID,ccf,dx,dy)) #rearrange columns in allResults

  if(missing(centerCell)){
    centerCell <- "none"
  }
  if(missing(scaleCell)){
    scaleCell <- "none"
  }

  return(list(
    "params" = list("theta" = theta,
                    "cellNumHoriz" = cellNumHoriz,
                    "cellNumVert" = cellNumVert,
                    "minObservedProp" = minObservedProp,
                    "centerCell" = centerCell,
                    "mat1Shift" = m1,
                    "mat2Shift" = m2,
                    "scaleCell" = scaleCell,
                    "mat1ScaleFactor" = sd1,
                    "mat2ScaleFactor" = sd2),
    "ccfResults" = allResults
  ))
}

#' Wrapper for applying the cmcR::cellCCF function to x3p1 vs. x3p2 and again
#' for x3p2 vs. x3p1. See cellCCF function documentation for more details.
#' @name cellCCF_bothDirections
#'
#' @seealso cellCCF function
#' @export

cellCCF_bothDirections <- function(x3p1,
                                   x3p2,
                                   thetas = seq(-30,30,by = 3),
                                   cellNumHoriz = 7,
                                   cellNumVert = cellNumHoriz,
                                   centerCell,
                                   scaleCell,...){



  comparison_1to2 <- cmcR::cellCCF(x3p1 = x3p1,
                                   x3p2 = x3p2,
                                   thetas = thetas,
                                   cellNumHoriz = cellNumHoriz,
                                   cellNumVert = cellNumVert,
                                   centerCell = centerCell,
                                   scaleCell = scaleCell)

  comparison_2to1 <- cmcR::cellCCF(x3p1 = x3p2,
                                   x3p2 = x3p1,
                                   thetas = thetas,
                                   cellNumHoriz = cellNumHoriz,
                                   cellNumVert = cellNumVert,
                                   centerCell = centerCell,
                                   scaleCell = scaleCell)

  return(list("comparison_1to2" = comparison_1to2,
              "comparison_2to1" = comparison_2to1))
}