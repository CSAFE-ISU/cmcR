#' Divides an image into cells
#' @name cellDivision
#'
#' @param surfaceMat surface matrix of a breech face impression
#' @param horizSplits number of splits along horizontal axis
#' @param vertSplits number of splits along vertical axis
#'
#' @description This is a helper for the cellCCF function.

cellDivision <- function(surfaceMat,
                         horizSplits = 8,
                         vertSplits = 8){

  splitSurfaceMat <- surfaceMat %>%
    imager::as.cimg() %>%
    imager::imsplit(axis = "x",
                    nb = horizSplits) %>%
    purrr::map(.f = ~ imager::imsplit(.x,
                                      axis = "y",
                                      nb = vertSplits)) %>%
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
                                    rotatedSurfaceMat){
  #perform the appropriate subsetting of image A to create a list of larger
  #cells than those in image B
  splitRotatedSurfaceMat <- rotatedSurfaceMat[cornerLocs[["top.y"]]:cornerLocs[["bottom.y"]],
                                              cornerLocs[["left.x"]]:cornerLocs[["right.x"]]]

  return(splitRotatedSurfaceMat)
}

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

checkForBreechface <- function(cell,
                               bfMinimumProp = .15){
  containsBreechfaceBool <- (sum(!is.na(cell)) > bfMinimumProp*length(as.vector(cell)))
  return(containsBreechfaceBool)
}

#' Subtract the average pixel value from an image and replace NAs with 0
#' @name shiftSurfaceMatbyMean
#'
#' @param surfaceMat surface matrix
#' @param mean average pixel value to shift an image by
#'
#' @description This is a helper for the cellCCF function.
shiftSurfaceMatbyMean <- function(surfaceMat,m){
  surfaceMat <- surfaceMat - m
  surfaceMat[is.na(surfaceMat)] <- 0
  return(surfaceMat)
}

splitSurfaceMat1 <- function(surfaceMat,horizSplits,vertSplits){

  surfaceMat_split <- cellDivision(surfaceMat,
                                   horizSplits = horizSplits,
                                   vertSplits = vertSplits) #split image 1 into cells

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
                     checkForBreechface) %>%
    unlist()

  return(list(surfaceMat_split = surfaceMat_split,
              cellIDs = cellIDs,
              cellSideLengths = cellSideLengths,
              mat1PixCounts = matPixCounts))
}

getMat2SplitLocations <- function(cellIDs,
                                  cellSideLengths,
                                  mat2Dim){
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
                    c("left" = floor(xyLoc["x"] - sideLength["col"]),
                      "right" = ceiling(xyLoc["x"] + sideLength["col"]),
                      "top" = floor(xyLoc["y"] - sideLength["row"]),
                      "bottom" = ceiling(xyLoc["y"] + sideLength["row"]))

                  #replace negative indices with 0 (left/upper-most cells):
                  expandedCellCorners[expandedCellCorners <= 0] <- 1
                  #replace indices greater than the maximum index with the maximum index
                  #(right/bottom-most cells):
                  #Note that imager treats the rows of a matrix as the "x" axis and the columns as the "y" axis, contrary to intuition. As such, we need to swap the dimensions for when we subset the image further down in the function
                  if(expandedCellCorners[c("right.x")] > mat2Dim[1]){
                    expandedCellCorners[c("right.x")] <- mat2Dim[1]
                  }
                  if(expandedCellCorners[c("bottom.y")] > mat2Dim[2]){
                    expandedCellCorners[c("bottom.y")] <- mat2Dim[2]
                  }

                  return(expandedCellCorners)
                }) %>%
    setNames(cellIDs)

  return(mat2_splitCorners)
}

#' Calculate the maximum correlation between two breech face impressions split
#' into cells for range of rotation values
#' @name cellCCF
#' @export
#'
#' @param mat1 a matrix representing the surface matrix of a breech face
#'   impression
#' @param mat2 a matrix representing the surface matrix of a breech face
#'   impression to be compared to mat1
#' @param thetas rotation values (in degrees) for which mat2 will be rotated,
#'   split into cells, and compared to mat1
#' @param horizSplits number of splits along horizontal axis to divide mat1 into
#' @param vertSplits number of splits along vertical axis to divide mat1 into
#'
#' @return The list allResults contains the CCF values, horizontal, and vertical
#'   translations calculated for each cell, for each rotation value. The data
#'   frame topResults contains the CCF value, associated horizontal/vertical
#'   translation, and rotation value at which each cell in mat1 achieved its
#'   highest CCF with its paired cell in mat2.
#'
#' @description (Move this description elsewhere) This is an R implementation of the Congruent Matching Cells
#'   algorithm as described in 'Proposed "Congruent Matching Cells (CMC) Method
#'   for Ballistic Identification and Error Rate Estimation' by John Song
#'   (2015). The method works by first dividing mat1, the surface matrix
#'   representing the height values of a breech face impression microscopy scan,
#'   into pairwise disjoint cells. mat2, the surface matrix of a breech face
#'   impression scan to be compared to mat1, is also broken up into cells
#'   centered at the same location as those in mat1. However, the cells in mat2
#'   are larger than those in mat1 (typically 4 times the size except on the
#'   border of mat2). The cross-correlation function (CCF) is then calculated
#'   between each cell in mat1 and its larger, paired cell in mat2 using a Fast
#'   Fourier Transform. For each mat1, mat2 cell pair, the maximum CCF value as
#'   well as the necessary horizontal and vertical translations to align the mat1
#'   cell with the mat2 cell to achieve this max CCF are recorded. mat2 is then
#'   rotated by some amount (say, 2.5 degrees) and the process of dividing mat2
#'   into cells and calculating the max CCF for each cell in mat1 is repeated.
#'   The rotations performed on mat2 are dictated by the value(s) passed to the
#'   theta argument.
#'
#' @seealso
#' \url{https://pdfs.semanticscholar.org/4bf3/0b3a23c38d8396fa5e0d116cba63a3681494.pdf}
#'

cellCCF <- function(mat1,
                    mat2,
                    thetas = seq(-30,30,by = 3),
                    horizSplits = 7,
                    vertSplits = 7){
  #Needed tests:
  # mat1 and mat2 must be matrices
  # thetas, horizsplits, and vertsplits should be integers (horizsplits and vertsplits would optimally be equal - maybe print a warning if not?)

  mat1_mean <- mean(as.vector(mat1),na.rm = TRUE)
  mat2_mean <- mean(as.vector(mat2),na.rm = TRUE)

  #initialize list to contain all CCF values for each theta value
  allResults <- purrr::map(thetas,
                           function(theta) theta = data.frame(cell_ID = NA,
                                                              corr = NA,
                                                              dx = NA,
                                                              dy = NA)) %>%
    setNames(thetas)

  mat1_split <- splitSurfaceMat1(surfaceMat = mat1,
                                 horizSplits = horizSplits,
                                 vertSplits = vertSplits)

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
      rotateSurfaceMatrix(theta) %>%
      t() #because cellDivision starts in the left hand corner of image A and
    #moves down in tiling, we need to transpose image B so that the
    #correct cells match up (extractCellbyCornerLocs starts in the
    #left-hand corner but moves right.)

    #Now that we've rotated image B, we want to create a list consisting of the
    #cells that we are to compare to image A's cells. These will be 4 times
    #larger than Image A's cells (up to image B's boundary) and defined by the expandedCellCorners x,y
    #pairs calculated above.
    mat2_splitRotated <-
      purrr::map(mat2_splitCorners,
                 ~ extractCellbyCornerLocs(cornerLocs = .,
                                           rotatedSurfaceMat = mat2_rotated))

    #Determine whether a cell contains more than 10 pixels of breechface
    #information (for some reason checking for more than 0 pixels of
    #breechface was causing the function to fail)
    mat2PixCounts <- mat2_splitRotated %>%
      purrr::map(checkForBreechface) %>%
      unlist()

    #remove cells that don't include enough breechface:
    mat1_splitFiltered <- purrr::flatten(mat1_split$surfaceMat_split)[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]
    mat2_splitFiltered <- mat2_splitRotated[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #grab the cell IDs for each cell not removed above. This is used to update
    #topResults below
    filteredCellID <- mat1_split$cellIDs[mat1_split$mat1PixCounts == TRUE & mat2PixCounts == TRUE]

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs)
    im1_splitShifted <- mat1_splitFiltered %>%
      purrr::map(~ shiftSurfaceMatbyMean(.,mat1_mean)) %>%
      setNames(filteredCellID) #Need to set the names to avoid repeated list
    #labels

    #shift the pixel values in each image so that they both have 0 mean. Then
    #replace the NA values with 0 (FFTs can't deal with NAs)
    im2_splitShifted <- mat2_splitFiltered %>%
      purrr::map(~ shiftSurfaceMatbyMean(.,mat2_mean)) %>%
      setNames(filteredCellID)

    #calculate the correlation of each cell pair
    corrValues <- purrr::map2_dfr(.x = im1_splitShifted,
                                  .y = im2_splitShifted,
                                  .f = ~ data.frame(purrr::flatten(comparison(.x,t(.y))))) %>% #returns a nested list of corr,dx,dy value
      dplyr::mutate(cellID = filteredCellID)

    allResults[paste0(theta)][[1]] <- corrValues
  }

  allResults <- allResults %>%
    purrr::map(~ dplyr::select(.,cellID,corr,dx,dy)) #rearrange columns in allResults)

  return(allResults)
}

# cellCCF_topResults <- function(ccfAllResults){
#Initialize df to contain maximum CCF values per cell
# topResults <- data.frame(cellID = mat1_split$cellIDs,
#                          CCFmax = 0,
#                          dx = NA,
#                          dy = NA,
#                          theta = NA)
#   for(iden in filteredCellID){
#     oldRow <- dplyr::filter(topResults,cellID == iden)
#     comparisonRow <- dplyr::filter(corrValues,cellID == iden) %>%
#       dplyr::select(-cellID) #corr,dx,dy values for a given cell on this rotation value
#
#     if(!purrr::is_empty(comparisonRow$corr) &
#        oldRow$CCFmax < comparisonRow$corr){
#       #replace CCFmax,dx,dy,theta with new maximum values for each cell
#       topResults[which(topResults$cell_ID == iden),2:5] <- c(comparisonRow, theta)
#     }
#   }
# }