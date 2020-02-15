#plot the CCF distribution

#plot the selected regions in each cell indicated by the CCF


#' @name getHighCCFPairs
#'
#' Must provide the parameters under which the comparison was made in the form
#' of a list containing, in this order, (1) the grid of theta values, (2) the
#' number of horizontal splits, and (3) the number of vertical splits. For
#' example, list(theta = seq(-30,30,by = 3), horizSplits = 7,vertSplits = 7)

getHighCCFPairs <- function(x3p1,
                            x3p2,
                            cellCCF_output,
                            initialCMCsOnly = FALSE,...){

  cellCCF_results <- cellCCF_output$ccfResults
  params <- cellCCF_output$params

  #Get the initial "CMC candidate" cells from the CCF results
  if(initialCMCsOnly){
    topResults <- cellCCF_results %>%
      cmcR:::countInitialCMCs(...)
  }
  else{
    topResults <- cellCCF_results %>%
      cmcR:::topResultsPerCell()
  }

  if(nrow(topResults) == 0){
    return(list("topResults" = topResults,
                "highCCFPairs" = NA))
  }

  mat1_split <- x3p1$surface.matrix %>%
    cmcR:::splitSurfaceMat1(horizSplits = params[[2]],
                            vertSplits = params[[3]])

  mat2_splitCorners <- cmcR:::getMat2SplitLocations(cellIDs = mat1_split$cellIDs,
                                                    cellSideLengths = mat1_split$cellSideLengths,
                                                    mat2Dim = dim(x3p2$surface.matrix))

  topResults_perTheta <- topResults %>%
    split(f = .$theta)

  #Now we want to pair cells/regions from the two x3ps by the rotation theta value by which they achieved highest CCFs
  mat1_highCCFCells <- mat1_split$surfaceMat_split %>%
    purrr::flatten() %>%
    setNames(mat1_split$cellIDs)

  mat1_highCCFCells <- mat1_highCCFCells[names(mat1_highCCFCells) %in% topResults$cellID]

  mat2_highCCFRegions <- purrr::map(topResults_perTheta,
                                    .f = function(ccfResults){
                                      theta <- as.numeric(unique(ccfResults$theta))

                                      mat2_rotated <- cmcR:::rotateSurfaceMatrix(surfaceMat = x3p2$surface.matrix,
                                                                                 theta = theta) %>%
                                        t()

                                      mat2_splitRotated <-
                                        purrr::map(mat2_splitCorners,
                                                   ~ cmcR:::extractCellbyCornerLocs(cornerLocs = .,
                                                                                    rotatedSurfaceMat = mat2_rotated)) %>%
                                        setNames(mat1_split$cellIDs)

                                      return(mat2_splitRotated[names(mat2_splitRotated) %in% ccfResults$cellID])
                                    }) %>%
    purrr::flatten()

  #the order of the cmcCandidate cells in im2 will likely be out of the order we want them to match up with the cells in im1, so we need to reorder them:
  mat2_highCCFRegions <- mat2_highCCFRegions[match(names(mat1_highCCFCells),
                                                   names(mat2_highCCFRegions))]

  #Return the "CMC candidate" pairs with the theta values at which they were most similar
  highCCFPairs <- purrr::map2(mat1_highCCFCells,
                              mat2_highCCFRegions,
                              list) %>%
    setNames(paste0("cellID = ",topResults$cellID,", theta = ",topResults$theta))

  return(list("topResults" = topResults,
              "highCCFPairs" = highCCFPairs))
}

#' @name stampOutMat1ShapeFromMat2

stampOutMat1ShapeFromMat2 <- function(pair,
                                      dx,
                                      dy){

  cell1 <- pair[[1]]
  cell2 <- pair[[2]]

  cell1 <- as.matrix(cell1)

  centerCells <- floor(dim(cell2)/2)

  alignmentLoc <- centerCells - c(dy,dx)

  halfDim <-  cartridges3D:::round2((dim(cell1))/2, 0)

  alignedCell_topLeft <- alignmentLoc - halfDim + 1
  alignedCell_bottomRight <- alignmentLoc + halfDim - 1

  if(dim(cell1)[1] %% 2 == 0){ #number of rows is always off by 1 when nrow of cell1 is even
    alignedCell_bottomRight[1] <- alignedCell_bottomRight[1] + 1
  }

  if(dim(cell1)[2] %% 2 == 0){ #number of cols is always off by 1 when ncol of cell1 is even
    alignedCell_bottomRight[2] <- alignedCell_bottomRight[2] + 1
  }

  #The alignment location may cause parts of im1 to fall outside of im2, so we need to pad im2 with NAs so that the correlation between im2 and the aligned im1 can be succesfully calculated:
  if(alignedCell_topLeft[1] <= 0){
    topRowsToPad <- 1 - alignedCell_topLeft[1]

    rowPad <- matrix(rep(rep(NA,
                             times = ncol(cell2)),
                         times = topRowsToPad),
                     ncol = ncol(cell2))

    cell2 <- rbind(rowPad,cell2)

    alignedCell_topLeft[1] <- 1
    alignedCell_bottomRight[1] <- alignedCell_bottomRight[1] + topRowsToPad
  }

  if(alignedCell_topLeft[2] <= 0){
    leftColsToPad <- 1 - alignedCell_topLeft[2]

    colPad <- matrix(rep(rep(NA,
                             times = nrow(cell2)),
                         times = leftColsToPad),
                     nrow = nrow(cell2))

    cell2 <- cbind(colPad,cell2)

    alignedCell_topLeft[2] <- 1
    alignedCell_bottomRight[2] <- alignedCell_bottomRight[2] + leftColsToPad
  }

  if(alignedCell_bottomRight[1] > nrow(cell2)){
    bottomRowsToPad <- alignedCell_bottomRight[1] - nrow(cell2)

    rowPad <- matrix(rep(rep(NA,
                             times = ncol(cell2)),
                         times = bottomRowsToPad),
                     ncol = ncol(cell2))

    cell2 <- rbind(cell2,rowPad)
  }

  if(alignedCell_bottomRight[2] > ncol(cell2)){
    rightColsToPad <- alignedCell_bottomRight[2] - ncol(cell2)

    colPad <- matrix(rep(rep(NA,
                             times = nrow(cell2)),
                         times = rightColsToPad),
                     nrow = nrow(cell2))

    cell2 <- cbind(cell2,colPad)
  }

  cell2Cropped <- cell2[alignedCell_topLeft[1]:alignedCell_bottomRight[1],
                        alignedCell_topLeft[2]:alignedCell_bottomRight[2]]

  return(list(cell1,cell2Cropped))
}

#' @name extractMatchedCells

extractMatchedCells <- function(x3p1,
                                x3p2,
                                cellCCF_results,
                                initialCMCsOnly = FALSE,...){

  highCCFPairs_results <- cmcR:::getHighCCFPairs(x3p1,
                                                 x3p2,
                                                 cellCCF_results,
                                                 initialCMCsOnly = FALSE)

  topResults <- highCCFPairs_results$topResults
  highCCFPairs <- highCCFPairs_results$highCCFPairs

  alignedCells <- purrr::pmap(list(highCCFPairs,
                                   topResults$dx,
                                   topResults$dy),
                              cmcR:::stampOutMat1ShapeFromMat2)

  return(alignedCells)
}

#' @name calcRawCorr

calcRawCorr <- function(x3p1,
                        x3p2,
                        cellCCF_results,
                        initialCMCsOnly = FALSE,
                        use = "pairwise.complete.obs",...){

  highCCFCells <- cmcR:::extractMatchedCells(x3p1,
                                             x3p2,
                                             cellCCF_results,
                                             initialCMCsOnly = FALSE,...)

  rawCorr <- purrr::map_dbl(highCCFCells,
                            function(pair){
                              if(any(
                                all(is.na(as.vector(pair[[1]]))) | all(is.na(as.vector(pair[[2]])))
                                )){
                                return(NA)
                              }
                              else{
                                cor(as.vector(pair[[1]]),
                                    as.vector(pair[[2]]),
                                    use = use)
                              }
                            })

  return(rawCorr)
}

