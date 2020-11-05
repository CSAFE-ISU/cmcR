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
    dplyr::mutate(cellIndex = linear_to_matrix(index = (cellNum %% ceiling(sqrt(max(cellNum)))) +
                                                 floor((ceiling(sqrt(max(cellNum)))^2 - cellNum)/ceiling(sqrt(max(cellNum))))*ceiling(sqrt(max(cellNum))) +
                                                 ifelse(cellNum %% ceiling(sqrt(max(cellNum))) == 0,ceiling(sqrt(max(cellNum))),0),
                                               nrow = ceiling(sqrt(max(cellNum))),
                                               byrow = TRUE)) %>%
    dplyr::select(cellIndex,cellHeightValues)

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
#' @export

comparison_calcPropMissing <- function(heightValues){
  heightValues %>%
    purrr::map_dbl(~ sum(is.na(.$surface.matrix))/length(.$surface.matrix))
}

#' Extract regions from a target scan based on associated cells in reference
#' scan
#'
#' @name comparison_getTargetRegions
#' @param cellHeightValues list/tibble column of x3p objects containing a
#'   reference scan's cells (as returned by comparison_cellDivision)
#' @param target x3p object containing a breech face scan to be compared to
#'   the reference cell.
#' @param theta degrees that the target scan is to be rotated prior
#'   extracting regions.
#' @param regionSizeMultiplier ratio between the area of each target scan
#'   regions and the reference scan cells (e.g., 9 means that the regions'
#'   surface matrices will have thrice the number of rows and columns as the
#'   cells' surface matrices, 4 means twice the number rows and columns, etc.)
#'@examples
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
#' @export

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
                 regionMatrix <- extractCellbyCornerLocs(cornerLocs = cornerIndices,
                                                                rotatedSurfaceMat = target_surfaceMat_rotated,
                                                                mat2Dim = dim(target$surface.matrix))

                 region_x3p <- x3ptools::df_to_x3p(data.frame(x = 1,y = 1,value = NA))

                 region_x3p$surface.matrix <- regionMatrix

                 #update metainformation
                 region_x3p$header.info <- target$header.info
                 region_x3p$header.info$sizeY <- ncol(regionMatrix)
                 region_x3p$header.info$sizeX <- nrow(regionMatrix)

                 return(region_x3p)
               } )
}

#' Standardize height values of a scan by centering/scaling by desired
#' statistics and replacing missing values
#'
#' @name comparison_standardizeHeights
#'
#' @param heightValues list/tibble column of x3p objects
#' @param withRespectTo currently ignored
#' @param centerBy statistic by which to center (i.e., subtract from) the height
#'   values
#' @param scaleBy statistic by which to scale (i.e., divide) the height values
#'
#' @note this function adds information to the metainformation of the x3p scan
#'   it is given that is required for calculating, for example, the
#'   pairwise-complete correlation using the comparison_cor function.
#'@examples
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
#' @export

comparison_standardizeHeights <- function(heightValues,
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
#' @name comparison_replaceMissing
#' @param heightValues list/tibble column of x3p objects
#' @param replacement value to replace NAs
#'@examples
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
#' @export

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
#'@note The FFT is not defined for matrices containing missing values. The
#'  missing values in the cell and region need to be replaced before using this
#'  function. See the \link[cmcR]{comparison_replaceMissing} function to replace
#'  missing values after standardization.
#'@seealso \url{(https://mathworld.wolfram.com/Cross-CorrelationTheorem.html)}
#'
#'@return a data frame containing the translation (x,y) at which the CCF was
#'  maximized in aligning a target scan region to its associated reference scan
#'  cell.
#'
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
#'@export
comparison_fft_ccf <- function(cellHeightValues,regionHeightValues){
  ccfList <- purrr::map2(cellHeightValues,
                         regionHeightValues,
                         ~ ccfComparison(mat1 = .x$surface.matrix,
                                                mat2 = .y$surface.matrix,
                                                ccfMethod = "fft"))

  return(ccfList)
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
#'data(fadul1.1_processed,fadul1.2_processed)
#'
#'cellTibble <- comparison_allTogether(reference = fadul1.1_processed,target = fadul1.2_processed)
#'
#'head(cellTibble)
#'
#'@export

comparison_allTogether <- function(reference,
                                   target,
                                   theta = 0,
                                   numCells = 64,
                                   maxMissingProp = .85){

  reference %>%
    comparison_cellDivision(numCells) %>%
    dplyr::mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = cellHeightValues,
                                                            target = target,
                                                            theta = theta)) %>%
    dplyr::mutate(cellPropMissing = comparison_calcPropMissing(cellHeightValues),
           regionPropMissing = comparison_calcPropMissing(regionHeightValues)) %>%
    dplyr::filter(cellPropMissing <= maxMissingProp & regionPropMissing <= maxMissingProp) %>%
    dplyr::mutate(cellHeightValues = comparison_standardizeHeights(cellHeightValues),
           regionHeightValues = comparison_standardizeHeights(regionHeightValues)) %>%
    dplyr::mutate(cellHeightValues_replaced = comparison_replaceMissing(cellHeightValues),
           regionHeightValues_replaced = comparison_replaceMissing(regionHeightValues)) %>%
    dplyr::mutate(fft_ccf_df = comparison_fft_ccf(cellHeightValues = cellHeightValues_replaced,
                                           regionHeightValues = regionHeightValues_replaced)) %>%
    dplyr::mutate(pairwiseCompCor = comparison_cor(cellHeightValues,regionHeightValues,fft_ccf_df)) %>%
    tidyr::unnest(fft_ccf_df) %>%
    dplyr::select(cellIndex,x,y,fft_ccf,pairwiseCompCor) %>%
    dplyr::mutate(theta = theta)

}