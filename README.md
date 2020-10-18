
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cmcR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jzemmels/cmcR.svg?branch=master)](https://travis-ci.com/jzemmels/cmcR)
[![Codecov test
coverage](https://codecov.io/gh/jzemmels/cmcR/branch/master/graph/badge.svg)](https://codecov.io/gh/jzemmels/cmcR?branch=master)
<!-- badges: end -->

The cmcR package provides an open-source implementation of the Congruent
Matching Cells method for cartridge case identification as proposed by
[Song
(2013)](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911193)
as well as the “High CMC” method proposed by [Tong et
al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4730689/pdf/jres.120.008.pdf).

## Installation

<!-- You can install the released version of cmcR from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("cmcR") -->

<!-- ``` -->

Install the development version from
[GitHub](https://github.com/jzemmels/cmcR) with:

``` r
# install.packages("devtools")
devtools::install_github("jzemmels/cmcR")
```

Cartridge case scan data can be accessed at the [NIST Ballisitics and
Toolmarks Research
Database](https://tsapps.nist.gov/NRBTD/Studies/Search)

## Example

We will illustrate the package’s functionality here. This is intended to
be a concise demonstration, so please refer to the package vignettes
available under the “Articles” tab of the [package
website](https://csafe-isu.github.io/cmcR/index.html) for more detailed
information.

``` r
library(cmcR)
library(magrittr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

Consider the known match cartridge case pair Fadul 1-1 and Fadul 1-2.
The `read_x3p` function from the
[x3ptools](https://github.com/heike/x3ptools) package can read scans
from the [NBTRD](https://tsapps.nist.gov/NRBTD/Studies/Search) given the
appropriate address. The two scans are read below and visualized using
the
[`x3pListPlot`](https://csafe-isu.github.io/cmcR/reference/x3pListPlot.html)
function.

``` r
fadul1.1 <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/2d9cc51f-6f66-40a0-973a-a9292dbee36d")

fadul1.2 <- x3ptools::read_x3p("https://tsapps.nist.gov/NRBTD/Studies/CartridgeMeasurement/DownloadMeasurement/cb296c98-39f5-46eb-abff-320a2f5568e8")

cmcR::x3pListPlot(list("Fadul 1-1" = fadul1.1,
                       "Fadul 1-2" = fadul1.2),
                  type = "faceted")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

### Preprocessing

To perform a proper comparison of these two cartridge cases, we need to
remove regions that do not come into uniform or consistent contact with
the breech face of the firearm. These include the small clusters of
pixels in the corners of the two scans, caused by the staging area in
which the scans are taken, and the plateaued region of points around the
firing pin impression hole near the center of the scan. A variety of
processing procedures are implemented in the cmcR package. Consider the
[funtion
reference](https://csafe-isu.github.io/cmcR/reference/index.html) of the
cmcR package for more information regarding these procedures. As is
commonly done when comparing cartridge cases, we first downsample each
scan (by a factor of 4, selecting every other row/column) using the
`sample_x3p` function.

``` r
fadul1.1_processed <- fadul1.1 %>%
  cmcR::preProcess_cropBFExterior(radiusOffset = -30) %>%
  cmcR::preProcess_filterBFInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeBFTrend() %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

fadul1.2_processed <- fadul1.2 %>%
  cmcR::preProcess_cropBFExterior(radiusOffset = -30) %>%
  cmcR::preProcess_filterBFInterior(radiusOffset = 200) %>%
  cmcR::preProcess_removeBFTrend() %>%
  cmcR::preProcess_gaussFilter() %>%
  x3ptools::sample_x3p()

cmcR::x3pListPlot(list("Processed Fadul 1-1" = fadul1.1_processed,
                       "Processed Fadul1-2" = fadul1.2_processed),
                  type = "faceted")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Cell-based comparison procedure

Functions of the form `comparison_*` perform the steps of the cell-based
comparison procedure. The data generated from the cell-based comparison
procedure are kept in a [`tibble`](https://tibble.tidyverse.org/) where
one row represents a single cell/region pairing.

The `comparison_cellDivision` function divides a scan up into a grid of
cells. The `cellIndex` column represents the `row,col` location in the
original scan each cell inhabits. Each cell is stored as an `.x3p`
object in the `cellHeightValues` column. We will want to remove cells
that contain too many missing values. The benefit of using a `tibble`
structure is that processes such as removing rows can be accomplished
using simple `dplyr` commands such as `filter`.

``` r
cellTibble <- fadul1.1_processed %>%
  comparison_cellDivision(numCells = 64)

cellTibble
#> # A tibble: 64 x 2
#>    cellIndex cellHeightValues
#>    <chr>     <named list>    
#>  1 1, 1      <x3p>           
#>  2 1, 2      <x3p>           
#>  3 1, 3      <x3p>           
#>  4 1, 4      <x3p>           
#>  5 1, 5      <x3p>           
#>  6 1, 6      <x3p>           
#>  7 1, 7      <x3p>           
#>  8 1, 8      <x3p>           
#>  9 2, 1      <x3p>           
#> 10 2, 2      <x3p>           
#> # ... with 54 more rows
```

The `comparison_getTargetRegions` function extracts a region from a
target scan (in this case Fadul 1-2) to be paired with each cell in the
reference scan. The `regionSizeMultiplier` argument controls how much
larger the region is than its associated cell (without exceeding the
boundaries of the target scan).

``` r
cellTibble <- cellTibble %>%
  mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = cellHeightValues,
                                                          target_x3p = fadul1.2_processed,
                                                          regionSizeMultiplier = 9))

cellTibble
#> # A tibble: 64 x 3
#>    cellIndex cellHeightValues regionHeightValues
#>    <chr>     <named list>     <named list>      
#>  1 1, 1      <x3p>            <x3p>             
#>  2 1, 2      <x3p>            <x3p>             
#>  3 1, 3      <x3p>            <x3p>             
#>  4 1, 4      <x3p>            <x3p>             
#>  5 1, 5      <x3p>            <x3p>             
#>  6 1, 6      <x3p>            <x3p>             
#>  7 1, 7      <x3p>            <x3p>             
#>  8 1, 8      <x3p>            <x3p>             
#>  9 2, 1      <x3p>            <x3p>             
#> 10 2, 2      <x3p>            <x3p>             
#> # ... with 54 more rows
```

We will want to exclude cells and regions that contain few observations.
The `comparison_calcPropMissing` function calculates the proportion of
missing values in a surface matrix.

``` r
cellTibble <- cellTibble %>%
  mutate(cellPropMissing = comparison_calcPropMissing(cellHeightValues),
         regionPropMissing = comparison_calcPropMissing(regionHeightValues)) %>%
  filter(cellPropMissing <= .85 & regionPropMissing <= .85)

cellTibble %>%
  select(cellIndex,cellPropMissing,regionPropMissing)
#> # A tibble: 25 x 3
#>    cellIndex cellPropMissing regionPropMissing
#>    <chr>               <dbl>             <dbl>
#>  1 1, 6               0.833              0.808
#>  2 2, 7               0.657              0.700
#>  3 2, 8               0.834              0.757
#>  4 3, 8               0.353              0.638
#>  5 4, 8               0.153              0.576
#>  6 5, 1               0.117              0.767
#>  7 5, 8               0.0441             0.511
#>  8 6, 1               0.305              0.687
#>  9 6, 2               0.368              0.605
#> 10 6, 7               0.243              0.390
#> # ... with 15 more rows
```

We can also standardize the surface matrix height values by
centering/scaling by desired statistics. Also, to apply frequency-domain
techniques in comparing each cell and region, the missing values in each
scan need to be replaced. These operations are performed in the
`comparison_standardizeHeightValues` and
`comparison_replacingMissingValues` functions.

Then, the `comparison_fft.ccf` function estimates the translations
required to align the cell and region using the
[https://mathworld.wolfram.com/Cross-CorrelationTheorem.html](Cross-Correlation%20Theorem).
The `comparison_fft.ccf` function returns a data frame of 3 `x`, `y`,
and `fft.ccf` values: the \(x,y\) estimated translation values at which
the CCF\(_\max\) value is attained between the cell and region. The
`tidyr::unnest` function can unpack the data frame into 3 separate
columns, if desired.

``` r
cellTibble <- cellTibble  %>%
  mutate(cellHeightValues = comparison_standardizeHeightValues(cellHeightValues),
         regionHeightValues = comparison_standardizeHeightValues(regionHeightValues)) %>%
  mutate(cellHeightValues_replaced = comparison_replacingMissingValues(cellHeightValues),
         regionHeightValues_replaced = comparison_replacingMissingValues(regionHeightValues)) %>%
  mutate(fft.ccf_df = comparison_fft.ccf(cellHeightValues = cellHeightValues_replaced,
                                         regionHeightValues = regionHeightValues_replaced))

cellTibble %>%
  tidyr::unnest(cols = fft.ccf_df) %>%
  select(cellIndex,fft.ccf,x,y)
#> # A tibble: 25 x 4
#>    cellIndex fft.ccf     x     y
#>    <chr>       <dbl> <dbl> <dbl>
#>  1 1, 6        0.244    -7   -23
#>  2 2, 7        0.190    55    31
#>  3 2, 8        0.189    53    53
#>  4 3, 8        0.168   -23   -60
#>  5 4, 8        0.164    -3    13
#>  6 5, 1        0.213   -14   -49
#>  7 5, 8        0.139     7   -48
#>  8 6, 1        0.327   -53   -79
#>  9 6, 2        0.283     8   -45
#> 10 6, 7        0.169    13   -85
#> # ... with 15 more rows
```

Because so many missing values need to be replaced, the CCF\(_{\max}\)
value calculated in the `fft.ccf` column using frequency-domain
techniques is not a very good similarity score (doesn’t differentiate
matches from non-matches well). However, the `x` and `y` estimated
translations are good estimates of the “true” translation values needed
to align the cell and region. To calculate a more accurate similarity
score, we can use the pairwise-complete correlation in which only pairs
of non-missing pixels are considered in the correlation calculation.
This provides a better similarity metric. The pairwise-complete
correlation can be calculated with the `comparison_cor` function.

``` r
cellTibble %>%
  mutate(pairwiseCompCor = comparison_cor(cellHeightValues,regionHeightValues,fft.ccf_df)) %>%
  tidyr::unnest(fft.ccf_df) %>%
  select(cellIndex,x,y,pairwiseCompCor)
#> # A tibble: 25 x 4
#>    cellIndex     x     y pairwiseCompCor
#>    <chr>     <dbl> <dbl>           <dbl>
#>  1 1, 6         -7   -23           0.510
#>  2 2, 7         55    31           0.336
#>  3 2, 8         53    53           0.474
#>  4 3, 8        -23   -60           0.344
#>  5 4, 8         -3    13           0.336
#>  6 5, 1        -14   -49           0.430
#>  7 5, 8          7   -48           0.271
#>  8 6, 1        -53   -79           0.616
#>  9 6, 2          8   -45           0.581
#> 10 6, 7         13   -85           0.453
#> # ... with 15 more rows
```

Finally, this entire comparison procedure is to be repeated over a
number of rotations of the target scan. The resulting data frame below
contains the features that are used in the decision-rule procedure

``` r
cellTibble <- fadul1.1_processed %>%
  comparison_cellDivision(numCells = 64)

comparison_allTogether <- function(theta){
  
  cellTibble %>%
  mutate(regionHeightValues = comparison_getTargetRegions(cellHeightValues = cellHeightValues,
                                                          target_x3p = fadul1.2_processed,
                                                          regionSizeMultiplier = 9,
                                                          rotation = theta)) %>%
    mutate(cellPropMissing = comparison_calcPropMissing(cellHeightValues),
           regionPropMissing = comparison_calcPropMissing(regionHeightValues)) %>%
    filter(cellPropMissing <= .85 & regionPropMissing <= .85) %>%
    mutate(cellHeightValues = comparison_standardizeHeightValues(cellHeightValues),
           regionHeightValues = comparison_standardizeHeightValues(regionHeightValues)) %>%
    mutate(cellHeightValues_replaced = comparison_replacingMissingValues(cellHeightValues),
           regionHeightValues_replaced = comparison_replacingMissingValues(regionHeightValues)) %>%
    mutate(fft.ccf_df = comparison_fft.ccf(cellHeightValues = cellHeightValues_replaced,
                                           regionHeightValues = regionHeightValues_replaced)) %>%
    mutate(pairwiseCompCor = comparison_cor(cellHeightValues,regionHeightValues,fft.ccf_df)) %>%
    tidyr::unnest(fft.ccf_df) %>%
    select(cellIndex,x,y,pairwiseCompCor) %>%
    mutate(theta = theta)
  
}

kmComparisonFeatures <- purrr::map_dfr(seq(-30,30,by = 3),comparison_allTogether)

kmComparisonFeatures
#> # A tibble: 527 x 5
#>    cellIndex     x     y pairwiseCompCor theta
#>    <chr>     <dbl> <dbl>           <dbl> <dbl>
#>  1 2, 7        -28   -54           0.336   -30
#>  2 3, 1        -21    37           0.732   -30
#>  3 3, 8        -17   -22           0.503   -30
#>  4 4, 1         -9    35           0.675   -30
#>  5 4, 8         -8   -26           0.529   -30
#>  6 5, 1          1    34           0.561   -30
#>  7 5, 8          0   -22           0.499   -30
#>  8 6, 1          8    33           0.765   -30
#>  9 6, 2         65    66           0.524   -30
#> 10 6, 7          7   -12           0.328   -30
#> # ... with 517 more rows
```

### Decision rule

Under construction…
