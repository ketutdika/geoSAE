
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geoSAE

<!-- badges: start -->
<!-- badges: end -->

This function is an extension of the Small Area Estimation (SAE) model.
Geoadditive Small Area Model is a combination of the geoadditive model
with the Small Area Estimation (SAE) model, by adding geospatial
information to the SAE model.

## Authors

Ketut Karang Pradnyadika, Ika Yuni Wulansari

## Maintainer

Ketut Karang Pradnyadika <221709776@stis.ac.id>

## Installation

You can install the released version of geoSAE from
[CRAN](https://CRAN.R-project.org) or find my github repository
[Github](https://github.com/ketutdika)

## Example

``` r
library(geoSAE)

#Load the dataset for unit level
data(dataUnit)

#Load the dataset for spline-2
data(zspline)

#Load the dataset for area level
data(dataArea)

#Construct the data frame
y       <- dataUnit$y
x1      <- dataUnit$x1
x2      <- dataUnit$x2
x3      <- dataUnit$x3
formula <- y~x1+x2+x3
zspline <- as.matrix(zspline[,1:10])
dom     <- dataUnit$area
xmean   <- cbind(1,dataArea[,2:4])
zmean   <- dataArea[,5:14]
number  <- dataUnit$number
area    <- dataUnit$area
data    <- data.frame(number, area, y, x1, x2, x3)

#Estimate EBLUP
eblup_geosae <- eblupgeo(formula, zspline, dom, xmean, zmean, data)
eblup_geosae$eblup
#>           [,1]
#>  [1,] 28.80487
#>  [2,] 34.27508
#>  [3,] 35.32142
#>  [4,] 33.04335
#>  [5,] 23.91203
#>  [6,] 21.68321
#>  [7,] 20.40650
#>  [8,] 19.82616
#>  [9,] 34.74351
#> [10,] 39.99695
#> [11,] 35.20031
#> [12,] 28.46516
#> [13,] 31.94645
#> [14,] 24.12550
#> [15,] 30.97060
 
#Estimate MSE
mse_geosae <- pbmsegeo(formula,zspline,dom,xmean,zmean,data,B=100)
#> 
#> Bootstrap procedure with B = 100 iterations starts.
mse_geosae$mse
#>  [1]  3.9030125  3.6740805  2.4060875 11.4583413  1.8236213  5.8924137
#>  [7]  1.9357287  1.3473704  3.9215650  0.7674256  2.4733603 15.6073802
#> [13] 10.4960697  1.8853202  4.1279470
 
## eblup_geosae$eblup        #to see the result of EBLUPs with Geoadditive Small Area Model each area
## mse_geosae$mse            #to see the result of MSE with Geoadditive Small Area Model each area
```

## References

-   Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
-   Bocci, C., & Petrucci, A. (2016). Spatial information and
    geoadditive small area models. Analysis of poverty data by small
    area estimation, 245-259.
-   Ardiansyah, M., Djuraidah, A., & Kurnia, A. (2018). PENDUGAAN AREA
    KECIL DATA PRODUKTIVITAS TANAMAN PADI DENGAN GEOADDITIVE SMALL AREA
    MODEL. Jurnal Penelitian Pertanian Tanaman Pangan, 2(2), 101-110.
