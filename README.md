## DRIMSeq - differential transcript usage and transcript usage QTL analyses with a Dirichlet-multinomial model in RNA-Seq


`DRIMSeq` package provides two frameworks. One for the differential transcript usage (DTU) analysis between different designs and one for the tuQTL analysis. Both are based on modeling transcript counts with the Dirichlet-multinomial distribution. DTU analysis can be performed at the gene and/or transcript level. The package also makes available functions for visualization and exploration of the data and results.


# Bioconductor installation 

`DRIMSeq` is available on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/DRIMSeq.html) and can be installed with the command:

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DRIMSeq")
```

The vignette containing all the instructions on how to use the DRIMSeq package for DTU and tuQTL analyses can be accessed by entering:

``` r
browseVignettes("DRIMSeq")
```

or on the [Bioconductor website](https://www.bioconductor.org/packages/release/bioc/vignettes/DRIMSeq/inst/doc/DRIMSeq.pdf).


# Devel installation from Github


To install the latest development version, use the `devtools` package (available [here](https://github.com/hadley/devtools)):

``` r
devtools::install_github("markrobinsonuzh/DRIMSeq")
```

To install it with the vignette type:

``` r
devtools::install_github("markrobinsonuzh/DRIMSeq", build_vignettes = TRUE)
```

The vignette can be accessed from the R console by typing:

``` r
vignette("DRIMSeq")
```

or

``` r
browseVignettes("DRIMSeq")
```


# Data packages 

In order to run the examples from the vignette and the manual, you need to install two data packages [PasillaTranscriptExpr](https://www.bioconductor.org/packages/release/data/experiment/html/PasillaTranscriptExpr.html) and [GeuvadisTranscriptExpr](https://www.bioconductor.org/packages/release/data/experiment/html/GeuvadisTranscriptExpr.html).

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("PasillaTranscriptExpr")
biocLite("GeuvadisTranscriptExpr")
```


