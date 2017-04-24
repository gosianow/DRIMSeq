## DRIMSeq package - Differential transcript usage and transcript usage QTL analyses with Dirichlet-multinomial model in RNA-Seq


`DRIMSeq` package provides two frameworks. One for the differential transcript usage (DTU) analysis between different designs and one for the tuQTL analysis. Both are based on modeling transcript counts with Dirichlet-multinomial distribution. DTU analysis can be performed at the gene and/or transcript level. The package also makes available functions for visualization and exploration of the data and results.


# Installation 

To install the latest development version, use the `devtools` package (available [here](https://github.com/hadley/devtools)):

```
devtools::install_github("markrobinsonuzh/DRIMSeq")
```

# Vignette

The vignette contains all the instructions on how to use the DRIMSeq package for DTU and tuQTL analyses. After installation, the vignette can be accessed from the R console by typing:

```
browseVignettes("DRIMSeq")
```

In order to run the examples from the vignette and the manual, you need to install two data packages [PasillaTranscriptExpr](https://github.com/markrobinsonuzh/PasillaTranscriptExpr) and [GeuvadisTranscriptExpr](https://github.com/markrobinsonuzh/GeuvadisTranscriptExpr).

