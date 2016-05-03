## DRIMSeq package - Differential splicing and sQTL analyses with Dirichlet-multinomial model in RNA-Seq


`DRIMSeq` package provides two frameworks. One for the differential splicing analysis between different conditions and one for the sQTL analysis. Both are based on modeling counts of genomic features (i.e., transcripts, exons or exonic bins) with Dirichlet-multinomial distribution. The package also makes available functions for visualization and exploration of the data and results.


# Installation 

To install the latest development version, use the `devtool` package (available [here](https://github.com/hadley/devtools))

```
devtools::install_github("markrobinsonuzh/DRIMSeq")
```

# Vignette

The vignette contains all the instructions on how to use the DRIMSeq package for differential and sQTL analyses. It can be found in the "vignettes/DRIMSeq.pdf" file. After installation, the vignette can be accessed from the R console by typing

```
browseVignettes("DRIMSeq")
```

In order to run the examples from the vignette and the manual, you need to install two data packages [PasillaTranscriptExpr](https://github.com/markrobinsonuzh/PasillaTranscriptExpr) and [GeuvadisTranscriptExpr](https://github.com/markrobinsonuzh/GeuvadisTranscriptExpr).

