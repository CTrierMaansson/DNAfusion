# DNAfusion <img src="vignettes/logo.png" width="150" align="right">

An R/Bioconductor package to evaluate BAM files and identify genefusions, such as EML4-ALK

> This package was created in order to increase the sensitivity of EML4-ALK detection from commercially available NGS products such as the AVENIO (Roche) pipeline. Paired-end sequencing of cfDNA generated BAM files can be used as input to discover EML4-ALK variants. This package was developed using position deduplicated BAM files generated with the AVENIO Oncology Analysis Software. These files are made using the AVENIO ctDNA surveillance kit and Illumina Nextseq 500 sequencing. This is a targeted hybridization NGS approach and includes ALK-specific but not EML4-specific probes. The package includes six functions. The output of the first function, `EML4_ALK_detection()`, is used to determine whether EML4-ALK is detected and serves as input for the next four  exploratory functions characterizing the EML4-ALK variant. The last function `EML4_ALK_analysis()` combines the output of the exploratory functions. To serve as examples, this package includes BAM files representing the EML4-ALK positive cell line H3122 and the EML4-ALK negative cell line, HCC827.

## Highlights
DNAfusion is under active development. In the future the package will include more genefusions such as ROS1 and RET, as well as other ALK fusionpartners.

## Installation

Use devtools to install the most recent version of DNAfusion from the GitHub repository.

```R
if (!require(devtools)) install.packages('devtools')
library(devtools)

devtools::install_github("CTrierMaansson/DNAfusion", build_vignettes = TRUE)
library(DNAfusion)

```

Alternatively, install DNAfusion published at Bioconductor.

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DNAfusion")
library(DNAfusion)

```

## Vignettes

See the documentation at https://github.com/CTrierMaansson/DNAfusion/tree/main/vignettes

Alternatively, the vignettes can be browsed in Rstudio with 

 ```R
browseVignettes("DNAfusion")
```
This requires `build_vignettes = TRUE` during installation with github

## Citation

Please cite package ‘DNAfusion’ in publications:

Christoffer Trier Maansson, Emma Roger Andersen, Maiken Parm Ulhoi, Peter Meldgaard, Boe Sandahl Sorensen (2022). DNAfusion: An R/Bioconductor package to identify gene fusions. R package version 1.1.2 
<https://github.com/CTrierMaansson/DNAfusion>
