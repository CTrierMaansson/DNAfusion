# DNAfusion <img src="vignettes/logo.png" width="150" align="right">
[![Anaconda-Server Badge](
https://anaconda.org/bioconda/bioconductor-dnafusion/badges/downloads.svg)](
https://anaconda.org/bioconda/bioconductor-dnafusion/badges/downloads.svg)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FCTrierMaansson%2FDNAfusion&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

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

If you use [DNAfusion](https://bioconductor.org/packages/release/bioc/html/DNAfusion.html) 
in published research please cite:

1.  Christoffer Trier Maansson, Emma Roger Andersen, Maiken Parm Ulhoi, Peter Meldgaard, Boe Sandahl Sorensen (2022). DNAfusion: An R/Bioconductor package to identify gene fusions. R package version 1.2.0 <https://github.com/CTrierMaansson/DNAfusion>

2.  Maansson CT, Andersen ER, Ulhoi MP, Meldgaard P, Sorensen BS. “DNAfusion: an
R/Bioconductor package for increased sensitivity of detecting gene fusions in liquid
biopsies.” BMC Bioinformatics, 2023, 24:131,
[doi:10.1186/s12859-023-05259-3](https://doi.org/10.1186/s12859-023-5259-3).
