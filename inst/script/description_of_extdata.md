# A description of example files in this package.

## BAM files in inst/extdata

This folder contains PosDeduped BAM files and bam.bai files from two separate cell types (H3122, HCC827). 
These files are made using the AVENIO ctDNA surveillance kit and Illumina Nextseq 500 sequencing.
The BAM files were generated using the AVENIO Oncology Analysis Software with hg38 as reference genome. 
In order to reduce the file size, only reads in EML4 are included.
EML4 reads were filtered with samtools view input.bam "Chr2:42169353-42332548" > output.bam

