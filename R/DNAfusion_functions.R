#' @noRd
#' @importFrom GenomicAlignments cigar
#' @importFrom GenomicRanges mcols
index_helper <- function(input){
    cigars <- cigar(input)
    result <- sub("\\M.*", "", cigars)
    splits <- strsplit(result, split = "S")
    fun <- function(str){
        if(length(str) > 1){
            str <- NA
        }
        else{
            str <- str
        }
        splits <- lapply(str, as.numeric)
        return(splits)
    }
    index2 <- vapply(splits, FUN = fun, FUN.VALUE = list(1))
    GenomicRanges::mcols(input)["index"] <- unlist(index2)
    input <- input[!is.na(mcols(input)["index"])[,1],]
    return(input)
}
#' Detection of EML4-ALK variants
#'
#' This function looks for EML4-ALK mate pair reads in the BAM file.
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom GenomicAlignments readGAlignments cigar seqnames
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom IRanges IRanges
#' @importFrom BiocBaseUtils isScalarNumber isScalarCharacter
#' @param file The name of the file which the data are to be read from.
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param mates `Interger`, the minimum number EML4-ALK mate pairs
#' needed to be detected in order to call a variant. Default=2.
#' @return A `GAlignments` object with soft-clipped reads representing
#'  EML4-ALK is returned. If no EML4-ALK is detected the `GAlignments`
#'  is empty.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' EML4_ALK_detection(file=H3122_bam,
#'                     genome="hg38",
#'                     mates=2)
#' EML4_ALK_detection(file=HCC827_bam,
#'                     genome="hg38",
#'                     mates=2)
#' @export
EML4_ALK_detection <- function(file, genome="hg38", mates=2){
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!isScalarNumber(mates)){
        stop("mates has to be a numeric")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    what <- c("mpos", "seq")
    if (genome =="hg38"){
        which <- GRanges(seqnames="chr2",
                            IRanges(start=42169353, end=42332548))
        param <- ScanBamParam(which=which, what=what)
        reads <- readGAlignments(file=file, param=param)
        reads <- reads[(29192774 < mcols(reads)[,1] &
                            mcols(reads)[,1] < 29921586 &
                            !is.na(mcols(reads)[,1])),]
    }
    else{
        which <- GRanges(seqnames="chr2",
                            IRanges(start=42396490, end=42559688))
        param <- ScanBamParam(which=which, what=what)
        reads <- readGAlignments(file=file, param=param)
        reads <- reads[(29415640 < mcols(reads)[,1] &
                            mcols(reads)[,1] < 30144477 &
                            !is.na(mcols(reads)[,1])),]
    }
    if (length(seqnames(reads))<mates){
        res <- GenomicAlignments::GAlignments()
        return(res)
    }
    clip_reads <- reads[cigar(reads) != "96M",]
    clip_reads <- clip_reads[!grepl("D", cigar(clip_reads)),]
    clip_reads <- clip_reads[!grepl("I", cigar(clip_reads)),]
    if (length(seqnames(clip_reads))<mates){
        res <- GenomicAlignments::GAlignments()
        return(res)
    }
    return(clip_reads)
}

#' Identification of EML4 breakpoint bases
#'
#' This function identifies the basepairs leading up to the EML4 breakpoint.
#'
#' @importFrom BiocBaseUtils isScalarNumber
#' @importFrom S4Vectors isEmpty
#' @param reads `GAlignments` object returned by `EML4_ALK_detection()`.
#' @param basepairs `Integer`, number of basepairs identified
#' from the EML4-ALK fusion. Default=20.
#' @return If EML4-ALK is detected, returns a `table` of identified
#' EML4 basepairs with the number of corresponding reads for each sequence.
#' If no EML4-ALK is detected "No EML4-ALK was detected" is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' EML4_sequence(EML4_ALK_detection(file=H3122_bam,
#'                                     genome="hg38",
#'                                     mates=2),
#'                 basepairs=20)
#' EML4_sequence(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2),
#'                 basepairs=20)
#' @export
EML4_sequence <- function(reads, basepairs=20){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarNumber(basepairs)){
        stop("basepairs has to be a numeric")
    }
    reads <- index_helper(reads)
    EML4_fun <- function(inp){
        return(substring(inp["seq"], (as.numeric(inp["index"]))-(basepairs-1),
                            as.numeric(inp["index"])))
    }
    EML4_seq <- apply(as.data.frame(reads), FUN=EML4_fun, MARGIN=1)
    EML4_tab <- table(EML4_seq)
    return(EML4_tab)
}

#' Identification of ALK breakpoint bases
#'
#' This function identifies the basepairs following the ALK breakpoint.
#'
#' @importFrom BiocBaseUtils isScalarNumber
#' @importFrom S4Vectors isEmpty
#' @param reads `GAlignments` returned by `EML4_ALK_detection()`.
#' @param basepairs `integer`, number of basepairs identified
#' from the EML4-ALK fusion. Default=20.
#' @return If EML4-ALK is detected, returns a `table` of identified
#' ALK basepairs with the number of corresponding reads for each sequence.
#' If no EML4-ALK is detected "No EML4-ALK was detected" is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' ALK_sequence(EML4_ALK_detection(file=H3122_bam,
#'                                 genome="hg38",
#'                                 mates=2),
#'                 basepairs=20)
#' ALK_sequence(EML4_ALK_detection(file=HCC827_bam,
#'                                 genome="hg38",
#'                                 mates=2),
#'                 basepairs=20)
#' @export
ALK_sequence <- function(reads, basepairs=20){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarNumber(basepairs)){
        stop("basepairs has to be a numeric")
    }
    reads <- index_helper(reads)
    ALK_fun <- function(inp){
        return(substring(inp["seq"], (as.numeric(inp["index"])+1),
                            (as.numeric(inp["index"])+basepairs)))
    }
    ALK_seq <- apply(as.data.frame(reads), FUN=ALK_fun, MARGIN=1)
    ALK_tab <- table(ALK_seq)
    return(ALK_tab)
}

#' EML4-ALK breakpoint
#'
#' This function identifies the genomic position in EML4
#' where the breakpoint has happened.
#'
#' @importFrom GenomicAlignments start
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors isEmpty
#' @param reads `GAlignments` object returned by `EML4_ALK_detection()`.
#' @return If EML4-ALK is detected, returns a `table` of genomic positions
#' with the number of corresponding reads for each sequence.
#' If no EML4-ALK is detected "No EML4-ALK was detected" is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' break_position(EML4_ALK_detection(file=H3122_bam,
#'                                     genome="hg38",
#'                                     mates=2))
#' break_position(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2))
#' @export
break_position <- function(reads){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    reads <- index_helper(reads)
    break_pos <- start(reads) + (mcols(reads)["index"][,1]-1)
    break_pos_tab <- table(break_pos)
    return(break_pos_tab)
}


#' Read depth at breakpoint
#'
#' This function identifies the read depth at the basepair
#' before the breakpoint in EML4.
#'
#' @importFrom bamsignals bamCoverage
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors isEmpty
#' @param file The name of the file which the data are to be read from.
#' @param reads `GAlignments` object returned by `EML4_ALK_detection()`.
#' @return If EML4-ALK is detected a single `integer` corresponding
#' to the read depth at the breakpoint is returned.
#' If no EML4-ALK is detected "No EML4-ALK was detected" is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' break_position_depth(file=H3122_bam,
#'                         EML4_ALK_detection(file=H3122_bam,
#'                                             genome="hg38",
#'                                             mates=2))
#' break_position_depth(file=HCC827_bam,
#'                         EML4_ALK_detection(file=HCC827_bam,
#'                                             genome="hg38",
#'                                             mates=2))
#' @export
break_position_depth <- function(file, reads){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    break_pos_tab <- break_position(reads)
    stop_pos <- as.numeric(names(which.max(break_pos_tab)))
    depth <- bamCoverage(file,
                            GRanges(seqnames="chr2",
                            IRanges(start=(stop_pos), end=stop_pos+1)),
                                    mapqual=0, verbose=FALSE)
    return(max(depth[1]))
}
#' Complete EML4-ALK analysis
#'
#' This functions collects the results from the other functions of the package.
#'
#' @importFrom BiocBaseUtils isScalarNumber isScalarCharacter
#' @importFrom S4Vectors isEmpty
#' @param file The name of the file which the data are to be read from.
#' @param genome `character` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param mates `interger`, the minimum number EML4-ALK mate pairs needed
#' to be detected in order to call a variant. Default=2.
#' @param basepairs `integer`, number of basepairs identified
#' from the EML4-ALK fusion. Default=20.
#' @return A `list` object with
#' clipped_reads corresponding to `EML4_ALK_detection()`,
#' last_EML4 corresponding to `EML4_sequence()`,
#' first_ALK corresponding to `ALK_sequence()`,
#' breakpoint corresponding to `break_position()`,
#' and read_depth corresponding to `break_position_depth()`.
#' If no EML4-ALK is detected an empty `GAlignments` is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' EML4_ALK_analysis(file=H3122_bam,
#'                     genome="hg38",
#'                     mates=2,
#'                     basepairs=20)
#' EML4_ALK_analysis(file=HCC827_bam,
#'                     genome="hg38",
#'                     mates=2,
#'                     basepairs=20)
#' @export
EML4_ALK_analysis <- function(file, genome="hg38", mates=2, basepairs=20){
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!isScalarNumber(mates)){
        stop("mates has to be a numeric")
    }
    if (!(genome %in% c("hg38", "hg19"))){
        stop("the reference genome has to be hg38 or hg19")
    }
    res <- EML4_ALK_detection(file=file, genome=genome, mates=mates)
    if(isEmpty(res)){
        return(res)
    }
    if(!isScalarNumber(basepairs)){
        stop("basepairs has to be a numeric")
    }
    EML4 <- EML4_sequence(reads=res, basepairs=basepairs)
    ALK <- ALK_sequence(reads=res, basepairs=basepairs)
    position <- break_position(reads=res)
    position_depth <- break_position_depth(file=file, reads=res)
    return(list(clipped_reads=res,
                last_EML4=EML4,
                first_ALK=ALK,
                breakpoint=position,
                read_depth=position_depth))
}



#' Detection of ALK breakpoint
#' This function identifies the genomic position in ALK
#' where the breakpoint has happened.
#' This function looks for ALK-EML4 mate pair reads in the BAM file.
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom GenomicAlignments readGAlignments cigar seqnames
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom IRanges IRanges
#' @importFrom BiocBaseUtils isScalarNumber isScalarCharacter
#' @param file The name of the file which the data are to be read from.
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param mates `Interger`, the minimum number ALK-EML4 mate pairs
#' needed to be detected in order to call a variant. Default=2.
#' @return A `GAlignments` object with soft-clipped reads representing
#'  ALK-EML4 is returned. If no ALK-EML4 is detected the `GAlignments`
#'  is empty.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' ALK_EML4_detection(file=H3122_bam,
#'                     genome="hg38",
#'                     mates=2)
#' ALK_EML4_detection(file=HCC827_bam,
#'                     genome="hg38",
#'                     mates=2)
#' @export
ALK_EML4_detection <- function(file, genome="hg38", mates=2){ 
  if(!isScalarCharacter(genome)){ 
    stop("genome has to be a character")
  }
  if(!isScalarNumber(mates)){
    stop("mates has to be a numeric")
  }
  if(!(genome %in% c("hg38", "hg19"))){
    stop("The reference genome has to be hg38 or hg19")
  }
  what <- c("mpos", "seq")
  if (genome =="hg38"){
    which <- GRanges(seqnames="chr2",
                     IRanges(start=29192774, end=29921586))
    param <- ScanBamParam(which=which, what=what)
    reads <- readGAlignments(file=file, param=param) #reads a file containing aligned reads
    reads <- reads[(42169353 < mcols(reads)[,1] & 
                      mcols(reads)[,1] < 42332548 &
                      !is.na(mcols(reads)[,1])),]
  }
  else{ #genome=hg19
    which <- GRanges(seqnames="chr2",
                     IRanges(start=29415640, end=30144477))
    param <- ScanBamParam(which=which, what=what)
    reads <- readGAlignments(file=file, param=param)
    reads <- reads[(42396490 < mcols(reads)[,1] &
                      mcols(reads)[,1] < 42559688 &
                      !is.na(mcols(reads)[,1])),]
  }
  if (length(seqnames(reads))<mates){
    res <- GenomicAlignments::GAlignments()
    return(res)
  }
  clip_reads <- reads[cigar(reads) != "96M",] #not equal to 96M
  clip_reads <- clip_reads[!grepl("D", cigar(clip_reads)),]
  clip_reads <- clip_reads[!grepl("I", cigar(clip_reads)),]
  if (length(seqnames(clip_reads))<mates){
    res <- GenomicAlignments::GAlignments()
    return(res)
  }
  return(clip_reads)
}



#' Detect ALK and EML4 introns of the breakpoint
#' This function identifies the introns in ALK and EML4
#' where the breakpoint has happened.
#' @importFrom GenomicFeatures
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom IRanges IRanges

introns_ALK_EML4 <-function(breakpoint_ALK, breakpoint_EML4){
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ALK_txs <- transcriptsBy(txdb_ALK, by="gene")[["238"]]
  ALK_tx_names <- mcols(ALK_txs)$tx_name
  ALK_intron <- suppressWarnings(intronsByTranscript(txdb,use.names=TRUE)[ALK_tx_names])
  EML4_txs <- transcriptsBy(txdb, by="gene")[["27436"]]
  EML4_tx_names <- mcols(EML4_txs)$tx_name
  EML4_intron <- suppressWarnings(intronsByTranscript(txdb,use.names=TRUE)[EML4_tx_names])
  breakpoint_ALK <- as.numeric(names(which.max(breakpoint_ALK)))
  break_gr_ALK <- GRanges(seqnames="chr2",
                          IRanges(start = breakpoint_ALK,end = breakpoint_ALK),
                          strand="*")
  res_ALK <- findOverlaps(break_gr,rev(ALK_intron$ENST00000389048.8),ignore.strand=TRUE)
  intron_ALK <- subjectHits(res_ALK)
  breakpoint_EML4 <- as.numeric(names(which.max(breakpoint_EML4)))
  break_gr_EML4 <- GRanges(seqnames="chr2",
                           IRanges(start = breakpoint_EML4,end = breakpoint_EML4),
                           strand="*")
  res_EML4 <- findOverlaps(break_gr_EML4,EML4_intron$ENST00000318522.10,ignore.strand=TRUE)
  intron_EML4 <- subjectHits(res_EML4)
  df <- data.frame(intron_ALK=intron_ALK, intron_EML4 = intron_EML4)
  return(df)
}



#' Detect the variants of ALK-EML4
#' This function identifies ALK-EML4 variants using the intron of the breakpoint
#' of EML4
#' @importFrom dplyr

find_variants <- function(EML4intron){
  df <- data.frame(Variant=c('Variant 1 (E13,A20)', 'Variant 2 (E20,A20)','Variant 3a/b(E6,A20)', 'variant 4(E15,A20)','variant 5a/b(E2,A20)','variant 5(E18,A20)','variant 7(E14,A20)','variant 8a/b(E17,A20)'),
                   Intron_EML4=c(13,20,6,15,2,18,14,17),
                   Intron_ALK= c(19,19,19,19,19,19,19,19)) 
  return(df %>% filter_all(any_vars(Intron_EML4 %in% c(EML4intron))))
}



























