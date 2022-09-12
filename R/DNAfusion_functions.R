#' @keywords internal
#' @importFrom GenomicAlignments cigar
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
    mcols(input)$indeces <- unlist(index2)
    input <- input[!is.na(mcols(input[,3]))[,1],]
    return(input)
}
#' Detection of EML4-ALK variants
#' 
#' This function looks for EML4-ALK mate pair reads in the BAM file.
#' @import dplyr
#' @importFrom GenomicRanges GRanges 
#' @importFrom GenomicAlignments readGAlignments cigar seqnames
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom IRanges IRanges
#' @param file The name of the file which the data are to be read from.
#' @param genome `Character string` representing the reference genome. 
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param mates `Interger`, the minimum number EML4-ALK mate pairs
#' needed to be detected in order to call a variant. Default=2.
#' @return If EML4-ALK is detected a `GAlignments` object with soft-clipped
#' reads representing EML4-ALK is returned. 
#' Otherwise "No EML4-ALK was detected" is returned.
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
    if(!isa(genome, "character")){
        stop("genome has to be a character")
    }
    if(!isa(mates, "numeric")){
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
        res <- "No EML4-ALK was detected"
        return(res)
    }
    clip_reads <- reads[cigar(reads) != "96M",]
    clip_reads <- clip_reads[!grepl("D", cigar(clip_reads)),]
    clip_reads <- clip_reads[!grepl("I", cigar(clip_reads)),]
    if (length(seqnames(reads))<mates){
        res <- "No EML4-ALK was detected"
        return(res)
    }
    return(clip_reads)
}

#' Identification of EML4 breakpoint bases
#'
#' This function identifies the basepairs leading up to the EML4 breakpoint.
#'
#' @import dplyr
#' @param reads `GAlignments` object returned by EML4_ALK_detection().
#' @param basepairs `Integer`, number of basepairs identified 
#' from the EML4-ALK fusion. Default=20.
#' @return If EML4-ALK is detected, returns a `table` of identified
#' EML4 basepairs with the number of corresponding reads for each sequence.
#' Otherwise "No EML4-ALK was detected" is returned.
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
        if(reads == "No EML4-ALK was detected"){
            return("No EML4-ALK was detected")
        }
        else{
            stop("reads must be a GAlignments object")
        }
    }
    if(!isa(basepairs, "numeric")){
        stop("basepairs has to be a numeric")
    }
    reads <- index_helper(reads)
    EML4_fun <- function(inp){
        return(substring(inp[10], (as.numeric(inp[11]))-(basepairs-1),
                         as.numeric(inp[11])))
    }
    EML4_seq <- apply(as.data.frame(reads), FUN=EML4_fun, MARGIN=1)
    EML4_tab <- table(EML4_seq)
    return(EML4_tab)
}

#' Identification of ALK breakpoint bases
#'
#' This function identifies the basepairs following the ALK breakpoint.
#'
#' @import dplyr
#' @param reads `GAlignments` returned by EML4_ALK_detection().
#' @param basepairs `integer`, number of basepairs identified
#' from the EML4-ALK fusion. Default=20.
#' @return If EML4-ALK is detected, returns a `table` of identified
#' ALK basepairs with the number of corresponding reads for each sequence.
#' Otherwise "No EML4-ALK was detected" is returned.
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
        if(reads == "No EML4-ALK was detected"){
            return("No EML4-ALK was detected")
        }
        else{
            stop("reads must be a GAlignments object")
        }
    }
    if(!isa(basepairs, "numeric")){
        stop("basepairs has to be a numeric")
    }
    reads <- index_helper(reads)
    ALK_fun <- function(inp){
        return(substring(inp[10], (as.numeric(inp[11])+1),
                            (as.numeric(inp[11])+basepairs)))
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
#' @import dplyr
#' @importFrom GenomicAlignments start 
#' @param reads `GAlignments` object returned by EML4_ALK_detection().
#' @return If EML4-ALK is detected, returns a `table` of genomic positions
#' with the number of corresponding reads for each sequence. 
#' Otherwise "No EML4-ALK was detected" is returned.
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
        if(reads == "No EML4-ALK was detected"){
            return("No EML4-ALK was detected")
        }
        else{
            stop("reads must be a GAlignments object")
        }
    }
    reads <- index_helper(reads)
    break_pos <- start(reads) + (mcols(reads)[,3]-1)
    break_pos_tab <- table(break_pos)
    return(break_pos_tab)    
}


#' Read depth at breakpoint
#'
#' This function identifies the read depth at the basepair
#' before the breakpoint in EML4.
#'
#' @import dplyr
#' @importFrom bamsignals bamCoverage
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @param file The name of the file which the data are to be read from.
#' @param reads `GAlignments` object returned by EML4_ALK_detection().
#' @return If EML4-ALK is detected a single integer corresponding
#' to the read depth at the breakpoint is returned. 
#' Otherwise "No EML4-ALK was detected" is returned
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
        if(reads == "No EML4-ALK was detected"){
            return("No EML4-ALK was detected")
            }
        else{
            stop("reads must be a GAlignments object")
            }
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
#' @import dplyr
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
    if(!isa(genome, "character")){
        return("ERROR: genome has to be a character")
    }
    if(!isa(mates, "numeric")){
        return("ERROR: mates has to be a numeric")
    }
    if (!(genome %in% c("hg38", "hg19"))){
        return("ERROR: The reference genome has to be hg38 or hg19")
    }
    res <- EML4_ALK_detection(file=file, genome=genome, mates=mates)
    if(!isa(res, "GAlignments")){
        return(res)
    }
    if(!isa(basepairs, "numeric")){
        return("ERROR: basepairs has to be a numeric")
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

