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
#' Detection of ALK and EML4 breakpoint
#' 
#' This function identifies the genomic position in ALK and EML4
#' where the breakpoint has happened.
#' This function looks for ALK-EML4 and EML4-ALK mate pair reads 
#' in the BAM file.
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom GenomicAlignments readGAlignments cigar seqnames
#' @importFrom Rsamtools ScanBamParam scanBam
#' @importFrom IRanges IRanges
#' @importFrom BiocBaseUtils isScalarNumber isScalarCharacter
#' @importFrom BiocGenerics append
#' @param file The name of the file which the data are to be read from.
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param mates `Interger`, the minimum number ALK-EML4 mate pairs
#' needed to be detected in order to call a variant. Default=2.
#' @return A `GAlignments` object with soft-clipped reads representing
#'  ALK-EML4/EML4_ALK is returned. If no ALK-EML4 or EML4-ALK is detected the 
#'  `GAlignments`is empty.
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
    clip_reads_EML4 <- reads[cigar(reads) != "96M",]
    clip_reads_EML4 <- clip_reads_EML4[!grepl("D", cigar(clip_reads_EML4)),]
    clip_reads_EML4 <- clip_reads_EML4[!grepl("I", cigar(clip_reads_EML4)),]
    if (length(seqnames(clip_reads_EML4))<mates){
        res <- GenomicAlignments::GAlignments()
        return(res)
    }
    what <- c("mpos", "seq")
    if (genome =="hg38"){
        which <- GRanges(seqnames="chr2",
                            IRanges(start=29192774, end=29921586))
        param <- ScanBamParam(which=which, what=what)
        reads <- readGAlignments(file=file, param=param)
        reads <- reads[(42169353 < mcols(reads)[,1] & 
                            mcols(reads)[,1] < 42332548 &
                            !is.na(mcols(reads)[,1])),]
    }
    else{ 
        which <- GRanges(seqnames="chr2",
                            IRanges(start=29415640, end=30144477))
        param <- ScanBamParam(which=which, what=what)
        reads <- readGAlignments(file=file, param=param)
        reads <- reads[(42396490 < mcols(reads)[,1] &
                            mcols(reads)[,1] < 42559688 &
                            !is.na(mcols(reads)[,1])),]
    }
    clip_reads_ALK <- reads[cigar(reads) != "96M",]
    clip_reads_ALK <- clip_reads_ALK[!grepl("D", cigar(clip_reads_ALK)),]
    clip_reads_ALK <- clip_reads_ALK[!grepl("I", cigar(clip_reads_ALK)),]
    return(append(clip_reads_EML4,clip_reads_ALK))
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
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
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
#'                 basepairs=20,
#'                 genome="hg38")
#' EML4_sequence(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2),
#'                 basepairs=20,
#'                 genome="hg38")
#' @export
EML4_sequence <- function(reads,basepairs=20,genome="hg38"){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarNumber(basepairs)){
        stop("basepairs has to be a numeric")
    }
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    if (genome =="hg38"){
        reads <- reads[start(reads) > 42169353 & start(reads) < 42332548]
    }
    else{
        reads <- reads[start(reads) > 42396490 & start(reads) < 42559688]
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
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @return If EML4-ALK is detected, returns a `table` of identified
#' ALK basepairs with the number of corresponding reads for each sequence.
#' If no spanning reads in ALK is detected an empty GAlignments 
#' object is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#'
#' ALK_sequence(EML4_ALK_detection(file=H3122_bam,
#'                                     genome="hg38",
#'                                     mates=2),
#'                 basepairs=20,
#'                 genome="hg38")
#' ALK_sequence(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2),
#'                 basepairs=20,
#'                 genome="hg38")
#' @export
ALK_sequence <- function(reads, basepairs=20,genome="hg38"){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarNumber(basepairs)){
        stop("basepairs has to be a numeric")
    }
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    if (genome =="hg38"){
        reads <- reads[start(reads) > 29192774 & start(reads) < 29921586]
    }
    else{
        reads <- reads[start(reads) > 29415640 & start(reads) < 30144477]
    }
    if(isEmpty(reads)){
        return(reads)
    }
    reads <- index_helper(reads)
    ALK_fun <- function(inp){
        return(substring(inp["seq"], (as.numeric(inp["index"]))-(basepairs-1),
                            as.numeric(inp["index"])))
    }
    ALK_seq <- apply(as.data.frame(reads), FUN=ALK_fun, MARGIN=1)
    ALK_tab <- table(ALK_seq)
    return(ALK_tab)
}

#' EML4-ALK breakpoint
#'
#' This function identifies the genomic position in EML4 or ALK,
#' where the breakpoint has happened.
#'
#' @importFrom GenomicAlignments start
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors isEmpty
#' @param reads `GAlignments` object returned by `EML4_ALK_detection()`.
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param gene `Character string` representing the gene.
#' Can be either "ALK" or "EML4". 
#' @return If EML4-ALK is detected, it returns a `table` of genomic positions
#' with the number of corresponding reads for each sequence.
#' If no spanning reads in EML4 or ALK is detected
#' an empty GAlignments object is returned.
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
#'                                     mates=2),gene="EML4",genome="hg38")
#' break_position(EML4_ALK_detection(file=H3122_bam,
#'                                     genome="hg38",
#'                                     mates=2),gene="ALK",genome="hg38")
#' break_position(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2),gene="EML4",genome="hg38")
#' break_position(EML4_ALK_detection(file=HCC827_bam,
#'                                     genome="hg38",
#'                                     mates=2),gene="ALK",genome="hg38")
#' @export
break_position <- function(reads,gene,genome="hg38"){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    if(!isScalarCharacter(gene)){
        stop("gene has to be a character")
    }
    if(!(gene %in% c("ALK", "EML4"))){
        stop("The gene has to be ALK or EML4")
    }
    if(gene == "EML4"){
        if (genome =="hg38"){
            reads <- reads[start(reads) > 42169353 & start(reads) < 42332548]
        }
        else{
            reads <- reads[start(reads) > 42396490 & start(reads) < 42559688]
        }
    }
    if (gene == "ALK"){
        if (genome =="hg38"){
            reads <- reads[start(reads) > 29192774 & start(reads) < 29921586]
        }
        else{
            reads <- reads[start(reads) > 29415640 & start(reads) < 30144477]
        }
    }
    if(isEmpty(reads)){
        return(reads)
    }
    reads <- index_helper(reads)
    break_pos <- start(reads) + (mcols(reads)["index"][,1]-1)
    break_pos_tab <- table(break_pos)
    return(break_pos_tab)
}


#' Read depth at breakpoint
#'
#' This function identifies the read depth at the basepair
#' before the breakpoint in EML4 or ALK
#'
#' @importFrom bamsignals bamCoverage
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors isEmpty
#' @param file The name of the file which the data are to be read from.
#' @param reads `GAlignments` object returned by `EML4_ALK_detection()`.
#' @param genome `Character string` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @param gene `Character string` representing the gene.
#' Can be either "ALK" or "EML4".
#' @return If EML4-ALK is detected a single `integer` corresponding
#' to the read depth at the breakpoint is returned.
#' If no spanning reads in EML4 or ALK is detected
#' an empty GAlignments object is returned.
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
#'                                             mates=2),
#'                                             gene="ALK",genome="hg38")
#' break_position_depth(file=H3122_bam,
#'                         EML4_ALK_detection(file=H3122_bam,
#'                         genome="hg38",
#'                         mates=2),
#'                         gene="EML4",genome="hg38")              
#' break_position_depth(file=HCC827_bam,
#'                         EML4_ALK_detection(file=HCC827_bam,
#'                                             genome="hg38",
#'                                             mates=2),
#'                                             gene="ALK",genome="hg38")
#' break_position_depth(file=H3122_bam,
#'                         EML4_ALK_detection(file=H3122_bam,
#'                                             genome="hg38",
#'                                             mates=2),
#'                                             gene="EML4",genome="hg38")
#' @export
break_position_depth <- function(file,reads,gene,genome="hg38"){
    if(!isa(reads, "GAlignments")){
        stop("reads must be a GAlignments object")
    }
    if(isEmpty(reads)){
        return("No EML4-ALK was detected")
    }
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    if(!isScalarCharacter(gene)){
        stop("gene has to be a character")
    }
    if(!(gene %in% c("ALK", "EML4"))){
        stop("The gene has to be ALK or EML4")
    }
    break_pos_tab <- break_position(reads,gene,genome)
    if(isEmpty(break_pos_tab)){
        return(break_pos_tab)
    }
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
#' breakpoint_ALK corresponding to `break_position()`, gene = "ALK",
#' breakpoint_EML4 corresponding to `break_position()`,gene = "EML4",
#' read_depth_ALK corresponding to `break_position_depth()`.gene = "ALK",
#' and read_depth_EML4 corresponding to `break_position_depth()`
#' ,gene = "EML4".
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
    EML4 <- EML4_sequence(reads=res, basepairs=basepairs,genome)
    ALK <- ALK_sequence(reads=res, basepairs=basepairs,genome)
    if(isEmpty(ALK)){
        ALK <- NA
    }
    position_ALK <- break_position(reads=res,gene="ALK",genome)
    if(isEmpty(position_ALK)){
        position_ALK <- NA
    }
    position_EML4 <- break_position(reads=res,gene="EML4",genome)
    position_depth_ALK <- break_position_depth(file=file, reads=res,gene="ALK"
                                                ,genome)
    if(isEmpty(position_depth_ALK)){
        position_depth_ALK <- NA
    }
    position_depth_EML4 <- break_position_depth(file=file, reads=res,gene="EML4"
                                                ,genome)
    return(list(clipped_reads=res,
                last_EML4=EML4,
                first_ALK=ALK,
                breakpoint_ALK=position_ALK,
                breakpoint_EML4=position_EML4,
                read_depth_ALK=position_depth_ALK,
                read_depth_EML4=position_depth_EML4))
}



#' Detect ALK and EML4 introns of the breakpoint
#' 
#' This function identifies the introns of ALK and EML4
#' where the breakpoint has happened.
#' 
#' @importFrom GenomicFeatures intronsByTranscript transcriptsBy
#' @importFrom S4Vectors subjectHits
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene 
#' TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom GenomicRanges GRanges mcols 
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom BiocBaseUtils isScalarCharacter
#' @param file The name of the file which the data are to be read from.
#' @param genome `character` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @return A`dataframe`of the ALK- and EML4-intron of the breakpoint is returned
#' corresponding to the transcript ENST00000389048.8 for ALK and 
#' ENST00000318522.10 for EML4.
#' If the breakpoint is not located in introns of ALK or EML4, 
#' "Breakpoint not located in intron of ALK" or 
#' "Breakpoint not located in intron of EML4" is returned.
#' @examples
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#' introns_ALK_EML4(file=H3122_bam,genome="hg38")
#' introns_ALK_EML4(file=HCC827_bam,genome="hg38")
#' @export
introns_ALK_EML4 <-function(file, genome="hg38"){
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    ALK_txs <- transcriptsBy(txdb, by="gene")[["238"]]
    ALK_tx_names <- mcols(ALK_txs)$tx_name
    ALK_intron <- suppressWarnings(intronsByTranscript(txdb,use.names=TRUE)
                                    [ALK_tx_names])
    EML4_txs <- transcriptsBy(txdb, by="gene")[["27436"]]
    EML4_tx_names <- mcols(EML4_txs)$tx_name
    EML4_intron <- suppressWarnings(intronsByTranscript(txdb, use.names=TRUE)
                                    [EML4_tx_names])
    breakpoint_ALK <-break_position(EML4_ALK_detection(file=file,genome=genome),
                                    gene="ALK",genome=genome)
    breakpoint_EML4 <-break_position(EML4_ALK_detection(file=file,
                                    genome=genome),
                                    gene="EML4",genome=genome)
    if (isScalarCharacter(breakpoint_ALK)){
        return("No ALK-EML4 was detected")
    }
    breakpoint_ALK <- as.numeric(names(which.max(breakpoint_ALK)))
    break_gr_ALK <- GRanges(seqnames="chr2",
                            IRanges(start = breakpoint_ALK,
                                    end = breakpoint_ALK),strand="*")
    res_ALK <- findOverlaps(break_gr_ALK,rev(ALK_intron$ENST00000389048.8),
                            ignore.strand=TRUE)
    intron_ALK <- subjectHits(res_ALK)
    if (length(intron_ALK)==0){
        return("Breakpoint not located in intron of ALK")
    }
    breakpoint_EML4 <- as.numeric(names(which.max(breakpoint_EML4)))
    break_gr_EML4 <- GRanges(seqnames="chr2",IRanges(start = breakpoint_EML4,
                            end = breakpoint_EML4),
                            strand="*")
    res_EML4 <- findOverlaps(break_gr_EML4,EML4_intron$ENST00000318522.10,
                            ignore.strand=TRUE)
    intron_EML4 <- subjectHits(res_EML4)
    if (length(intron_EML4)==0){
        return("Breakpoint not located in intron of EML4")
    }
    df <- data.frame(intron_ALK=intron_ALK, intron_EML4 = intron_EML4)
    return(df)
}



#' Detect the variants of ALK-EML4
#' 
#' This function identifies ALK-EML4 variants using the intron of the breakpoint
#' of EML4
#' 
#' @importFrom S4Vectors isEmpty
#' @param file The name of the file which the data are to be read from.
#' @param genome `character` representing the reference genome.
#' Can be either "hg38" or "hg19". Default="hg38".
#' @return A `dataframe`of the ALK-EML4 variant is returned.
#' If no variant is detected, "No ALK-EML4 was detected" is returned.
#' @examples 
#' H3122_bam <- system.file("extdata",
#' "H3122_EML4.bam",
#' package="DNAfusion")
#' HCC827_bam <-  system.file("extdata",
#' "HCC827_EML4.bam",
#' package="DNAfusion")
#' find_variants(file=H3122_bam,genome="hg38")
#' find_variants(file=HCC827_bam,genome="hg38")
#' @export
find_variants <- function(file, genome="hg38"){
    if(!isScalarCharacter(genome)){
        stop("genome has to be a character")
    }
    if(!(genome %in% c("hg38", "hg19"))){
        stop("The reference genome has to be hg38 or hg19")
    }
    introns <-introns_ALK_EML4(file=file,genome=genome)
    if (isScalarCharacter(introns)){
        if(isEmpty(EML4_ALK_detection(file,genome = genome))){
            return("No ALK-EML4 was detected")
        }
        else{
            return("Breakpoint not located in intron of ALK or EML4")
        }
    }
    EML4intron <-introns_ALK_EML4(file=file,genome=genome)$intron_EML4
    ALKintron <-introns_ALK_EML4(file=file,genome=genome)$intron_ALK
    df <- data.frame(Variant=c('Variant 1 (E13,A20)', 'Variant 2 (E20,A20)',
                                'Variant 3a/b(E6,A20)', 'variant 4(E15,A20)',
                                'variant 5a/b(E2,A20)','variant 5(E18,A20)',
                                'variant 7(E14,A20)','variant 8a/b(E17,A20)'),
                                Intron_EML4=c(13,20,6,15,2,18,14,17),
                                Intron_ALK= c(19,19,19,19,19,19,19,19))
    if (!(EML4intron %in% c(13,20,6,15,2,18,14,17))){
        return(list("Variant of breakpoint-introns not classified",
                    EML4_intron = EML4intron,
                    ALK_intron = ALKintron))
    }
    res <- df[df[,2] %in% EML4intron,]
    return(res)
}





