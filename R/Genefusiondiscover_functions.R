#' Detection of EML4-ALK variants
#'
#' This function looks for EML4-ALK mate pair reads in the BAM file.
#'
#' @param file The name of the file which the data are to be read from.
#' @param genome `Character string` representing the reference genome. Can be either "hg38" or "hg19". Default = "hg38".
#' @param mates `Interger`, the minimum number EML4-ALK mate pairs needed to be detected in order to call a variant. Default = 2.
#' @return If EML4-ALK is detected a `data.frame` with soft-clipped reads representing EML4-ALK is returned. Otherwise "No EML4-ALK was detected" is returned.
#' @examples
#' EML4_ALK_detection(file = H3122_bam,
#'  genome = "hg38",
#'   mates = 2)
#' EML4_ALK_detection(file = HCC827_bam,
#'  genome = "hg38",
#'   mates = 2)
#' @export
EML4_ALK_detection <- function(file, genome = "hg38", mates = 2){
  `%ni%` <- Negate(`%in%`)
  if(!isa(genome, "character")){
    return("ERROR: genome has to be a character")
  }
  if(!isa(mates, "numeric")){
    return("ERROR: mates has to be a numeric")
  }
  if (genome %ni% c("hg38", "hg19")){
    return("ERROR: The reference genome has to be hg38 or hg19")
  }
  what <- c("mpos", "pos", "seq","cigar")
  if (genome =="hg38"){
    which <- GenomicRanges::GRanges(seqnames="chr2", IRanges::IRanges(start = 42169353, end = 42332548))
    param <- Rsamtools::ScanBamParam(which = which, what = what)
    bam <- Rsamtools::scanBam(file = file, param = param)
    reads <- data.frame(sequences = bam$`chr2:42169353-42332548`$seq,
                        mate = bam$`chr2:42169353-42332548`$mpos,
                        position = bam$`chr2:42169353-42332548`$pos,
                        cigar = bam$`chr2:42169353-42332548`$cigar)
    reads <- reads %>% dplyr::filter(mate < 29921586 & mate > 29192774)
  }
  else{
    which <- GenomicRanges::GRanges(seqnames="chr2", IRanges::IRanges(start = 42396490, end = 42559688))
    param <- Rsamtools::ScanBamParam(which = which, what = what)
    bam <- Rsamtools::scanBam(file = file, param = param)
    reads <- data.frame(sequences = bam$`chr2:42396490-42559688`$seq,
                        mate = bam$`chr2:42396490-42559688`$mpos,
                        position = bam$`chr2:42396490-42559688`$pos,
                        cigar = bam$`chr2:42396490-42559688`$cigar)
    reads <- reads %>% dplyr::filter(mate < 30144477 & mate > 29415640)
  }
  if (length(reads$mate)<mates){
    res <- "No EML4-ALK was detected"
    return(res)
  }
  clip_reads <- reads %>% dplyr::filter(cigar != "96M")
  clip_reads <- clip_reads %>% dplyr::filter(!grepl("D",cigar))
  clip_reads <- clip_reads %>% dplyr::filter(!grepl("I",cigar))
  if (length(clip_reads$mate)<mates){
    res <- "No EML4-ALK was detected"
    return(res)
  }
  return(clip_reads)
}

#' Identification of EML4 breakpoint bases
#'
#' This function identifies the basepairs leading up to the EML4 breakpoint.
#'
#' @param reads `Data.frame` returned by EML4_ALK_detection().
#' @param basepairs `Integer`, number of basepairs identified from the EML4-ALK fusion. Default = 20.
#' @return If EML4-ALK is detected, returns a table of identified EML4 basepairs with the number of corresponding reads for each sequence. Otherwise "No EML4-ALK was detected" is returned.
#' @examples
#' EML4_sequence(EML4_ALK_detection(file = H3122_bam,
#'  genome = "hg38",
#'   mates = 2),
#'    basepairs = 20)
#' EML4_sequence(EML4_ALK_detection(file = HCC827_bam,
#'  genome = "hg38",
#'   mates = 2),
#'    basepairs = 20)
#' @export
EML4_sequence <- function(reads, basepairs = 20){
  if(!isa(reads, "data.frame")){
    if(reads == "No EML4-ALK was detected"){
      return("No EML4-ALK was detected")
    }
    else{
      return("ERROR: reads must be a data.frame")
    }
  }
  if(!isa(basepairs, "numeric")){
    return("ERROR: basepairs has to be a numeric")
  }
  fun <- function(str) sub("\\M.*", "",str)
  index1 <- sapply(reads$cigar,FUN = fun)
  fun1 <- function(ind1){
    if (length(ind1)>1){
      return(NA)
    }
    else{
      return(ind1)
    }
  }
  fun <- function(ind){
    splits <- strsplit(ind,split = "S")
    splits <- lapply(splits, as.numeric)
    splits_ind <- lapply(splits,FUN = fun1)
    return(splits_ind)
  }
  index2 <- sapply(index1,FUN = fun)
  reads$indeces <- index2
  reads <- reads %>% dplyr::filter(!is.na(indeces))
  EML4_fun <- function(inp){
    return(substring(inp$sequences,(inp$indeces-(basepairs-1)),inp$indeces))
  }
  char_df <-reads %>% dplyr::select(sequences,indeces)
  EML4_seq <- apply(char_df, FUN = EML4_fun, MARGIN = 1)
  EML4_tab <- table(EML4_seq)
  return(EML4_tab)
}

#' Identification of ALK breakpoint bases
#'
#' This function identifies the basepairs following the ALK breakpoint.
#'
#' @param reads `data.frame` returned by EML4_ALK_detection().
#' @param basepairs `integer`, number of basepairs identified from the EML4-ALK fusion. Default = 20.
#' @return If EML4-ALK is detected, returns a `table` of identified ALK basepairs with the number of corresponding reads for each sequence. Otherwise "No EML4-ALK was detected" is returned.
#' @examples
#' ALK_sequence(EML4_ALK_detection(file = H3122_bam,
#'  genome = "hg38",
#'   mates = 2),
#'    basepairs = 20)
#' ALK_sequence(EML4_ALK_detection(file = HCC827_bam,
#'  genome = "hg38",
#'   mates = 2),
#'    basepairs = 20)
#' @export
ALK_sequence <- function(reads, basepairs = 20){
  if(!isa(reads, "data.frame")){
    if(reads == "No EML4-ALK was detected"){
      return("No EML4-ALK was detected")
    }
    else{
      return("ERROR: reads must be a data.frame")
    }
  }
  if(!isa(basepairs, "numeric")){
    return("ERROR: basepairs has to be a numeric")
  }
  fun <- function(str) sub("\\M.*", "",str)
  index1 <- sapply(reads$cigar,FUN = fun)
  fun1 <- function(ind1){
    if (length(ind1)>1){
      return(NA)
    }
    else{
      return(ind1)
    }
  }
  fun <- function(ind){
    splits <- strsplit(ind,split = "S")
    splits <- lapply(splits, as.numeric)
    splits_ind <- lapply(splits,FUN = fun1)
    return(splits_ind)
  }
  index2 <- sapply(index1,FUN = fun)
  reads$indeces <- index2
  reads <- reads %>% dplyr::filter(!is.na(indeces))
  ALK_fun <- function(inp){
    return(substring(inp$sequences,(inp$indeces+1),(inp$indeces+basepairs)))
  }
  char_df <- reads %>% dplyr::select(sequences,indeces)
  ALK_seq <- apply(char_df, FUN = ALK_fun, MARGIN = 1)
  ALK_tab <- table(ALK_seq)
  return(ALK_tab)
}

#' EML4-ALK breakpoint
#'
#' This function identifies the genomic position in EML4 where the breakpoint has happened.
#'
#' @param reads Data.frame returned by EML4_ALK_detection().
#' @return If EML4-ALK is detected, returns a `table` of genomic positions with the number of corresponding reads for each sequence. Otherwise "No EML4-ALK was detected" is returned.
#' @examples
#' break_position(EML4_ALK_detection(file = H3122_bam,
#'  genome = "hg38",
#'   mates = 2))
#' break_position(EML4_ALK_detection(file = HCC827_bam,
#'  genome = "hg38",
#'   mates = 2))
#' @export
break_position <- function(reads){
  if(!isa(reads, "data.frame")){
    if(reads == "No EML4-ALK was detected"){
      return("No EML4-ALK was detected")
    }
    else{
      return("ERROR: reads must be a data.frame")
    }
  }
  fun <- function(str) sub("\\M.*", "",str)
  index1 <- sapply(reads$cigar,FUN = fun)
  fun1 <- function(ind1){
    if (length(ind1)>1){
      return(NA)
    }
    else{
      return(ind1)
    }
  }
  fun <- function(ind){
    splits <- strsplit(ind,split = "S")
    splits <- lapply(splits, as.numeric)
    splits_ind <- lapply(splits,FUN = fun1)
    return(splits_ind)
  }
  index2 <- sapply(index1,FUN = fun)
  reads$indeces <- index2
  reads <- reads %>% dplyr::filter(!is.na(indeces))
  reads$indeces <- as.numeric(reads$indeces)
  break_pos <- reads$position + reads$indeces-1
  break_pos_tab <- table(break_pos)
  return(break_pos_tab)
}


#' Read depth at breakpoint
#'
#' This function identifies the read depth at the basepair before the breakpoint in EML4.
#'
#' @param file The name of the file which the data are to be read from.
#' @param reads `data.frame` returned by EML4_ALK_detection().
#' @return If EML4-ALK is detected a single integer corresponding to the read depth at the breakpoint is returned. Otherwise "No EML4-ALK was detected" is returned
#' @examples
#' break_position_depth(file = H3122_bam,
#'  EML4_ALK_detection(file = H3122_bam,
#'   genome = "hg38",
#'    mates = 2))
#' break_position_depth(file = HCC827_bam,
#'  EML4_ALK_detection(file = HCC827_bam,
#'   genome = "hg38",
#'    mates = 2))
#' @export
break_position_depth <- function(file, reads){
  if(!isa(reads, "data.frame")){
    if(reads == "No EML4-ALK was detected"){
      return("No EML4-ALK was detected")
    }
    else{
      return("ERROR: reads must be a data.frame")
    }
  }
  fun <- function(str) sub("\\M.*", "",str)
  index1 <- sapply(reads$cigar,FUN = fun)
  fun1 <- function(ind1){
    if (length(ind1)>1){
      return(NA)
    }
    else{
      return(ind1)
    }
  }
  fun <- function(ind){
    splits <- strsplit(ind,split = "S")
    splits <- lapply(splits, as.numeric)
    splits_ind <- lapply(splits,FUN = fun1)
    return(splits_ind)
  }
  index2 <- sapply(index1,FUN = fun)
  reads$indeces <- index2
  reads <- reads %>% dplyr::filter(!is.na(indeces))
  reads$indeces <- as.numeric(reads$indeces)
  break_pos <- reads$position + reads$indeces-1
  break_pos_tab <- table(break_pos)
  stop_pos <- as.numeric(names(which.max(break_pos_tab)))
  depth <- bamsignals::bamCoverage(file,
                                   GenomicRanges::GRanges(seqnames = "chr2",IRanges::IRanges(start=(stop_pos),end=stop_pos+1)),
                                   mapqual=0,verbose=F)
  return(max(depth[1]))
}

#' Complete EML4-ALK analysis
#'
#' This functions collects the results from the other functions of the package.
#'
#' @param file The name of the file which the data are to be read from.
#' @param genome `character` representing the reference genome. Can be either "hg38" or "hg19". Default = "hg38".
#' @param mates `interger`, the minimum number EML4-ALK mate pairs needed to be detected in order to call a variant. Default = 2.
#' @param basepairs `integer`, number of basepairs identified from the EML4-ALK fusion. Default = 20.
#' @return A `list` object with clipped_reads corresponding to `EML4_ALK_detection()`, last_EML4 corresponding to `EML4_sequence()`, first_ALK corresponding to `ALK_sequence()`, breakpoint corresponding to `break_position()`, and read_depth corresponding to `break_position_depth()`.
#' @examples
#' EML4_ALK_analysis(file = H3122_bam,
#'  genome = "hg38",
#'   mates = 2,
#'    basepairs = 20)
#' EML4_ALK_analysis(file = HCC827_bam,
#'  genome = "hg38",
#'   mates = 2,
#'    basepairs = 20)
#' @export
EML4_ALK_analysis <- function(file, genome = "hg38", mates = 2, basepairs = 20){
  if(!isa(genome, "character")){
    return("ERROR: genome has to be a character")
  }
  if(!isa(mates, "numeric")){
    return("ERROR: mates has to be a numeric")
  }
  if (genome %ni% c("hg38", "hg19")){
    return("ERROR: The reference genome has to be hg38 or hg19")
  }
  res <- EML4_ALK_detection(file = file, genome = genome, mates = mates)
  if(!isa(res, "data.frame")){
    return(res)
  }
  if(!isa(basepairs, "numeric")){
    return("ERROR: basepairs has to be a numeric")
  }
  EML4 <- EML4_sequence(reads = res, basepairs = basepairs)
  ALK <- ALK_sequence(reads = res, basepairs = basepairs)
  position <- break_position(reads = res)
  position_depth <- break_position_depth(file = file, reads = res)
  return(list(clipped_reads = res,
              last_EML4 = EML4,
              first_ALK = ALK,
              breakpoint = position,
              read_depth = position_depth))
}

