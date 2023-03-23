testthat::context("Test of DNAfusion functions")
library(DNAfusion)

file <- system.file("extdata","HCC827_EML4.bam", package = "DNAfusion")
gene <-  "EML4"
detection <- EML4_ALK_detection(file = file)
EML4_seq <- EML4_sequence(EML4_ALK_detection(file = file))
ALK_seq <- ALK_sequence(EML4_ALK_detection(file = file))
pos <- break_position(EML4_ALK_detection(file = file), gene=gene)
pos_depth <- break_position_depth(file = file, EML4_ALK_detection(file = file),
                                  gene=gene)
analysis <- EML4_ALK_analysis(file = file)

testthat::expect_true(is.character(file))
testthat::expect_true(is.character(gene))
testthat::expect_true(isa(detection, "GAlignments"))
testthat::expect_true(is.table(EML4_seq) | is.character(EML4_seq))
testthat::expect_true(is.table(ALK_seq) | is.character(ALK_seq))
testthat::expect_true(is.table(pos) | is.character(pos))
testthat::expect_true(is.integer(pos_depth) | is.character(pos_depth))
testthat::expect_true(is.list(analysis) | isa(analysis, "GAlignments"))


file <- system.file("extdata","H3122_EML4.bam", package = "DNAfusion")
file <- file.path("C:/Users/emma-/OneDrive - Aarhus Universitet/6. Semester/Bachelorprojekt/Rstudio/ALK og raske/Sorted_ALK1.bam")
gene <-  "EML4"
detection <- EML4_ALK_detection(file = file)
EML4_seq <- EML4_sequence(EML4_ALK_detection(file = file))
ALK_seq <- ALK_sequence(EML4_ALK_detection(file = file))
pos <- break_position(EML4_ALK_detection(file = file), gene=gene)
pos_depth <- break_position_depth(file = file, EML4_ALK_detection(file = file),
                                  gene=gene)
analysis <- EML4_ALK_analysis(file = file)

testthat::expect_true(is.character(file))
testthat::expect_true(is.character(gene))
testthat::expect_true(isa(detection, "GAlignments"))
testthat::expect_true(is.table(EML4_seq) | is.character(EML4_seq))
testthat::expect_true(is.table(ALK_seq) | is.character(ALK_seq))
testthat::expect_true(is.table(pos) | is.character(pos))
testthat::expect_true(is.integer(pos_depth) | is.character(pos_depth))
testthat::expect_true(is.list(analysis) | isa(analysis, "GAlignments"))
testthat::expect_error(break_position(detection,gene = "character"))
testthat::expect_error(break_position(detection,gene = 1))
testthat::expect_error(break_position(detection,genome = "38",gene=1))
testthat::expect_error(break_position(detection,genome = "38",gene="character"))
testthat::expect_error(break_position(detection,genome = 38,gene=1))
testthat::expect_error(break_position(detection,genome = 38,gene="character"))
testthat::expect_error(break_position_depth(file=file,detection,gene=1))
testthat::expect_error(break_position_depth(file=file, detection,
                                            gene="character"))
testthat::expect_error(break_position_depth(file=file, detection,genome = 38))
testthat::expect_error(break_position_depth(file=file, detection,genome = "38"))                       
testthat::expect_error(ALK_sequence(detection, basepairs = "20"))
testthat::expect_error(EML4_sequence(detection, basepairs = "20"))

testthat::expect_error(EML4_ALK_detection(file = file, genome = 38))
testthat::expect_error(EML4_ALK_detection(file = file, mates = "38"))
testthat::expect_error(EML4_ALK_detection(file = file, genome = "hg37"))


testthat::expect_error(ALK_sequence(detection,genome="38"))
testthat::expect_error(ALK_sequence(detection,genome=38))
testthat::expect_error(EML4_sequence(detection,genome="38"))
testthat::expect_error(EML4_sequence(detection,genome=38))

detection <- EML4_ALK_detection(file = file, genome = "hg19")
testthat::expect_true(isa(detection, "GAlignments"))

testthat::expect_error(EML4_sequence(EML4_ALK_detection(file=file),
                                     genome = "hg19"))
testthat::expect_error(ALK_sequence(EML4_ALK_detection(file=file),
                                     genome = "hg19"))
testthat::expect_error(break_position(EML4_ALK_detection(file=file),gene="EML4",
                                      genome = "hg19"))
testthat::expect_error(break_position(EML4_ALK_detection(file=file),gene="ALK",
                                      genome = "hg19"))


testthat::expect_error(EML4_sequence(file))
testthat::expect_error(ALK_sequence(file))
testthat::expect_error(break_position(file))
testthat::expect_error(break_position_depth(reads = file))
testthat::expect_error(EML4_ALK_analysis(file, genome = 38))
testthat::expect_error(EML4_ALK_analysis(file, mates = "38"))
testthat::expect_error(EML4_ALK_analysis(file, genome = "37"))
testthat::expect_error(EML4_ALK_analysis(file, basepairs = "20"))




#EMMAS KODE
file <- system.file("extdata","HCC827_EML4.bam", package = "DNAfusion")
introns <- introns_ALK_EML4(file = file)
variants <- find_variants(file = file)

testthat::expect_true(is.character(file))
testthat::expect_true(is.data.frame(introns) | is.character(introns))
testthat::expect_true(is.data.frame(variants) | is.character(variants))

testthat::expect_error(introns_ALK_EML4(file = file, genome = 38))
testthat::expect_error(introns_ALK_EML4(file = file, genome = "hg37"))
testthat::expect_error(find_variants(file = file, genome = 38))
testthat::expect_error(find_variants(file = file, genome = "hg37"))

file <- system.file("extdata","H3122_EML4.bam", package = "DNAfusion")
introns <- introns_ALK_EML4(file = file)
variants <- find_variants(file = file)

testthat::expect_true(is.character(file))
testthat::expect_true(is.data.frame(introns) | is.character(introns))
testthat::expect_true(is.data.frame(variants) | is.character(variants))







