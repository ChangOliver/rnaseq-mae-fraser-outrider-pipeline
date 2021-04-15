# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input vcf file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character"),
  make_option(c("-d", "--dataset"), type="numeric", default=NULL, 
              help="37 for hs37d5 and 38 for GRCh38", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Input & output must be supplied.", call.=FALSE)
}

# define dataset
if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("Specify a dataset to use.", call.=FALSE)
} else if (opt$dataset == 37){
  library(MafDb.gnomAD.r2.1.hs37d5)
  mafdb <- MafDb.gnomAD.r2.1.hs37d5 
} else if (opt$dataset == 38){
  library(MafDb.gnomAD.r2.1.GRCh38)
  mafdb <- MafDb.gnomAD.r2.1.GRCh38
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))

# main --------------------------------------------------------------------
library(tMAE)
library(dplyr)
library(stringr) # str_match

# extract sample name
sampleID <- str_match(basename(opt$input), "(.*?).vcf")[2] # sampleID.sorted.dedup.SplitNCigar.recal.bam.ASEcalling

# allelic counts generated using ASEReadCounter
allelicCounts <- fread(opt$input)

print(paste0("Processing ", sampleID))
resMAE <- DESeq4MAE(allelicCounts, minCoverage = 10)
resMAE[, signif := padj < 0.05]
resMAE[, signif_ALT := signif == TRUE & altRatio >= 0.8]

# convert results table into GRanges object
rng <- GRanges(seqnames = resMAE$contig, 
               ranges = IRanges(start = resMAE$position, end = resMAE$position), 
               strand = '*')
resMAE$gnomadAF <- gscores(mafdb, rng)$AF
resMAE[, rare := (gnomadAF <= 0.01 | is.na(gnomadAF))]
resMAE <- mutate(resMAE, sampleID = sampleID) %>%
          relocate(sampleID, .before=contig) %>% 
          replace(is.na(.), ".")
write.csv(resMAE, paste0(opt$output, sampleID, '.MAE.result.csv'), row.names=FALSE)
