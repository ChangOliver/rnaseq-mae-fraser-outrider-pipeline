# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input vcf directory path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory path", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Input & output must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$input)){
  print_help(opt_parser)
  stop("Input must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$input <- ifelse(substr(opt$input, nchar(opt$input), nchar(opt$input))=='/', opt$input, paste0(opt$input,'/'))
opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))

# main --------------------------------------------------------------------
library(ggplot2)
library(data.table)
library(tMAE)
library(dplyr)
library(stringr)
library(MafDb.gnomAD.r2.1.hs37d5)
# library(MafDb.gnomAD.r2.1.GRCh38)

# extract vcf files
mae.all <- list.files(path = opt$input, pattern="vcf$")

# define dataset
# mafdb <- MafDb.gnomAD.r2.1.GRCh38
mafdb <- MafDb.gnomAD.r2.1.hs37d5 

# MAE analysis per sample
for (i in c(1:length(mae.all))){
  # extract sample name
  sampleID <- str_match(mae.all[i], "(.*?).sorted.dedup.SplitNCigar.recal.bam.ASEcalling.vcf")[2]

  # allelic counts generated using ASEReadCounter
  allelicCountsFile <- paste0(opt$input, mae.all[i])
  allelicCounts <- fread(allelicCountsFile)
  
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
}
