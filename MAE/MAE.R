suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(tMAE))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(MafDb.gnomAD.r2.1.hs37d5))

# take arguments (dataDir)
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[[1]]

# extract bam files
mae.all <- list.files(path = dataDir,
                      pattern = "*.vcf$",
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)

mafdb <- MafDb.gnomAD.r2.1.hs37d5 
# mafdb <- MafDb.gnomAD.r2.1.GRCh38
dir.create("MAE_results", showWarnings = FALSE)

for (i in c(1:length(mae.all))){
  sampleID <- str_match(mae.all[i], "(.*?).sorted.dedup.SplitNCigar.recal.bam.ASEcalling.vcf")[2]

  # allelic counts generated using ASEReadCounter
  allelicCountsFile <- paste0(dataDir, '/', mae.all[i])
  allelicCounts <- fread(allelicCountsFile)
  
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
  write.csv(resMAE, paste0('MAE_results/', sampleID, '.MAE.result.csv'), row.names=FALSE)
}
