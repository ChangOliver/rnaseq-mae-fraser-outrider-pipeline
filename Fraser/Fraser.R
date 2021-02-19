suppressMessages(library(FRASER))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# take arguments (dataDir)
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[[1]]

# extract bam files
bam.all <- list.files(path = dataDir,
                      pattern = "*.bam$",
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)

# create file table
dataInfo <- data.table(sampleID = str_match(bam.all, "(.*?).sorted.dedup.bam")[,2], bamFile = paste0(dataDir, '/', bam.all), pairedEnd=TRUE)

settings <- FraserDataSet(colData=dataInfo, workingDir='all_wkdir')
register(MulticoreParam(workers=min(10, multicoreWorkers())))
fds <- countRNAData(settings)

# compute stats
fds <- calculatePSIValues(fds)

# filtering junction with low expression
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20,
                                      minDeltaPsi=0.0, filter=TRUE)

# use PCA to speed up the tutorial
fds <- FRASER(fds, q=bestQ(fds, type="psi5"), implementation="PCA")

# alternatively, we also provide a way to use biomart for the annotation:
fds <- annotateRanges(fds)
# get results: we recommend to use an FDR cutoff 0.05
res <- results(fds, zScoreCutoff=NA, padjCutoff=0.05, deltaPsiCutoff=0.3)

samples <- unique(res$sampleID)
dir.create("Fraser_results", showWarnings = FALSE)
for (s in c(1:length(samples))){
  res_s <- res %>% filter(grepl(samples[s], sampleID)) %>% replace(is.na(.), ".")
  write.csv(res_s, paste0("Fraser_results/", samples[s], '.FRASER.result.csv'), row.names=FALSE)
}

saveFraserDataSet(fds, dir='all_wkdir', name="Data_Analysis")

# images
samples = bam.all
types = c("theta", "psi5", "psi3")

for (sample in samples){
  path = paste0("./results/", sample, "/Volcano/")
  dir.create(path, showWarnings = FALSE)
  for (type in types){
    png(file=paste0(path, sample, "-Volcano-", type, ".png"))
    print(plotVolcano(fds, sampleID=sample, type=type))
    dev.off()
  }
}

for (sample in samples){
  extract <- res[res$sampleID == sample]
  path = paste0("./results/", sample, "/QQ/")
  dir.create(path, showWarnings = FALSE)
  for (i in c(1:length(extract))){
    if (is.na(extract[i]$hgncSymbol@values)){
      png(file=paste0(path, sample, "-QQ-", i, "-", extract[i]$type, ".png"))
    }
    else{
      png(file=paste0(path, sample, "-QQ-", i, "-", extract[i]$hgncSymbol, "-", extract[i]$type, ".png"))
    }
    print(plotQQ(fds, result=extract[i]))
    dev.off()
  }
}

for (sample in samples){
  extract <- res[res$sampleID == sample]
  write.csv(extract, paste0("./results/", sample, "/", sample, "-result.csv"))
  path = paste0("./results/", sample, "/Expression/")
  dir.create(path, showWarnings = FALSE)
  for (i in c(1:length(extract))){
      if (is.na(extract[i]$hgncSymbol@values)){
        png(file=paste0(path, "/", sample, "-EX-", i, "-", extract[i]$type, ".png"))
      }
      else{
        png(file=paste0(path, "/", sample, "-EX-", i, "-", extract[i]$hgncSymbol, "-", extract[i]$type, ".png"))
      }
      print(plotExpression(fds, result=extract[i]))
      dev.off()
  }  
}

for (sample in samples){
  extract <- res[res$sampleID == sample]
  path = paste0("./results/", sample, "/Pred_vs_Obs")
  dir.create(path, showWarnings = FALSE)
  for (i in c(1:length(extract))) {
    if (is.na(extract[i]$hgncSymbol@values)){
      png(file=paste0(path, "/", sample, "-PvO-", i, "-", extract[i]$type, ".png"))
    }
    else{
      png(file=paste0(path, "/", sample, "-PvO-", i, "-", extract[i]$hgncSymbol, "-", extract[i]$type, ".png"))
    }
    print(plotExpectedVsObservedPsi(fds, result=extract[i]))
    dev.off()
  }      
}
