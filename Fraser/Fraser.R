# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input bam directory path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory path", metavar="character")
  make_option(c("-t", "--tmp"), type="character", default=NULL, 
              help="tmp directory path", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$input <- ifelse(substr(opt$input, nchar(opt$input), nchar(opt$input))=='/', opt$input, paste0(opt$input,'/'))
opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))
opt$tmp <- ifelse(substr(opt$tmp, nchar(opt$tmp), nchar(opt$tmp))=='/', opt$tmp, paste0(opt$tmp,'/'))

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$tmp)){
  print_help(opt_parser)
  stop("Arguments must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$input)){
  print_help(opt_parser)
  stop("Input must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

if (!file_test("-d", opt$tmp)){
  dir.create(opt$tmp)
}

# main --------------------------------------------------------------------
library(FRASER)
library(dplyr)
library(stringr)

# extract bam files
bam.all <- list.files(path = opt$input, pattern = "dedup.bam$", all.files = TRUE, recursive = TRUE)

# create file table
dataInfo <- data.table(sampleID = str_match(bam.all, "(.*?).sorted.dedup.bam")[,2], bamFile = paste0(opt$input, bam.all), pairedEnd=TRUE)

settings <- FraserDataSet(colData=dataInfo, workingDir=opt$tmp)
register(MulticoreParam(workers=min(10, multicoreWorkers())))
fds <- countRNAData(settings)

# compute stats
fds <- calculatePSIValues(fds)

# filtering junction with low expression
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20, minDeltaPsi=0.0, filter=TRUE)

# use PCA to speed up the tutorial
fds <- FRASER(fds, q=bestQ(fds, type="psi5"), implementation="PCA")

# alternatively, we also provide a way to use biomart for the annotation:
fds <- annotateRanges(fds)
# get results: we recommend to use an FDR cutoff 0.05
res <- results(fds, zScoreCutoff=NA, padjCutoff=0.05, deltaPsiCutoff=0.3)

# output results by sample
samples <- unique(res$sampleID)
samples <- samples[!grepl("Control", samples)]
for (s in c(1:length(samples))){
  res_s <- res[res$sampleID == samples[s]] %>% replace(is.na(.), ".")
  write.csv(res_s, paste0(opt$output, samples[s], '.FRASER.result.csv'), row.names=FALSE)
}

# output graphs --------------------------------------------------------------------
types = c("theta", "psi5", "psi3")
dir.create(paste0(opt$output, "graphs"))

for (s in c(1:length(samples))){
  path = paste0(opt$output, "graphs/Volcano/")
  dir.create(path)
  for (type in types){
    png(file=paste0(path, samples[s], ".Volcano.", type, ".png"))
    print(plotVolcano(fds, sampleID=samples[s], type=type))
    dev.off()
  }
}

for (s in c(1:length(samples))){
  extract <- res[res$sampleID == samples[s]]
  QQpath = paste0(opt$output, "graphs/QQ/")
  EXpath = paste0(opt$output, "graphs/Expression/")
  POpath = paste0(opt$output, "graphs/Prediction_vs_Observation/")
  dir.create(QQpath)
  dir.create(EXpath)
  dir.create(POpath)

  for (i in c(1:length(extract))){
    if (is.na(extract[i]$hgncSymbol@values)){
      png(file=paste0(QQpath, samples[s], ".QQ-", i, ".", extract[i]$type, ".png"))
      print(plotQQ(fds, result=extract[i]))
      dev.off()

      png(file=paste0(EXpath, samples[s], ".EX-", i, ".", extract[i]$type, ".png"))
      print(plotExpression(fds, result=extract[i]))
      dev.off()

      png(file=paste0(POpath, samples[s], ".PvO-", i, ".", extract[i]$type, ".png"))
      print(plotExpectedVsObservedPsi(fds, result=extract[i]))
      dev.off()
    }
    else{
      png(file=paste0(QQpath, samples[s], ".QQ-", i, ".", extract[i]$type, ".", extract[i]$hgncSymbol, ".png"))
      print(plotQQ(fds, result=extract[i]))
      dev.off()

      png(file=paste0(EXpath, samples[s], ".EX-", i, ".", extract[i]$type, ".", extract[i]$hgncSymbol, ".png"))
      print(plotQQ(fds, result=extract[i]))
      dev.off()    

      png(file=paste0(POpath, samples[s], ".PvO-", i, ".", extract[i]$type, ".", extract[i]$hgncSymbol, ".png"))
      print(plotExpectedVsObservedPsi(fds, result=extract[i]))
      dev.off()      
    }
  }
}
