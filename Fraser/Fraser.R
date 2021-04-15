# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input directory of bam files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character"),
  make_option(c("-w", "--work"), type="character", default=NULL, 
              help="working directory for intermediate results", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", default=10, 
              help="number of cores to use (default is 10 or maximum cores available, whichever is smaller)", metavar="integer"),
  make_option(c("-f", "--FDR"), type="double", default=0.05, 
              help="FDRcutoffs (default is 0.05)", metavar="double"),
  make_option(c("-g", "--graphs"), type="logical", default=TRUE, 
              help="TRUE to output graphs or FALSE to turn off (default is TRUE)", metavar="logical")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$work)){
  print_help(opt_parser)
  stop("Input, workdir & output must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$input)){
  print_help(opt_parser)
  stop("Input must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$work)){
  dir.create(opt$work)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$input <- ifelse(substr(opt$input, nchar(opt$input), nchar(opt$input))=='/', opt$input, paste0(opt$input,'/'))
opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))
opt$work <- ifelse(substr(opt$work, nchar(opt$work), nchar(opt$work))=='/', opt$work, paste0(opt$work,'/'))

if (opt$graphs){
  dir.create(paste0(opt$output, "graphs"))
  VCpath = paste0(opt$output, "graphs/Volcano/")
  QQpath = paste0(opt$output, "graphs/QQ/")
  EXpath = paste0(opt$output, "graphs/Expression/")
  dir.create(VCpath)
  dir.create(QQpath)
  dir.create(EXpath)
}

# main --------------------------------------------------------------------
library(FRASER)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

register(MulticoreParam(workers=min(opt$cores, multicoreWorkers())))

# extract bam files
bam.all <- list.files(path = opt$input, pattern = "sorted.dedup.bam$", all.files = TRUE, recursive = TRUE)

# create file table
dataInfo <- data.table(sampleID = str_match(basename(bam.all), "(.*?).sorted.dedup.bam")[,2], bamFile = paste0(opt$input, bam.all), pairedEnd=TRUE)

settings <- FraserDataSet(colData=dataInfo, workingDir=opt$work)
fds <- countRNAData(settings, NcpuPerSample = min(opt$cores, multicoreWorkers()))

# compute stats
fds <- calculatePSIValues(fds)

# filtering junction with low expression
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample=20, minDeltaPsi=0.0, filter=TRUE)

q <- c(bestQ(fds, type="psi5"), bestQ(fds, type="psi3"), bestQ(fds, type="theta"))
names(q) <- c("psi5", "psi3", "theta")
fds <- FRASER(fds, q=q, implementation="PCA", BPPARAM = bpparam())

# alternatively, we also provide a way to use biomart for the annotation:
fds <- annotateRanges(fds)
saveFraserDataSet(fds, dir=opt$work)

# get results: we recommend to use an FDR cutoff 0.05
res <- results(fds, zScoreCutoff=NA, padjCutoff=opt$FDRcutoff, deltaPsiCutoff=0.3)

# output results & graphs by sample
samples <- unique(res$sampleID)
samples <- samples[!grepl("Control", samples)]
types = c("theta", "psi5", "psi3")

cl = makeCluster(opt$core)
registerDoParallel(cl)

system.time( foreach (s = 1:length(samples) ) %dopar% {
  library(FRASER)
  library(dplyr)
  
  res_s <- res[res$sampleID == samples[s]] %>% replace(is.na(.), ".")
  write.csv(res_s, paste0(opt$output, samples[s], '.FRASER.result.csv'), row.names=FALSE)
  
  if (opt$graphs) {
    
    for (type in types){
      png(file=paste0(VCpath, samples[s], ".Volcano.", type, ".png"))
      print(plotVolcano(fds, sampleID=samples[s], type=type))
      dev.off()
    }
    
    extract <- res[res$sampleID == samples[s]]
    
    for (i in c(1:length(extract))){
      if (is.na(extract[i]$hgncSymbol@values)){
        png(file=paste0(QQpath, samples[s], ".QQ-", i, ".", extract[i]$type, ".png"))
        print(plotQQ(fds, result=extract[i]))
        dev.off()
        
        png(file=paste0(EXpath, samples[s], ".EX-", i, ".", extract[i]$type, ".png"))
        print(plotExpression(fds, result=extract[i]))
        dev.off()
      }
      else{
        png(file=paste0(QQpath, samples[s], ".QQ-", i, ".", extract[i]$type, ".", extract[i]$hgncSymbol, ".png"))
        print(plotQQ(fds, result=extract[i]))
        dev.off()
        
        png(file=paste0(EXpath, samples[s], ".EX-", i, ".", extract[i]$type, ".", extract[i]$hgncSymbol, ".png"))
        print(plotQQ(fds, result=extract[i]))
        dev.off()     
      }
    }  
  }
})

stopCluster(cl)
