# argument parsing --------------------------------------------------------
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--case"), type="character", default=NULL, 
              help="case name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character"),
  make_option(c("-w", "--work"), type="character", default=NULL, 
              help="working directory for intermediate results", metavar="character"),
  make_option(c("-f", "--FDR"), type="double", default=0.05, 
              help="FDRcutoffs (default is 0.05)", metavar="double"),
  make_option(c("-p", "--dPSI"), type="double", default=0.3, 
              help="deltaPsiCutoffs (default is 0.3)", metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$case) | is.null(opt$output) | is.null(opt$work)){
  print_help(opt_parser)
  stop("Case name, workdir & output must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$work)){
  stop("Workdir does not exist.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))
opt$work <- ifelse(substr(opt$work, nchar(opt$work), nchar(opt$work))=='/', opt$work, paste0(opt$work,'/'))

VCpath = paste0(opt$output, "Volcano/")
dir.create(VCpath)

# main --------------------------------------------------------------------
suppressMessages(library(FRASER))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

fds <- loadFraserDataSet(opt$work)

# get results: we recommend to use an FDR cutoff 0.05
res <- results(fds, zScoreCutoff=NA, padjCutoff=opt$FDR, deltaPsiCutoff=opt$dPSI)

# output results & graphs by sample
sample <- opt$case
types = c("theta", "psi5", "psi3")

if (sample %in% res$sampleID ){

  res_s <- res[res$sampleID == sample] %>% replace(is.na(.), ".")
  write.csv(res_s, paste0(opt$output, sample, '.FRASER.adjusted.result.csv'), row.names=FALSE)

  for (type in types){
    png(file=paste0(VCpath, sample, ".Volcano.", type, ".png"))
    print(plotVolcano(fds, sampleID=sample, type=type, padjCutoff=opt$FDR, deltaPsiCutoff=opt$dPSI))
    dev.off()

    png(file=paste0(VCpath, sample, ".Volcano.", type, ".label.png"))
    print(plotVolcano(fds, sampleID=sample, type=type, label="aberrant", padjCutoff=opt$FDR, deltaPsiCutoff=opt$dPSI))
    dev.off()
  }  
} else {
  message("Case name not found.")
}
