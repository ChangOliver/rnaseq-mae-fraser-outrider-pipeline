# argument parsing --------------------------------------------------------
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--case"), type="character", default=NULL, 
              help="input directory of case count files", metavar="character"),
  make_option(c("-c", "--control"), type="character", default=NULL, 
              help="input directory of control count files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character"),
  make_option(c("-d", "--dataset"), type="numeric", default=NULL, 
              help="37 for hs37d5 and 38 for GRCh38", metavar="numeric"),
  make_option(c("-t", "--cores"), type="integer", default=10, 
              help="number of cores to use (default is 10 or maximum cores available, whichever is smaller)", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$case) | is.null(opt$control) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Case, control & output must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$case)){
  print_help(opt_parser)
  stop("Case must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$control)){
  print_help(opt_parser)
  stop("Control must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

if (opt$dataset == 37){
  mart <- read.csv("mart37_export.csv")
} else if (opt$dataset == 38){
  mart <- read.csv("mart38_export.csv")
} else if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("Dataset must be provided.", call.=FALSE)
}

opt$case<- ifelse(substr(opt$case, nchar(opt$case), nchar(opt$case))=='/', opt$case, paste0(opt$case,'/'))
opt$control <- ifelse(substr(opt$control, nchar(opt$control), nchar(opt$control))=='/', opt$control, paste0(opt$control,'/'))
opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))

# Function declaration --------------------------------------------------------
combine.htseq <- function(caseDir, ctrlDir){
  # adapted from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4#Combine_individual_HTSeq_files_from_the_.27all.27_mapping_series
  
  cases <- list.files(path = caseDir)
  ctrls <- list.files(path = ctrlDir)
  
  # we choose the 'all' series
  DT <- list()
  
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  for (i in 1:length(cases) ) {
  	infile = paste0(caseDir, cases[i])
  	DT[[cases[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  	cnts <- sub("\\..*", "", basename(cases[i]))
  	colnames(DT[[cases[i]]]) <- c("ID", cnts)
  }
  
  for (i in 1:length(ctrls) ) {
    infile = paste0(ctrlDir, ctrls[i])
    DT[[ctrls[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
    cnts <- sub("\\..*", "", basename(ctrls[i]))
    colnames(DT[[ctrls[i]]]) <- c("ID", cnts)
  }
  
  # merge all elements based on first ID columns
  data <- DT[[cases[1]]]
  
  # we now add each other table with the ID column as key
  for (i in 2:length(cases)) {
	y <- DT[[cases[i]]]
	z <- merge(data, y, by = c("ID"))
	data <- z
  }
  
  for (i in 1:length(ctrls)) {
    y <- DT[[ctrls[i]]]
    z <- merge(data, y, by = c("ID"))
    data <- z
  }
  
  # ID column becomes rownames
  rownames(data) <- data$ID
  data <- data[,-1]
  
  ## add total counts per sample
  data <- rbind(data, tot.counts=colSums(data))
  
  # take all data rows to a new table
  data.all <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]
  
  # cleanup intermediate objects
  rm(y, z, i, DT, data)

  return(data.all)
}

run.OUTRIDER <- function(ctsTable){
  # Load data
  ctsTable <- tibble::rownames_to_column(ctsTable, "geneID")
  ctsTable <- mutate(ctsTable, geneID = sub("\\..*$", "", geneID))

  ctsTable <- merge(ctsTable, mart, by ='geneID', all.x=T) %>%
	mutate(geneID = ifelse(is.na(geneName), geneID, geneName)) %>%
	subset(select=-c(geneName))

  ctsMatrix <- data.matrix(ctsTable)[, -1]
  rownames(ctsMatrix) <- ctsTable[,1]

  # Create OUTRIDER object
  ods <- OutriderDataSet(countData=ctsMatrix)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- OUTRIDER(ods, BPPARAM = bpparam())
  
  return(ods)
}

plot_QQ_ExRank <- function(ods, res, resDir){
  
  # draw aberrant only
  res <- subset(res, res$aberrant)
  if (nrow(res) == 0){
  	return()
  }
  
  imgPath <- paste0(resDir, "graphs/")
  dir.create(paste0(imgPath, "QQ"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(imgPath, "Expression"), showWarnings = FALSE)
	
  for (i in 1:nrow(res)){
	  
    # QQ plot of a gene with outlier events
	  png(paste0(imgPath, "QQ/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 715, height = 494)
	  QQ <- plotQQ(ods, res[i, geneID])
	  dev.off()
	
	  # expression rank of a gene with outlier events
	  EX <- plotExpressionRank(ods, res[i, geneID], basePlot=TRUE)
	  ggsave(paste0(imgPath, "Expression/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 9.93, height = 6.86, dpi = 300)
  }
}

# main --------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(OUTRIDER))
suppressMessages(library(rlist)) #list.append
suppressMessages(library(ggplot2))

register(MulticoreParam(workers=min(opt$cores, multicoreWorkers())))

data.all <- combine.htseq(opt$case, opt$control)
  
# Run OUTRIDER
print("Running OUTRIDER analysis")
ods <- run.OUTRIDER(data.all)
odsResult <- results(ods, all = TRUE)

# output results
print("Dumping results & graphs")
samples <- sub("\\..*", "", basename(list.files(path = opt$case)))

for (s in 1:length(samples)){
  res_s <- odsResult %>% filter(grepl(samples[s], sampleID)) %>% replace(is.na(.), ".")
  write.csv(res_s, paste0(opt$output, samples[s], '.OUTRIDER.result.csv'), row.names=FALSE)
  plot_QQ_ExRank(ods, res_s, opt$output)
}
