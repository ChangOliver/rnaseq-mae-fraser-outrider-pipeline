# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input directory of count files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character")
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

# Function declaration --------------------------------------------------------
combine.htseq <- function(cntDir, pat){
  # adapted from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4#Combine_individual_HTSeq_files_from_the_.27all.27_mapping_series
  
  tophat.all <- list.files(path = cntDir, pattern = pat)
  
  # we choose the 'all' series
  myfiles <- tophat.all
  DT <- list()
  
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  for (i in 1:length(myfiles) ) {
  	infile = paste0(cntDir, myfiles[i])
  	DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  	cnts <- gsub("(.*).htseq-count.txt", "\\1", myfiles[i])
  	colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
  }
  
  # merge all elements based on first ID columns
  data <- DT[[myfiles[1]]]
  
  # we now add each other table with the ID column as key
  for (i in 2:length(myfiles)) {
	y <- DT[[myfiles[i]]]
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
  rm(y, z, i, DT)

  return(data.all)
}

run.OUTRIDER <- function(countDir, ctsTable){
  # Load data
  ctsTable <- tibble::rownames_to_column(ctsTable, "geneID")
  ctsTable <- mutate(ctsTable, geneID = sub("\\..*$", "", geneID))

  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  queryResult <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
					   filters = c('ensembl_gene_id'),
					   values = ctsTable$geneID, mart=ensembl) %>% #query
				  mutate(geneID = ensembl_gene_id) %>%
				  subset(select = -c(ensembl_gene_id)) #rename column

  ctsTable <- merge(ctsTable, queryResult, by ='geneID', all=T) %>%
	mutate(geneID = ifelse(is.na(external_gene_name), geneID, external_gene_name)) %>%
	subset(select=-c(external_gene_name))

  ctsMatrix <- data.matrix(ctsTable)[, -1]
  rownames(ctsMatrix) <- ctsTable[,1]

  # Create OUTRIDER object
  ods <- OutriderDataSet(countData=ctsMatrix)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- OUTRIDER(ods)
  
  return(ods)
}

plot_QQ_ExRank <- function(ods, res, resDir){
  
  # draw aberrant only
  res <- subset(res, res$aberrant)
  if (nrow(res) == 0){
  	return()
  }
  
  imgPath <- paste0(resDir, "graphs/")
  dir.create(paste0(imgPath, "QQ"), recursive = TRUE)
  dir.create(paste0(imgPath, "Expression"))
	
  for (i in c(1:nrow(res))){
	
	  png(paste0(imgPath, "QQ/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 715, height = 494)
	  QQ <- plotQQ(ods, res[i, geneID])
	  dev.off()
	
	  # expression rank of a gene with outlier events
	  EX <- plotExpressionRank(ods, res[i, geneID], basePlot=TRUE)
	  ggsave(paste0(imgPath, "Expression/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 9.93, height = 6.86, dpi = 300)
  }
}

# main --------------------------------------------------------
library(dplyr)
library(OUTRIDER)
library(biomaRt)
library(ggplot2)
library(ggrepel)

data.all <- combine.htseq(opt$input, "htseq-count.txt$")
  
# Run OUTRIDER
print("Running OUTRIDER analysis...")
ods <- run.OUTRIDER(opt$input, data.all)
odsResult <- results(ods, all = TRUE)

# output results
print("Dumping results & graphs")
samples <- unique(odsResult$sampleID)
samples <- samples[!grepl("Control", samples)]

for (s in c(1:length(samples))){
  res_s <- odsResult %>% filter(grepl(samples[s], sampleID)) %>% replace(is.na(.), ".")
  write.csv(res_s, paste0(opt$output, samples[s], '.OUTRIDER.result.csv'), row.names=FALSE)
  plot_QQ_ExRank(ods, res_s, opt$output)
}
