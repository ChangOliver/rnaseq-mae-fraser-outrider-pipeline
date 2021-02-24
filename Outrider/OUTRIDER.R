#############################
# Load libraries
#############################
print("Loading libraries...")
suppressMessages(library(dplyr))
suppressMessages(library(OUTRIDER))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

#############################
# Function declaration
#############################
combine.htseq <- function(cntDir, pat){
  # adapted from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4#Combine_individual_HTSeq_files_from_the_.27all.27_mapping_series
  
  tophat.all <- list.files(path = cntDir, pattern = pat, all.files = TRUE, recursive = TRUE)
  
  # we choose the 'all' series
  myfiles <- tophat.all
  DT <- list()
  
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  for (i in 1:length(myfiles) ) {
	infile = paste(cntDir, myfiles[i], sep = "/")
	DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
	cnts <- gsub("(.*)/6.STATS/(.*).htseq-count.txt", "\\1", myfiles[i])
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
  
  # # write data to file
  # write.csv(data.all, file = paste0(cntDir, outFile))
  
  # cleanup intermediate objects
  rm(y, z, i, DT)

  return(data.all)
}

run.OUTRIDER <- function(countDir, ctsTable){
  #### Load data
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

  #### Create OUTRIDER object
  ods <- OutriderDataSet(countData=ctsMatrix)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- OUTRIDER(ods)
  
  return(ods)
}

plot_QQ_ExRank <- function(res, resDir){
  
  res <- subset(res, res$aberrant)
  if (nrow(res) == 0){
	# print("No aberrant genes found")
	return()
  }
  
  imgPath <- paste0(resDir, "/graphs/")
  dir.create(paste0(imgPath, "/QQ"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(imgPath, "/Expression"), recursive = TRUE, showWarnings = FALSE)
	
  for (i in c(1:nrow(res))){
	
	  png(paste0(imgPath, "/QQ/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 715, height = 494)
	  QQ <- plotQQ(ods, res[i, geneID])
	  dev.off()
	
	  # expression rank of a gene with outlier events
	  EX <- plotExpressionRank(ods, res[i, geneID], basePlot=TRUE)
	  ggsave(paste0(imgPath, "/Expression/", res[i]$sampleID, ".", res[i]$geneID, ".png"), width = 9.93, height = 6.86, dpi = 300)
  }
}

#############################
# Main
#############################
print("Starting analysis")
args <- commandArgs(trailingOnly = TRUE)

dataPath <- args[[1]]

data.all <- combine.htseq(dataPath, "htseq-count.txt$")
  
#### Run OUTRIDER
print("Running OUTRIDER analysis...")
if (file_test("-d", dataPath)){
	ods <- run.OUTRIDER(dataPath, data.all)
} else{
	ods <- readRDS(file=dataPath)
}

print("extracting result")
odsResult <- OUTRIDER::results(ods, all = TRUE)

print("Dumping results & graphs")
samples <- unique(odsResult$sampleID)
for (s in c(1:length(samples))){
  res_s <- odsResult %>% 
		filter(grepl(samples[s], sampleID)) %>%
		replace(is.na(.), ".")
  savePath <-  paste0(dataPath, "/", samples[s], "/22.OUTRIDER/")
  dir.create(savePath, showWarnings = FALSE)
  write.csv(res_s, paste0(savePath, samples[s], '.OUTRIDER.result.csv'), row.names=FALSE)

  plot_QQ_ExRank(res_s, savePath)
}
