#############################
# Load libraries
#############################
print("Loading libraries...")
suppressMessages(library(dplyr))
suppressMessages(library(OUTRIDER))
suppressMessages(library(biomaRt))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel)

#############################
# Function declaration
#############################
combine.htseq <- function(cntDir, pat, outFile){
  # adapted from https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4#Combine_individual_HTSeq_files_from_the_.27all.27_mapping_series
  
  tophat.all <- list.files(path = cntDir,
                           pattern = pat,
                           all.files = TRUE,
                           recursive = FALSE,
                           ignore.case = FALSE,
                           include.dirs = FALSE)
  
  # we choose the 'all' series
  myfiles <- tophat.all
  DT <- list()
  
  # read each file as array element of DT and rename the last 2 cols
  # we created a list of single sample tables
  for (i in 1:length(myfiles) ) {
    infile = paste(cntDir, myfiles[i], sep = "/")
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
  
  # write data to file
  write.csv(data.all, file = paste0(cntDir, outFile))
  
  # cleanup intermediate objects
  rm(y, z, i, DT)
}

run.OUTRIDER <- function(countDir, resDir){
  #### Load data
  cts <- paste0(countDir, "/htseq-counts-all.csv")
  ctsTable <- read.csv(cts, check.names=FALSE, header=TRUE, stringsAsFactors = FALSE)
  colnames(ctsTable)[1] <- "geneID"
  ctsTable <- mutate(ctsTable, geneID = sub("\\..*$", "", geneID))

  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  queryResult <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                       filters = c('ensembl_gene_id'),
                       values = ctsTable$geneID, mart=ensembl) %>% #query
                  mutate(geneID = ensembl_gene_id) %>%
                  subset(select=-c(ensembl_gene_id)) #rename column

  ctsTable <- merge(ctsTable, queryResult, by ='geneID', all=T) %>%
    mutate(geneID = ifelse(is.na(external_gene_name), geneID, external_gene_name)) %>%
    subset(select=-c(external_gene_name))

  ctsMatrix <- data.matrix(ctsTable)[, -1]
  rownames(ctsMatrix) <- ctsTable[,1]

  #### Create OUTRIDER object
  ods <- OutriderDataSet(countData=ctsMatrix)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- OUTRIDER(ods)
  
  # saveRDS(ods, file = paste0(resDir, "/OUTRIDER_dataset.rds"))
  return(ods)
}

plot_QQ_ExRank <- function(res, resDir){
  
  res <- subset(res, res$aberrant)
  if (nrow(res) == 0){
    print("No aberrant genes found")
    return()
  }
  
  imgPath <- paste0(resDir, "/imgs/")
  dir.create(imgPath, showWarnings = FALSE)
    
  for (i in c(1:nrow(res))){
    
      png(paste(paste0(imgPath, "QQ"), res[i]$sampleID, res[i]$geneID, ".png", sep = "_"), width = 715, height = 494)
      QQ <- plotQQ(ods, res[i, geneID])
      dev.off()
    
      # expression rank of a gene with outlier events
      EX <- plotExpressionRank(ods, res[i, geneID], basePlot=TRUE)
      ggsave(paste(paste0(imgPath, "EX"), res[i]$sampleID, res[i]$geneID, ".png", sep = "_"), width = 9.93, height = 6.86, dpi = 300)
  }
}

#############################
# Main
#############################
print("Starting analysis")
args <- commandArgs(trailingOnly = TRUE)

dataPath <- args[[1]]

combine.htseq(dataPath, "*.htseq-count.txt", "/htseq-counts-all.csv")

outriderDir <- "Outrider_results"
dir.create(outriderDir, showWarnings = FALSE, recursive = TRUE)
  
#### Run OUTRIDER
print("Running OUTRIDER analysis...")
if (file_test("-d", dataPath)){
	ods <- run.OUTRIDER(dataPath, outriderDir)
} else{
	ods <- readRDS(file=dataPath)
}

message("extracting result")
odsResult <- OUTRIDER::results(ods, all = TRUE)

message("generating QQ & ExRank plot")
plot_QQ_ExRank(odsResult, outriderDir)

print("Dumping results...")
samples <- unique(odsResult$sampleID)
dir.create("Outrider_results", showWarnings = FALSE)
for (s in c(1:length(samples))){
  res_s <- odsResult %>% 
        filter(grepl(samples[s], sampleID)) %>%
        replace(is.na(.), ".")
  write.csv(res_s, paste0("Outrider_results/", samples[s], '.OUTRIDER.result.csv'), row.names=FALSE)
}
