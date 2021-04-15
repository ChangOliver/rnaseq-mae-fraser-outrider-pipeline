# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input directory of htseq-count files", metavar="character"),
  make_option(c("-c", "--case"), type="character", default=NULL, 
              help="case file to be processed", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character"),
  make_option(c("-f", "--FC"), type="double", default=2.5, 
              help="FCCutoffs (default is 2.5)", metavar="double"),
  make_option(c("-p", "--padj"), type="double", default=2, 
              help="padjCutoff (default is 2)", metavar="double")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output) | is.null(opt$case) ){
  print_help(opt_parser)
  stop("Input, case & output must be supplied.", call.=FALSE)
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

# functions --------------------------------------------------------
get.translation <- function(res, col){
  idSymbolTable <- read.csv(file = "./translation_table.csv")
  translationDict <- as.vector(idSymbolTable$Gene.name)
  names(translationDict) <- idSymbolTable$Ensembl.ID
  
  for (i in c(1:nrow(res))){
    res[i, col] <- translationDict[[res[[i, 1]]]]
  }
  
  return(res)
}
get.count <- function(ctrlDir, case){
  
  # Get htseq-count files
  ctrl <- grep("RNA-Control-*",list.files(ctrlDir),value=TRUE)
  case <- basename(case)
  files <- append(ctrl, case)
  
  # Extract id
  names <- append(sub("(.*).htseq-count.txt", "\\1", ctrl), sub("(.*).htseq-count.txt", "\\1", case))
  
  # Mark case or control
  conditions <- rep("Control", times=length(files))
  conditions[length(files)] <- "Case"
  
  # create file info table, set conditions as factor
  table <- data.frame(sampleName = names, fileName = files, condition = conditions)
  table$condition <- relevel(factor(table$condition), ref = "Control")
  
  return(list("table" = table, "condition" = conditions, "names" = names))
}
run.DESeq2 <- function(dir, table, condition){
  
  # Create DESeq object
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table, directory = dir, design= ~ condition)
  
  # Run differential expression and save dds object
  ddsHTSeq <- DESeq(ddsHTSeq, minReplicatesForReplace=Inf)
  # saveRDS(ddsHTSeq, file=paste0(getwd(), "/DESeq2_", case, "_dataset.rds")) 
  return(ddsHTSeq)
}

extract.result <- function(dds){
  
  # Extract result
  res <- DESeq2::results(dds, contrast=c("condition", "Case", "Control"), cooksCutoff=FALSE, independentFiltering=FALSE)
  res <- data.frame(res, row.names = res@rownames)
  res <- subset(res, !is.na(res$padj))
  res[ , "name"] <- rownames(res)
  res <- rownames_to_column(res[c(7, 1:6)]) %>% get.translation(2)

  return(res)
}

plot.volcano <- function(resList, resDir, case, fcCutoff, padjCutoff){
  
  case <- sub("(.*).htseq-count.txt", "\\1", basename(case))
  
  imgPath <- paste0(resDir, "/graphs/")
  dir.create(imgPath, showWarnings = FALSE)
  
  resList <- mutate(resList, 
                    group = ifelse( (log2FoldChange >= fcCutoff & -log10(padj) > padjCutoff), "R", 
                                    ifelse( (log2FoldChange <= (-fcCutoff) & -log10(padj) > padjCutoff), "G", "B")))
  
  resList <- resList[order(resList$group, decreasing = TRUE), ]
  volcano <- ggplot(resList, mapping = aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(aes(color = group)) + 
    scale_color_manual(values = c("R" = "#EC7063", "G" = "#48C9B0", "B" = "#85929E")) + #1:red 2:green 3:gray
    ggtitle(case) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") + 
    theme(legend.position = "none", 
          plot.title = element_text(size = rel(1.5), hjust = 0.5), 
          axis.title = element_text(size = rel(1.25)))
  ggsave(paste0(imgPath, case, ".volcano.png"), width = 9.93, height = 6.86, dpi = 300)
  
  label_volcano <- volcano + geom_label_repel(data = filter(resList, group != "B"), 
                                              aes(label = as.character(name)), 
                                              size = 3) 
  ggsave(paste0(imgPath, case, ".label.volcano.png"), width = 9.93, height = 6.86, dpi =300)
  
  return(resList)
}

# main --------------------------------------------------------
library(dplyr)
library(DESeq2)
library(tibble) #rownames_to_column
library(ggplot2)
library(ggrepel)

case <- sub("(.*).htseq-count.txt", "\\1", basename(opt$case))

# Fetch htseq data
print(paste0("Running DESeq2 analysis on ", case))
returnVal <- get.count(opt$input, opt$case)
dds <- run.DESeq2(opt$input, returnVal$table, returnVal$condition)

message("extracting result")
desResult <- extract.result(dds)

message("generating volcano plot & csv")
desResult <- plot.volcano(desResult, opt$output, opt$case, opt$FC, opt$padj)
write.csv(desResult, paste0(opt$output, case, ".DESeq2.result.csv"), row.names = FALSE)