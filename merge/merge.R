# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-m", "--mae"), type="character", default=NULL, 
              help="input MAE directory", metavar="character"),
  make_option(c("-f", "--fraser"), type="character", default=NULL, 
              help="input FRASER directory", metavar="character"),
  make_option(c("-u", "--outrider"), type="character", default=NULL, 
              help="input OUTRIDER directory", metavar="character"),                
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$mae) | is.null(opt$fraser) | is.null(opt$outrider) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Input & output must be supplied.", call.=FALSE)
}

if (!file_test("-d", opt$mae)){
  print_help(opt_parser)
  stop("Input MAE must be a directory.", call.=FALSE)
}
if (!file_test("-d", opt$fraser)){
  print_help(opt_parser)
  stop("Input FRASER must be a directory.", call.=FALSE)
}
if (!file_test("-d", opt$outrider)){
  print_help(opt_parser)
  stop("Input OUTRIDER must be a directory.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$mae <- ifelse(substr(opt$mae, nchar(opt$mae), nchar(opt$mae))=='/', opt$mae, paste0(opt$mae,'/'))
opt$fraser <- ifelse(substr(opt$fraser, nchar(opt$fraser), nchar(opt$fraser))=='/', opt$fraser, paste0(opt$fraser,'/'))
opt$outrider <- ifelse(substr(opt$outrider, nchar(opt$outrider), nchar(opt$outrider))=='/', opt$outrider, paste0(opt$outrider,'/'))
opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))

# functions --------------------------------------------------------
readInput <- function(maeDir, fraserDir, outriderDir, sample){
  # ==== Input FRASER ====
  fraser <- read.csv(paste0(fraserDir, sample, '.FRASER.result.csv')) %>% #read data table
    subset(select = -c(X, sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("FRASER_", colname)) %>% #add prefix
    na_if(".")
  # ==== Input MAE ====
  mae <- read.csv(paste0(maeDir, sample, '.MAE.result.csv')) %>% #read data table
    subset(select = -c(X, sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("MAE_", colname)) %>% #add prefix
    na_if(".")
  # ==== Input OUTRIDER ====
  outrider <- read.csv(paste0(outriderDir, sample, '.OUTRIDER.result.csv')) %>% #read data table
    subset(select = -c(X, sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("OUTRIDER_", colname)) #add prefix
  # subset(OUTRIDER_aberrant == TRUE)
  
  # query range with gene name
  queryResult <- getBM(attributes=c('start_position', 'end_position', 'external_gene_name'),
                       filters = c('external_gene_name'), 
                       values = outrider$OUTRIDER_geneID, mart=ensembl) %>% #query
    rename(OUTRIDER_geneID = external_gene_name) %>% #rename column
    unique()
  outrider <- merge(outrider, queryResult, by ='OUTRIDER_geneID', all=T) %>%
    rename(OUTRIDER_start = start_position) %>% #rename column
    rename(OUTRIDER_end = end_position) %>% #rename column
    relocate(OUTRIDER_start, OUTRIDER_end, .after=OUTRIDER_geneID) %>% #reorder column
    na_if(".")
  
  return(list("mae"=mae, "fraser"=fraser, "outrider"=outrider))
}
mergeAll <- function(mae, fraser, outrider){
  # create chromosome/contig list
  chromosomes <- sort(unique( append( unique(fraser$FRASER_seqnames), unique(mae$MAE_contig))))
  
  # create empty df to store merged data
  merged <- data.frame( matrix( ncol = ncol(mae)+ncol(fraser), nrow = 0 ) )
  colnames(merged) <- c( colnames(mae), colnames(fraser) )
  
  for ( c in c(1:length(chromosomes))){
    fraser_c<- filter(fraser, FRASER_seqnames==chromosomes[c])
    mae_c <- filter(mae, MAE_contig==chromosomes[c])
    merged_tmp <- fuzzy_full_join(mae_c, fraser_c, 
                                  by = c("MAE_position" = "FRASER_start", "MAE_position" = "FRASER_end"), 
                                  match_fun = list(`>=`, `<=`))
   
    merged_mae <- merged_tmp[!is.na(merged_tmp$MAE_position), ]
    merged_fraser <- merged_tmp[!is.na(merged_tmp$FRASER_start), ]
    rm(merged_tmp)
    
    merged_mae <- fuzzy_full_join(merged_mae, outrider, 
                                  by = c("MAE_position" = "OUTRIDER_start", "MAE_position" = "OUTRIDER_end"),
                                  match_fun = list(`>=`, `<=`))
    merged_fraser <- fuzzy_full_join(merged_fraser, outrider, 
                                  by = c("FRASER_start" = "OUTRIDER_start", "FRASER_end" = "OUTRIDER_end"),
                                  match_fun = list(`>=`, `<=`))
    
    merged <- rbind(merged, merged_mae, merged_fraser)
  }
  
  return(merged)
}
# main --------------------------------------------------------
library(plyr)
library(tidyverse)
library(stringr)
library(biomaRt)
library(gdata) #for cbindX()
library(fuzzyjoin)

# extract sample names
mae.all <- list.files(path = opt$mae, pattern = "*.MAE.result.csv$", full.name = FALSE, recursive = FALSE)
samples <- str_match(mae.all, "(.*?).MAE.result.csv")[, 2]

# define dataset for biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 
gene.list <- getBM(attributes=c('start_position', 'end_position', 'external_gene_name'), mart=ensembl) %>%
  arrange(start_position, external_gene_name)

# loop all samples
for (s in c(1:length(samples))){
  print(paste0("merging ", samples[s]))
  
  data <- readInput(opt$mae, opt$fraser, opt$outrider, samples[s])
  
  # merge
  merged_result <- mergeAll(data$mae, data$fraser, data$outrider)
  
  # Post processing of merged_result
  merged_result <- merged_result %>%
    replace(is.na(.), ".") %>% # NAs to "."
    mutate(sampleID = samples[s]) %>%
    mutate(seqnames = ifelse(MAE_contig == ".", FRASER_seqnames, MAE_contig)) %>%
    subset(select = -c(MAE_contig, FRASER_seqnames)) %>%
    mutate(geneID = ifelse(OUTRIDER_geneID == ".", FRASER_hgncSymbol, OUTRIDER_geneID)) %>%
    relocate(sampleID, geneID, seqnames) %>%
    arrange(seqnames, OUTRIDER_start, FRASER_start, MAE_position)
  
  # write output
  write_tsv(merged_result, paste0(opt$output, samples[s], '.rnaseq.merge.tsv'))
}
