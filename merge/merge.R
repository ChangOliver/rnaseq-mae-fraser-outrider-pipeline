# argument parsing --------------------------------------------------------
library(optparse)

option_list = list(
  make_option(c("-m", "--mae"), type="character", default=NULL, 
              help="input MAE.result.csv", metavar="character"),
  make_option(c("-f", "--fraser"), type="character", default=NULL, 
              help="input FRASER.result.csv", metavar="character"),
  make_option(c("-u", "--outrider"), type="character", default=NULL, 
              help="input OUTRIDER.result.csv", metavar="character"),
  make_option(c("-d", "--dataset"), type="numeric", default=NULL, 
              help="37 for GRCh37 and 38 for GRCh38", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory to store results", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$mae) | is.null(opt$fraser) | is.null(opt$outrider) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Input & output must be provided.", call.=FALSE)
}

if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("Choose a dataset to use.", call.=FALSE)
}

if (!file_test("-d", opt$output)){
  dir.create(opt$output)
}

opt$output <- ifelse(substr(opt$output, nchar(opt$output), nchar(opt$output))=='/', opt$output, paste0(opt$output,'/'))

# functions --------------------------------------------------------
readInput <- function(maeFile, fraserFile, outriderFile){
  
  # ==== Input mart ====
  mart <- read.csv(ifelse(opt$dataset==37, "mart37_export.csv", "mart38_export.csv")) %>%
          rename(OUTRIDER_geneID = GeneName, OUTRIDER_chr=Chromosome, OUTRIDER_start=GeneStart, OUTRIDER_end=GeneEnd)
  
  # ==== Input FRASER ====
  fraser <- read.csv(fraserFile) %>% #read data table
    subset(select = -c(sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("FRASER_", colname)) %>% #add prefix
    na_if(".")
  # ==== Input MAE ====
  mae <- read.csv(maeFile) %>% #read data table
    subset(select = -c(sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("MAE_", colname)) %>% #add prefix
    na_if(".")
  # ==== Input OUTRIDER ====
  outrider <- read.csv(outriderFile) %>% #read data table
    subset(select = -c(sampleID)) %>% #remove redundant columns
    rename_all( function(colname) paste0("OUTRIDER_", colname)) #add prefix

  outrider <- merge(outrider, mart, by ='OUTRIDER_geneID', all.x=T) %>%
    relocate(OUTRIDER_chr, OUTRIDER_start, OUTRIDER_end, .after=OUTRIDER_geneID) %>% #reorder column
    na_if(".")
  
  return(list("mae"=mae, "fraser"=fraser, "outrider"=outrider))
}
mergeAll <- function(mae, fraser, outrider){
  # create chromosome/contig list
  chromosomes <- sort(unique(append( append(fraser$FRASER_seqnames, mae$MAE_contig), outrider$OUTRIDER_chr )))
  
  # create empty df to store merged data
  merged <- data.frame( matrix( ncol = ncol(mae)+ncol(fraser)+ncol(outrider), nrow = 0 ) )
  colnames(merged) <- c( colnames(mae), colnames(fraser), colnames(outrider) )
  
  # sort
  mae <- arrange(mae, mae$MAE_contig, mae$MAE_position)
  fraser <- arrange(fraser, fraser$FRASER_seqnames, fraser$FRASER_start, fraser$FRASER_end)
  outrider <- arrange(outrider, outrider$OUTRIDER_chr, outrider$OUTRIDER_start, outrider$OUTRIDER_end)
  
  for ( c in c(1:length(chromosomes))){
    ptr <- vector(mode="list", length=5)
    names(ptr) <- c("m", "f", "o", "mrollback", "frollback")
    ptr[[1]] = ptr[[2]] = ptr[[3]] = 1
    ptr[[4]] = ptr[[5]] = 0
    
    fraser_c<- filter(fraser, FRASER_seqnames==chromosomes[c])
    mae_c <- filter(mae, MAE_contig==chromosomes[c])
    outrider_c <- filter(outrider, OUTRIDER_chr==chromosomes[c])
    
    mae_c.nrow <- nrow(mae_c)
    fraser_c.nrow <- nrow(fraser_c)
    outrider_c.nrow <- nrow(outrider_c)
    
    pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
    fstart <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_start[ptr$f] else Inf
    fend <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_end[ptr$f] else Inf
    ostart <- if (ptr$o <= outrider_c.nrow) outrider_c$OUTRIDER_start[ptr$o] else Inf
    oend <- if (ptr$o <= outrider_c.nrow) outrider_c$OUTRIDER_end[ptr$o] else Inf
    
    while (ptr$m <= mae_c.nrow | ptr$f <= fraser_c.nrow | ptr$o <= outrider_c.nrow){
      
      M.in.F <- if (ptr$m <= mae_c.nrow & ptr$f <= fraser_c.nrow) (pos > fstart & pos < fend) else FALSE
      M.in.O <- if (ptr$m <= mae_c.nrow & ptr$o <= outrider_c.nrow) (pos > ostart & pos < oend) else FALSE
      F.in.O <- if (ptr$f <= fraser_c.nrow & ptr$o <= outrider_c.nrow) (fstart > ostart & fend < oend) else FALSE
      
      if (M.in.F & M.in.O & F.in.O) {
        merged <- rbindlist( list(merged, cbind( mae_c[ ptr$m, ], fraser_c[ ptr$f, ], outrider_c[ ptr$o, ])), fill=TRUE)
      } 
      else {
        if (M.in.F) { merged <- rbindlist( list(merged, cbind( mae_c[ ptr$m, ], fraser_c[ ptr$f, ])), fill=TRUE) }
        if (M.in.O) { merged <- rbindlist( list(merged, cbind( mae_c[ ptr$m, ], outrider_c[ ptr$o, ])), fill=TRUE) }
        if (F.in.O) { merged <- rbindlist( list(merged, cbind( fraser_c[ ptr$f, ], outrider_c[ ptr$o, ])), fill=TRUE) }
      }
      
      min <- which.min( c(pos, fstart, ostart))
      
      if (min == 1) {
        merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, mae_c[ptr$m, ]), fill=TRUE) else merged
        ptr$m = ptr$m + 1
        pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
      } 
      else if (min == 2) {
        if (pos < fend) {
          merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, mae_c[ptr$m, ]), fill=TRUE) else merged
          ptr$mrollback = ptr$m
          ptr$m = ptr$m + 1 
          pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
        } 
        else { 
            merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, fraser_c[ptr$f, ]), fill=TRUE) else merged
            ptr$f = ptr$f + 1 
            ptr$m = if (ptr$mrollback != 0) ptr$mrollback else ptr$m
            ptr$mrollback = 0
            pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
            fstart <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_start[ptr$f] else Inf
            fend <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_end[ptr$f] else Inf
          }
      } 
      else if (min == 3) {
        ptr$mrollback = if (ptr$mrollback == 0) ptr$m else ptr$mrollback
        ptr$frollback = if (ptr$frollback == 0) ptr$f else ptr$frollback
        if (pos < oend) {
          merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, mae_c[ptr$m, ]), fill=TRUE) else merged
          ptr$m = ptr$m + 1 
          pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
        } 
        else if (fend < oend) {
          merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, fraser_c[ptr$f, ]), fill=TRUE) else merged
          ptr$f = ptr$f + 1 
          fstart <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_start[ptr$f] else Inf
          fend <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_end[ptr$f] else Inf
        } 
        else {
          merged <- if (!M.in.F & !M.in.O & !F.in.O) rbindlist(list(merged, outrider_c[ptr$o, ]), fill=TRUE) else merged
          ptr$o = ptr$o + 1
          ptr$m = if (ptr$mrollback != 0) ptr$mrollback else ptr$m
          ptr$f = if (ptr$frollback != 0) ptr$frollback else ptr$f
          ptr$mrollback = ptr$frollback = 0
          pos <- if (ptr$m <= mae_c.nrow) mae_c$MAE_position[ptr$m] else Inf
          fstart <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_start[ptr$f] else Inf
          fend <- if (ptr$f <= fraser_c.nrow) fraser_c$FRASER_end[ptr$f] else Inf
          ostart <- if (ptr$o <= outrider_c.nrow) outrider_c$OUTRIDER_start[ptr$o] else Inf
          oend <- if (ptr$o <= outrider_c.nrow) outrider_c$OUTRIDER_end[ptr$o] else Inf
        }
      }
    }
  }
  return(merged)
}

# main --------------------------------------------------------
library(plyr)
library(tidyverse)
library(stringr)
library(data.table) #for rbindlist()

# extract sample names
sample <- str_match(basename(opt$mae), "(.*?).MAE.result.csv")[, 2]

print(paste0("merging ", sample))
data <- readInput(opt$mae, opt$fraser, opt$outrider)

# merge
merged_result <- mergeAll(data$mae, data$fraser, data$outrider)

# Post processing of merged_result
merged_result <- merged_result %>%
  replace(is.na(.), ".") %>% # NAs to "."
  mutate(sampleID = sample) %>%
  mutate(seqnames = ifelse(MAE_contig == ".",
                           ifelse(FRASER_seqnames == ".", OUTRIDER_chr, FRASER_seqnames), MAE_contig)) %>%
  mutate(geneID = ifelse(OUTRIDER_geneID == ".", FRASER_hgncSymbol, OUTRIDER_geneID)) %>%
  subset(select = -c(MAE_contig, FRASER_seqnames, OUTRIDER_chr, OUTRIDER_geneID)) %>%
  relocate(sampleID, geneID, seqnames) %>%
  arrange(seqnames, MAE_position, FRASER_start, OUTRIDER_start)

# write output
write_tsv(merged_result, paste0(opt$output, sample, '.rnaseq.merge.tsv'))
