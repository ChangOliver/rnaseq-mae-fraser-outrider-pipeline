library(plyr)
library(tidyverse)
library(stringr)
library(biomaRt)
library(gdata) #for cbindX()

# take arguments
args <- commandArgs(trailingOnly = TRUE)
maeDir <- args[[1]]
fraserDir <- args[[2]]
outriderDir <- args[[3]]

# extract mae iles
mae.all <- list.files(path = maeDir,
                      pattern = "*.MAE.result.csv$",
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)
# extract fraser files
fraser.all <- list.files(path = fraserDir,
                      pattern = "*.FRASER.result.csv$",
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)
# extract outrider files
outrider.all <- list.files(path = outriderDir,
                         pattern = "*.OUTRIDER.result.csv$",
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# ==== Functions ====
mergeAll <- function(fraser_sample, mae_sample, outrider_sample){
  
  # create empty df to store merged data
  merged <- data.frame( matrix( ncol = 50, nrow = 0 ) )
  colnames(merged) <- c( colnames(mae), colnames(fraser) , colnames(outrider))
  outrider_c <- outrider_sample %>%
    mutate(OUTRIDER_start = as.numeric(OUTRIDER_start)) %>% 
    arrange(desc(OUTRIDER_start)) 
  
  # empty df for merging convenience
  empty.mae <- setNames(data.frame(matrix(ncol = 17, nrow = 0)), colnames(mae))
  empty.fraser <- setNames(data.frame(matrix(ncol = 17, nrow = 0)), colnames(fraser))
  empty.outrider <- setNames(data.frame(matrix(ncol = 16, nrow = 0)), colnames(outrider))
  
  # loop all chromosomes
  for ( c in c(1:length(chromosomes))){
    # subset chromosome c 
    fraser_c<- filter(fraser_sample, FRASER_seqnames==chromosomes[c]) %>% arrange(desc(FRASER_start))
    mae_c <- filter(mae_sample, MAE_contig==chromosomes[c]) %>% arrange(MAE_position)
    
    # count empty sets
    empty.cnt <- (nrow(fraser_c) == 0)+ (nrow(mae_c) == 0) + (nrow(outrider_c) == 0)
    
    if (empty.cnt == 3){
      next
      
    } else if (empty.cnt == 2){ 
      merged <- merged %>% 
        rbind( cbindX(mae_c, fraser_c, outrider_c) ) %>% # merge
        replace(is.na(.), ".") # NAs to "."
      
    } else if (empty.cnt == 1){
      if (nrow(outrider_c) == 0){
        # create empty df to store merged data
        full.join <- data.frame( matrix( ncol = 34, nrow = 0 ) )
        colnames(full.join) <- c( colnames(mae), colnames(fraser))
        
        for (i in c(1:nrow(mae_c))){
          pos <- mae_c$MAE_position[i]
          for (j in c(1:nrow(fraser_c))){
            # if pos between outrider start & end
            if (pos >= fraser_c$FRASER_start[j] & pos <= fraser_c$FRASER_end[j]){
              full.join <- rbind(full.join, cbindX(mae_c[i,], fraser_c[j,]))
            } else if (pos < fraser_c$FRASER_start[j])
            {
              break
            }
          }
        }
        
        #identify un-joined rows
        fraser_anti <- anti_join(fraser_c, full.join, by=character())
        mae_anti <- anti_join(mae_c, full.join, by=character())
        
        merged <- merged %>%
          rbind( cbindX(mae_anti, empty.fraser, empty.outrider), 
                 cbindX(empty.mae, fraser_anti, empty.outrider),
                 cbindX(full.join, empty.outrider)) %>% # merge
          replace(is.na(.), ".") # NAs to "."
        
      } else if (nrow(fraser_c) == 0){
        
        # create empty df to store merged data
        full.join <- data.frame( matrix( ncol = 33, nrow = 0 ) )
        colnames(full.join) <- c( colnames(mae), colnames(outrider))
        
        for (i in c(1:nrow(mae_c))){
          pos <- mae_c$MAE_position[i]
          for (j in c(1:nrow(outrider_c))){
            # if pos between outrider start & end
            if (pos >= outrider_c$OUTRIDER_start[j] & pos <= outrider_c$OUTRIDER_end[j]){
              full.join <- rbind(full.join, cbindX(mae_c[i,], outrider_c[j,]))
            } else if (pos < outrider_c$OUTRIDER_start[j] | is.na(outrider_c$OUTRIDER_start[j]))
            {
              break
            }
          }
        }
        
        #identify un-joined rows
        outrider_anti <- anti_join(outrider_c, full.join, by="OUTRIDER_start")
        mae_anti <- anti_join(mae_c, full.join, by="MAE_position")
        
        merged <- merged %>%
          rbind( cbindX(mae_anti, empty.fraser, empty.outrider), 
                 cbindX(empty.mae, empty.fraser, outrider_anti),
                 cbindX(full.join[,1:17], empty.fraser, full.join[,18:33])) %>% # merge
          replace(is.na(.), ".") # NAs to "."
        
      } else if (nrow(mae_c) == 0){
        
        # create empty df to store merged data
        full.join <- data.frame( matrix( ncol = 33, nrow = 0 ) )
        colnames(full.join) <- c( colnames(fraser), colnames(outrider))
        
        for (i in c(1:nrow(fraser_c))){
          start <- fraser_c$FRASER_start[i]
          end <- fraser_c$FRASER_end[i]
          for (j in c(1:nrow(outrider_c))){
            # if fraser start & end between outrider start & end
            if (start >= outrider_c$OUTRIDER_start[j] & end <= outrider_c$OUTRIDER_end[j]){
              full.join <- rbind(full.join, cbindX(fraser_c[i,], outrider_c[j,]))
            } else if (start < outrider_c$OUTRIDER_start[j] | is.na(outrider_c$OUTRIDER_start[j]))
            {
              break
            }
          }
        }
        
        #identify un-joined rows
        fraser_anti <- anti_join(fraser_c, full.join, by="FRASER_start")
        outrider_anti <- anti_join(outrider_c, full.join, by="OUTRIDER_start")
        
        merged <- merged %>%
          rbind( cbindX(empty.mae, empty.fraser, outrider_anti), 
                 cbindX(empty.mae, fraser_anti, empty.outrider),
                 cbindX(empty.mae, full.join)) %>% # merge
          replace(is.na(.), ".") # NAs to "."
      }
      
    } else if (empty.cnt == 0){
      # join mae & fraser first
      # create empty df to store merged data
      full.join <- data.frame( matrix( ncol = 34, nrow = 0 ) )
      colnames(full.join) <- c( colnames(mae), colnames(fraser))
      
      for (i in c(1:nrow(mae_c))){
        pos <- mae_c$MAE_position[i]
        for (j in c(1:nrow(fraser_c))){
          # if pos between outrider start & end
          if (pos >= fraser_c$FRASER_start[j] & pos <= fraser_c$FRASER_end[j]){
            full.join <- rbind(full.join, cbindX(mae_c[i,], fraser_c[j,]))
          } else if (pos < fraser_c$FRASER_start[j])
          {
            break
          }
        }
      }
      
      #identify un-joined rows
      fraser_anti <- anti_join(fraser_c, full.join, by=character())
      mae_anti <- anti_join(mae_c, full.join, by=character())
      # merge
      full.join <- rbind(full.join, cbindX(mae_anti, empty.fraser), cbindX(empty.mae, fraser_anti))
      
      # join outrider
      # create empty df to store merged data
      full.join.outrider <- data.frame( matrix( ncol = 50, nrow = 0 ) )
      colnames(full.join.outrider) <- c( colnames(full.join), colnames(outrider))
      
      for (i in c(1:nrow(full.join))){
        if (!is.na( full.join$MAE_position[i])){
          pos <- full.join$MAE_position[i]
          for (j in c(1:nrow(outrider_c))){
            # if pos between outrider start & end
            if (pos >= outrider_c$OUTRIDER_start[j] & pos <= outrider_c$OUTRIDER_end[j]){
              full.join.outrider <- rbind(full.join.outrider, cbindX(full.join[i,], outrider_c[j,]))
            } else if (pos < outrider_c$OUTRIDER_start[j] | is.na(outrider_c$OUTRIDER_start[j]))
            {
              break
            }
          }
        } else 
        {
          start <- full.join$FRASER_start[i]
          end <- full.join$FRASER_end[i]
          for (j in c(1:nrow(outrider_c))){
            # if fraser start & end between outrider start & end
            if (start >= outrider_c$OUTRIDER_start[j] & end <= outrider_c$OUTRIDER_end[j]){
              full.join.outrider <- rbind(full.join.outrider, cbindX(full.join[i,], outrider_c[j,]))
            } else if (start < outrider_c$OUTRIDER_start[j] | is.na(outrider_c$OUTRIDER_start[j]))
            {
              break
            }
          }
        }
      }
      
      #identify un-joined rows
      outrider_anti <- anti_join(outrider_c, full.join.outrider, by="OUTRIDER_start")
      full.join_anti <- anti_join(full.join, full.join.outrider, by=c("MAE_position", "FRASER_start"))
      
      merged <- merged %>%
        rbind( cbindX(full.join_anti, empty.outrider), 
               cbindX(empty.mae, empty.fraser, outrider_anti),
               full.join.outrider) %>% # merge
        replace(is.na(.), ".") # NAs to "."
    }
  }
  
  return(merged)
}



# ==== Input FRASER ====
fraser <- read.csv('result_tables/fraser-result.csv') %>% #read data table
          subset(select = -c(X)) %>% #remove first column
          rename_all( function(colname) paste0("FRASER_", colname)) %>% #add prefix
          filter(grepl("^P", FRASER_sampleID)) %>% #keep patient data
          arrange(FRASER_sampleID, FRASER_seqnames, FRASER_start) %>% #reorder
          replace_na(list(FRASER_hgncSymbol=".")) #replace NAs with "."

# ==== Input MAE ====
mae <- read.csv('result_tables/1126MAEanalysistest.raw.csv') %>% #read data table
      subset(select = -c(X)) %>% #remove first column
      rename_all( function(colname) paste0("MAE_", colname)) %>% #add prefix
      rename(MAE_sampleID = MAE_MAE_ID) %>% #rename column
      mutate(MAE_sampleID = sub( '0014A', 'P014', MAE_sampleID)) %>% #rename row (TODO:dynamic)
      mutate(MAE_sampleID = sub( '0026A', 'P026', MAE_sampleID)) %>% #rename row (TODO:dynamic)
      arrange(MAE_sampleID, MAE_contig, MAE_position) #reorder

# ==== Input OUTRIDER ====
outrider <- read.csv('result_tables/OUTRIDER_result.csv') %>% #read data table
            rename_all( function(colname) paste0("OUTRIDER_", colname)) %>% #add prefix
            # subset(OUTRIDER_aberrant == TRUE) %>%
            mutate(OUTRIDER_sampleID = str_replace(OUTRIDER_sampleID, 'RNA-Control-', 'C')) %>% #rename sample
            mutate(OUTRIDER_sampleID = str_replace(OUTRIDER_sampleID, 'RNA-', 'P')) %>% #rename sample
            filter(grepl("^P", OUTRIDER_sampleID)) %>% #keep patient data
            mutate(OUTRIDER_sampleID = str_sub(OUTRIDER_sampleID, end=-2)) #rename sample

# define dataset for biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 
# query range with gene name
outrider.gene.name <- filter(outrider, !grepl("^ENSG", OUTRIDER_geneID))
queryResult <- getBM(attributes=c('start_position', 'end_position', 'external_gene_name'),
                      filters = c('external_gene_name'), 
                      values = outrider.gene.name$OUTRIDER_geneID, mart=ensembl) %>% #query
              rename(OUTRIDER_geneID = external_gene_name) %>% #rename column
              unique()
outrider.gene.name <- merge(outrider.gene.name, queryResult, by ='OUTRIDER_geneID', all=T) %>%
                      rename(OUTRIDER_start = start_position) %>% #rename column
                      rename(OUTRIDER_end = end_position) %>% #rename column
                      relocate(OUTRIDER_start, OUTRIDER_end, .after=OUTRIDER_geneID) #reorder column
# query range with ensembl
outrider.ensg <- filter(outrider, grepl("^ENSG", OUTRIDER_geneID)) %>%
                 mutate(OUTRIDER_geneID = sub("\\..*$", "", OUTRIDER_geneID))
queryResult <- getBM(attributes=c('ensembl_gene_id', 'start_position', 'end_position', 'external_gene_name'),
                  filters = c('ensembl_gene_id'), 
                  values = outrider.ensg$OUTRIDER_geneID, mart=ensembl) %>% #query
              rename(OUTRIDER_geneID = ensembl_gene_id ) %>% #rename column
              unique()
outrider.ensg <- merge(outrider.ensg, queryResult, by ='OUTRIDER_geneID', all=T) %>%
                rename(OUTRIDER_start = start_position) %>% #rename column
                rename(OUTRIDER_end = end_position) %>% #rename column
                relocate(OUTRIDER_start, OUTRIDER_end, .after=OUTRIDER_geneID) %>% #reorder column
                mutate(OUTRIDER_geneID = ifelse(is.na(external_gene_name), OUTRIDER_geneID, external_gene_name)) %>% #replace ensembl if gene name exists
                subset(select = -c(external_gene_name)) #remove gene name column
# merge gene.name & ensg
outrider <- rbind(outrider.gene.name, outrider.ensg) %>% #merge
            arrange(OUTRIDER_sampleID, OUTRIDER_start, OUTRIDER_end) %>% #reorder
            replace_na(list(OUTRIDER_start=".", OUTRIDER_end=".")) #replace NAs with "."
# remove unneeded variables
rm(outrider.ensg, queryResult, outrider.gene.name)


# ==== Main ====
# create chromosome/contig list
chromosomes <- sort(unique( append( unique(fraser$FRASER_seqnames), unique(mae$MAE_contig))))
# create sample list
samples <- unique( c(unique(fraser$FRASER_sampleID), unique(mae$MAE_sampleID), unique(outrider$OUTRIDER_sampleID)))
# create directory to store merged results
dir.create("./merged_results", showWarnings = FALSE)

# loop all samples
for (s in c(1:length(samples))){
  if (s != 14 | s!= 26){
    next
  }
  print(s)
  
  # subset sample s
  fraser_sample <- filter(fraser, FRASER_sampleID==samples[s])
  mae_sample <- filter(mae, MAE_sampleID==samples[s])
  outrider_sample <- filter(outrider, OUTRIDER_sampleID==samples[s])
  
  # merge
  merged_result <- mergeAll(fraser_sample, mae_sample, outrider_sample)
  
  # Post processing of merged_result
  merged_result <- merged_result %>%
                    replace(is.na(.), ".") %>% # NAs to "."
                    mutate(sampleID = ifelse(MAE_sampleID == ".",
                                             ifelse(FRASER_sampleID == ".", OUTRIDER_sampleID, FRASER_sampleID),
                                             MAE_sampleID)) %>%
                    subset(select = -c(MAE_sampleID, FRASER_sampleID, OUTRIDER_sampleID)) %>%
                    mutate(seqnames = ifelse(MAE_contig == ".", FRASER_seqnames, MAE_contig)) %>%
                    subset(select = -c(MAE_contig, FRASER_seqnames)) %>%
                    relocate(sampleID, seqnames) %>%
                    arrange(sampleID, seqnames, OUTRIDER_start, FRASER_start, MAE_position) 
                    
  # write output
  write_tsv(merged_result, paste0('./merged_results/', samples[s], '.tsv'))
}
