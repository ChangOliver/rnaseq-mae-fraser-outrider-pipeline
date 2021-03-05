pkgs <- c('BiocManager', 'plyr', 'tidyverse', 'stringr', 'gdata', 'dplyr', 'data.table',
			'tibble', 'ggplot2', 'ggrepel', 'httr', 'jsonlite', 'rlist', 'magrittr', 'optparse')

install.packages(pkgs)

BiocManager::install("biomaRt")
BiocManager::install("OUTRIDER")
BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
BiocManager::install("FRASER")

install.packages("remotes")
remotes::install_github("mumichae/tMAE")
