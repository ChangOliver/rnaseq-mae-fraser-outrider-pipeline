pkgs <- c('BiocManager', 'plyr', 'tidyverse', 'gdata', 'stringr', 'dplyr', 'data.table',
			'tibble', 'ggplot2', 'ggrepel', 'optparse', 'remotes')

install.packages(pkgs)

BiocManager::install("biomaRt")
BiocManager::install("OUTRIDER")
BiocManager::install("MafDb.gnomAD.r2.1.hs37d5")
BiocManager::install("MafDb.gnomAD.r2.1.GRCh38")
BiocManager::install("FRASER")

remotes::install_github("mumichae/tMAE")
