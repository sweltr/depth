#!/usr/bin/env Rscript

getwd()

########################################
#
# Setup environment 5_tax_assign
#
########################################
set.seed(919191)
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(grid)
library(DECIPHER); packageVersion("DECIPHER")
library(magrittr)
library(dplyr)
library(stringr)

########################################
#
# TAXONOMY
#
#load("rdata/4.1_pre_tax_ssu_pseudo.rdata")
########################################

## ----read_combo, eval=FALSE-----------------
remove(list = ls())
seqtab.trim.nochim.consensus <- readRDS("rdata/seqtab.nochim_pseudo.rds")

########################################
# TAXONOMY = silva
########################################
seqtab.consensus <- seqtab.trim.nochim.consensus

tax_silva_v138.consensus <- assignTaxonomy(seqtab.consensus, "TAXONOMY_FILES/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v138.consensus, "5.tax_silva_v138.consensus.rds")

tax_silva_v132.consensus <- assignTaxonomy(seqtab.consensus, "TAXONOMY_FILES/silva_nr_v132_train_set.fa.gz", multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v132.consensus, "5.tax_silva_v132.consensus.rds")

########################################
# TAXONOMY = RDP
########################################

tax_rdp_v138.consensus <- assignTaxonomy(seqtab.consensus, "TAXONOMY_FILES/rdp_train_set_18.fa.gz", multithread = TRUE, verbose = TRUE)
saveRDS(tax_rdp_v138.consensus, "5.tax_rdp_v138.consensus.rds")

########################################
# TAXONOMY = ITGDB
########################################

tax_itgdb.consensus <- assignTaxonomy(seqtab.consensus, "TAXONOMY_FILES/itgdb_dada2.fa", multithread = TRUE, verbose = TRUE)
saveRDS(tax_itgdb.consensus, "5.tax_itgdb.consensus.rds")

########################################
# TAXONOMY = GSRDB
########################################

tax_gsrdb.consensus <- assignTaxonomy(seqtab.consensus, "TAXONOMY_FILES/gsrdb_dada2.fa", multithread = TRUE, verbose = TRUE)
saveRDS(tax_gsrdb.consensus, "5.tax_gsrdb.consensus.rds")

## ----save_image, eval=FALSE----------------------------------
save.image("rdata/5.1_dada2_wf_ssu_pseudo.rdata")


sessionInfo()
devtools::session_info()

quit()