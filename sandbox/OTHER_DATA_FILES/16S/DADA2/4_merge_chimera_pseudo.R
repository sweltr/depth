#!/usr/bin/env Rscript

getwd()

########################################
#
# Setup environment 4_merge_chimera
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
# MERGE & SEQ TABLE 
#
load("rdata/3.1_dada_ssu_pseudo.rdata")
########################################

## ----merge_paired_reads, eval=FALSE---------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

## ----head_file, eval = FALSE, echo = FALSE----
head(mergers[[1]])

## ----seq_table, eval=FALSE------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

read_length <-  data.frame(nchar(getSequences(seqtab)))
colnames(read_length) <- "length"
plot_a <- qplot(length, data = read_length, geom = "histogram", binwidth = 1, xlab = "read length", ylab = "total variants", xlim = c(240, 280)) 
ggsave("figures/ssu_read_length_before_pseudo.png", plot_a, width = 7, height = 3)
saveRDS(seqtab, "rdata/seqtab_merge_pseudo.rds")

seqtab.2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(251,255)]
dim(seqtab.2)
table(nchar(getSequences(seqtab.2)))

saveRDS(seqtab.2, "rdata/seqtab_pseudo.rds")

########################################
#
# COLLAPSE
# COMMENT chunk to skip COLLAPSE
########################################

#seqtab_collapse <- collapseNoMismatch(seqtab.2, minOverlap = 20, orderBy = "abundance",
#  identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)

#dim(seqtab_collapse)
#table(nchar(getSequences(seqtab_collapse)))

#read_length_all_collapse <-  data.frame(nchar(getSequences(seqtab_collapse)))
#colnames(read_length_all_collapse) <- "length"
#plot_all_collapse <- qplot(length, data = read_length_all_collapse, geom = "histogram", binwidth = 1, xlab = "read length", ylab = "total variants", xlim = c(240, 280)) 
#ggsave("figures/read_length_before_pseudo_all_collapse.png", plot_all_collapse, width = 7, height = 3)
#saveRDS(seqtab_collapse, "rdata/seqtab_collapse_pseudo.rds")
#seqtab.2 <- seqtab_collapse


########################################
#
# CHIMERA
#
########################################

## ----chimera-------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab.2,  method = "consensus", multithread = 20, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.2)

saveRDS(seqtab.nochim, "rdata/seqtab.nochim_pseudo.rds")

## ----seq_table2, eval=FALSE-----------------
table(nchar(getSequences(seqtab.nochim)))

table(colSums(seqtab.nochim>0))
table(rowSums(seqtab.nochim>0))


## ----save_image_pre_tax, eval=FALSE---------
save.image("rdata/4.1_pre_tax_ssu_pseudo.rdata")

## ----track_reads----------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(fnFs, getN), 
               sapply(fnFs.filtN, getN), 
               out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("raw", "pre_filt", "cut", "filtered", "denoisedF", 
                     "denoisedR", "merged",  "nonchim")

track <- data.frame(track)
track <- tibble::rownames_to_column(track, "SampleID")

track <- track %>% mutate(SampleID = str_replace_all(SampleID, "^.*/", "")) %>%
  mutate(SampleID = str_replace_all(SampleID, "_R[0-9].*", ""))


readr::write_delim(track, "ssu_wf_read_changes.txt", delim = "\t")

#rownames(track) <- sample.names
#track

## ----track_changes-------------------------------------------------------
#write.table(track, "tables/ssu_read_changes_pseudo.txt", sep = "\t", quote = FALSE, col.names=NA)

## ----save_image_pre_tax, eval=FALSE---------
save.image("rdata/4.1_pre_tax_ssu_pseudo.rdata")

quit()