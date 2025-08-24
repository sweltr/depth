#!/usr/bin/env Rscript

getwd()


########################################
#
# Setup environment 2_filt_err_derep
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
# Filter 
#
load("../CUTADAPT/rdata/1.2_cutadapt_ssu.rdata")
########################################


## ----set_paths_and_names--------------------
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
filt_path <- "/"

## ----filter_path, eval=FALSE---------------------------------
filtFs <- file.path("filtered", basename(cutFs))
filtRs <- file.path("filtered", basename(cutRs))


## ----filter, eval=FALSE--------------------------------------
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,  
                     maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, rm.phix = TRUE, 
                     compress = TRUE, multithread = 20) 
out


out <- data.frame(out)
out <- tibble::rownames_to_column(out, "SampleID")

out <- out %>% mutate(SampleID = str_replace_all(SampleID, "_R[0-9].*", ""))

readr::write_delim(out, "ssu_filter_track.txt", delim = "\t")

## -------------------------------------------
p1 <- plotQualityProfile(filtFs[1:41], aggregate = TRUE)
p2 <- plotQualityProfile(filtRs[1:41], aggregate = TRUE)

p3 <- grid.arrange(p1, p2, nrow = 1)
ggsave("figures/ssu_post_filt_plot_qscores.png", p3, width = 7, height = 3)

save.image("rdata/2.1_filt_ssu.rdata")

########################################
#
# Error 
#
########################################

## ----learn_errors_forward, eval=FALSE------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
## ----learn_errors_reverse, eval=FALSE------------------------
errR <- learnErrors(filtRs, multithread = TRUE)

## ----plot_errF------------------------------
p3 <- plotErrors(errF, nominalQ = TRUE)
ggsave("figures/plot_errorF_1.png", p3, width = 7, height = 5)
ggsave("figures/plot_errorF_2.png", p3)
## ----plot_errR------------------------------
p4 <- plotErrors(errR, nominalQ = TRUE)
ggsave("figures/plot_errorR_1.png", p4, width = 7, height = 5)
ggsave("figures/plot_errorR_2.png", p4)


########################################
#
# DEREP 
#
########################################

## ----derep----------------------------------

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

save.image("rdata/2.2_err_derep_ssu.rdata")

quit()