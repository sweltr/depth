#!/usr/bin/env Rscript

getwd()

########################################
#
# Setup environment 3_dada
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
# DADA 
#
load("rdata/2.2_err_derep_ssu.rdata")
########################################

## ----dadaF----------------------------------
dadaFs <- dada(derepFs, err = errF, pool = "pseudo", multithread = 20)

## ----inspect_f, eval=FALSE------------------
dadaFs[[2]]

## ----dadaR----------------------------------
dadaRs <- dada(derepRs, err = errR, pool = "pseudo", multithread = 20)


## ----inspect_r, eval=FALSE------------------
dadaRs[[2]]

save.image("rdata/3.1_dada_ssu_pseudo.rdata")

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), 
               sapply(dadaRs, getN))
colnames(track) <- c("denoisedF", "denoisedR")

track

track <- data.frame(track)
track <- tibble::rownames_to_column(track, "SampleID")

track <- track %>% mutate(SampleID = str_replace_all(SampleID, "_R.*", ""))

readr::write_delim(track, "ssu_dada_track.txt", delim = "\t")

save.image("rdata/3.1_dada_ssu_pseudo.rdata")

