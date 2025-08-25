#!/usr/bin/env Rscript

getwd()


########################################
#
# Setup environment 1_cut
#
########################################
set.seed(919191)
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(gridExtra)
library(grid)
library(DECIPHER); packageVersion("DECIPHER")
library(magrittr)
library(tidyverse) # contains the following packages: dplyr, readr, forcats, stringr, ggplot2, tibble, lubridate, tidyr, purrr
########################################
#
# PREPROCESSING
#
########################################

## ---- warning=FALSE------------------------------------------
path <- "/pool/genomics/stri_istmobiome/data/SWELTR/RAW_DATA/16S/2019"
#setwd("RAW")
#orig_fastq <- list.files(path = path, pattern = "*.fastq.gz")
#newname_fastq <- gsub("_S.*_L001", "", orig_fastq)
#newname_fastq <- gsub("_001", "", newname_fastq)
#file.rename(orig_fastq, newname_fastq)
#setwd("../")

## ----define path----------------------------
#path <- "RAW"  
head(list.files(path))


## ----sort_pairs-----------------------------
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

p1a <- plotQualityProfile(fnFs[1:41], aggregate = TRUE)
p2a <- plotQualityProfile(fnRs[1:41], aggregate = TRUE)

p3a <- grid.arrange(p1a, p2a, nrow = 1)
ggsave("CUTADAPT/figures/ssu_plot_qscores_raw.png", p3a, width = 7, height = 3)


## ----define_primers, eval=FALSE------------------------------
FWD <- "GTGCCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"


## ----prime_variants-------------------------
allOrients <- function(primer) {
    require(Biostrings)
    dna <- DNAString(primer) 
    orients <- c(Forward = dna, 
                 Complement = complement(dna), 
                 Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

save.image("CUTADAPT/rdata/1.1_setup_ssu.rdata")


########################################
#
# PREFILTER (OPTIONAL)
#
########################################
########## TO PREFILTER UNCOMMENT these Lines ##########

## ----filter_ns------------------------------
#fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
#fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 20)

#p1b <- plotQualityProfile(fnFs.filtN[1:41], aggregate = TRUE)
#p2b <- plotQualityProfile(fnRs.filtN[1:41], aggregate = TRUE)
#p3b <- grid.arrange(p1b, p2b, nrow = 1)
#ggsave("figures/ssu_plot_qscores_pre_filt.png", p3b, width = 7, height = 3)

########## AND COMMENT these Lines ##########
fnFs.filtN <- file.path(path, basename(fnFs)) 
fnRs.filtN <- file.path(path, basename(fnRs))

############################################################
########## END PREFILTERING OPTION #########################
############################################################

########################################
#
# CUTADAPT
#
########################################

## ----primer_assess--------------------------
sampnum <- 2
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

## ----primer_assessF-------------------------
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[sampnum]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[sampnum]]))

## ----primer_assessR-------------------------
rbind(REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[sampnum]]), 
       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[sampnum]]))


## ----cutadapt_setup-------------------------
cutadapt <- "/home/scottjj/miniconda3/envs/cutadapt/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


## ----cutadapt-------------------------------
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

for(i in seq_along(fnFs.filtN)) {system2(cutadapt,
                                   args = c(R1.flags, R2.flags, 
                                            "--times", 2, 
                                            "--minimum-length", 200,
                                            "--maximum-length", 300,
                                            "--error-rate", 0.10, 
                                            "--no-indels",
                                            "--discard-untrimmed",
                                            "--output", fnFs.cut[i], 
                                            "--paired-output", fnRs.cut[i],
                                            "--report full",
                                            "--cores", 10,
                                            fnFs.filtN[i], fnRs.filtN[i]))}

p1_cut <- plotQualityProfile(fnFs.cut[1:41], aggregate = TRUE)
p2_cut <- plotQualityProfile(fnRs.cut[1:41], aggregate = TRUE)

p3_cut <- grid.arrange(p1_cut, p2_cut, nrow = 1)
ggsave("CUTADAPT/figures/ssu_plot_qscores_cut.png", p3_cut, width = 7, height = 3)


## ----check_primers--------------------------
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[sampnum]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[sampnum]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[sampnum]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[sampnum]]))

save.image("CUTADAPT/rdata/1.2_cutadapt_ssu.rdata")

getN_in <- function(x) sum(getUniques(x))
track_in <- cbind( 
               sapply(fnFs, getN_in), 
               sapply(fnRs, getN_in), 
               sapply(fnFs.filtN, getN_in), 
               sapply(fnRs.filtN, getN_in), 
               sapply(fnFs.cut, getN_in),
               sapply(fnRs.cut, getN_in))
colnames(track_in) <- c("in_F", "in_R", "pre_filt_F", "pre_filt_R", "cut_F", "cut_R")

track_in <- data.frame(track_in)
track_in <- tibble::rownames_to_column(track_in, "SampleID")

track_in <- track_in %>% mutate(SampleID = str_replace_all(SampleID, "^.*/", "")) %>%
  mutate(SampleID = str_replace_all(SampleID, "_R.*", ""))

readr::write_delim(track_in, "CUTADAPT/cutadapt_track.txt", delim = "\t")

save.image("CUTADAPT/rdata/1.2_cutadapt_ssu.rdata")

sessionInfo()
devtools::session_info()

quit()
