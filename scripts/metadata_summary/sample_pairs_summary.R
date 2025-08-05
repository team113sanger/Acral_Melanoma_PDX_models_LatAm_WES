#!/usr/bin/env Rscript

library(here)
library(dplyr)
library(readr)
library(tidytable)

projdir<-here()
metadata_dir<-file.path(projdir,"metadata")
biosamp_manif_file<-file.path(metadata_dir,"6633-_2729_3248_sample_pair-manifest_complete_final_calling_list.tsv")

samp_pair_tab<- readr::read_delim(biosamp_manif_file,col_names = TRUE )

# Get the number of sample pair that were submitted for analysis
samp_pair_tab <- sample