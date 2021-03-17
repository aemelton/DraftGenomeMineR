###
#install.packages("Bioconductor")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("ORFik")
#install.packages("devtools")
#library(devtools)
#install_github("mhahsler/rBLAST")
#install.packages("remotes")
#remotes::install_github("GuillemSalazar/FastaUtils")

#
list.of.packages <- c("ape",
                      "Biostrings",
                      "dplyr",
                      "FastaUtils",
                      "ORFik",
                      "readr",
                      "tidyr",
                      "rBLAST",
                      "seqinr",
                      "stringr")
#

#
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
	install.packages(new.packages)
	}
#

#
lapply(list.of.packages, require, character.only = TRUE)
#

#                      
setwd("DraftGenomeMineR/Functions/")
files.sources <- list.files()
sapply(files.sources, source)
#