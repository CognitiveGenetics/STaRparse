#!/usr/bin/env Rscript

# IMPORT LIBRARIES
options(warn=-1)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(optparse))
options(warn=0)
# DEFINE ARGUMENTS FOR SCRIPT
option_list <- list(
  make_option(c("-d", "--dir"), type="character", action="store", default=NA, help="Script directories", metavar="character"),
  make_option(c("-r", "--input_reads"), type="character", action="store", default=NA, help="Read file path", metavar="character"),
  make_option(c("-c", "--input_coverage"), type="character", action="store", default=NA, help="Coverage file path", metavar="character"),
  make_option(c("-o", "--out"), type="character", action="store", default=NA, help="output file path and name", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ARGUMENT CHECK
if (length(opt) != 5){
  stop("Please ensure all aruments are supplied", call.=FALSE)
  } else {
	# DEFINE VARIABLES
	reads <- read.csv(file=opt$input_reads, sep="\t", header = TRUE)
	coverage <- read.csv(file=opt$input_coverage, sep="\t", header = TRUE)                                                                                   
	# RUN FILTER
	source(paste0(opt$dir, "/Source/filter/filter_function.R"))
        outputdf <- filter_by_coverage(reads, coverage)
	# REMOVE LARGE FILES
        rm(reads, coverage)
	# SAVE FILTERED CSV
        outputdf$Call_ID <- str_replace(outputdf$Call_ID, "_ExpansionHunter", "")
	write.table(outputdf, opt$out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
	}	


