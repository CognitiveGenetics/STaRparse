#!/usr/bin/env Rscript

#       IMPORT LIBRARIES
options(warn=-1)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
options(warn=0)
#       DEFINE ARGUMENTS FOR SCRIPT
option_list <- list(
  make_option(c("-d", "--dir"), type="character", action="store", default=NA, help="Path of global directory", metavar="character"),
  make_option(c("-r", "--reads"), type="character", action="store", default=NA, help="Read file path", metavar="character"),
  make_option(c("-o", "--out"), type="character", action="store", default=NA, help="output file path and name", metavar="character"),
  make_option(c("-s", "--savename"), type="character", action="store", help="Indicate job name", metavar="character"),
  make_option(c("-b", "--build"), type="character", action="store", help="Specify reference genome build (19, 37, or 38)", metavar="character"),
  make_option(c("-H", "--humandb"), type="character", action="store", help="Path to humandb", metavar="character")
 
)

opt <- parse_args(OptionParser(option_list=option_list))

#       ARGUMENT CHECK
if (length(opt) != 7){
  stop("Please ensure all aruments are supplied", call.=FALSE)
} else {
  
    #IMPORT FUNCTIONS
    source(paste0(opt$dir, "/Source/summaries/summaries_functions.R"))
    
    #IMPORT DATA
    reads <- read.csv(file=opt$reads, sep="\t", header = TRUE)
    # All1 and All2 are repeat length (alleles) observed in the current sample. 
    names(reads) <- c("JSON_Call_ID", "Sample_ID", "Chr", "Start", "End", "GT", "Ref_Units", "All1", "All2")
    outpath <- paste0(opt$out, opt$savename)
    # assign Call_ID as chr.start to
    reads$Call_ID=paste(gsub("chr", "", reads$Chr),reads$Start,sep='.')
 
    #RUN FUNCTIONS
    locusdf <- by_locus(reads, opt$out, opt$savename, opt$build, opt$humandb, opt$dir)
    ratio_of_stability(locusdf, outpath)
    by_chr(reads, locusdf, outpath)
    by_sample(reads, locusdf, outpath)
    
}
