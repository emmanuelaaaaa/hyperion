#!/usr/bin/env Rscript
##################
## Description 
##
## This script filters the sce object.
##
###########################
## libraries 

suppressPackageStartupMessages({
  require(optparse, quietly=T)
  require(SingleCellExperiment, quietly=T)
})

source("/t1-data/user/erepapi/Fellowship/COVID19_CyTOF/scripts/other_plots.R")

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the clustering", metavar="character"),
    make_option(c("--filters"), type="character", default=NULL, help="Comma separated filters to be evaluated. The variables need to be column names. (e.g. Area>30,Eccentricity<=0.4)", metavar="character"),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

dir.create(file.path(opt$outdir, "figures", opt$label), showWarnings = FALSE)

########### 
## loading the sce object from the initial QC analysis - no need to preprocess again
filters <- gsub(" ", replacement="", opt$filters)
filters <- unlist(strsplit(filters, split=","))

cat("Loading the sce obj...\n")
load(opt$infile)

for (i in filters) {
	var <- strsplit(i, split="<|>|<=|>=")[[1]][1]
	num <- gsub("[A-Za-z]", "", i)

	filter <- paste0("colData(sce)[['",var,"']]",num)
	index <- eval(parse(text=filter))
	cat("filtering for ", filter, " : removing ", length(which(index)), "numbers of cells. \n")
	sce <- sce[,index]
}

cat("The final sce object has dimensions:", dim(sce), "\n")

save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, "_", i, "filt.RData")))