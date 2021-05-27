#!/usr/bin/env Rscript
##################
## Description 
##
## This script extracts the region of interest (or reads in the IDS of the cells) and give the proportions of cells in each group of the annotations. 
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
  	require(SingleCellExperiment, quietly=T)
  	require(dplyr, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to extract the region and calculate the freqs.", metavar="character"),
    make_option(c("--outdir_table"), type="character", default=NULL,  help="output directory for figures (in datashare)", metavar="character"),
    make_option(c("--annot"), type="character", default=NULL,  help="File with annotations to be loaded over which to be calculating the freqs", metavar="character"),
    make_option(c("--ids"), type="character", default=NULL,  help="File with IDs of cells to be calculating the freqs (to be used instead of coordinates).", metavar="character"),
    make_option(c("--xlim"), type="character", default=NULL,  help="X coordinates (comma separated) over which to be calculating the freqs", metavar="character"),
    make_option(c("--ylim"), type="character", default=NULL,  help="Y coordinates (comma separated) over which to be calculating the freqs", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$annot),is.null(opt$xlim) & is.null(opt$ylim) & is.null(opt$ids),is.null(opt$outdir_table),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

###########################
# loading the files and returning the freqs
cat("Loading ",opt$infile, " \n")	
load(opt$infile)

annot <- read.table(opt$annot, header=T, stringsAsFactors=F, sep="\t")
annot <- annot[match(colData(sce)$cellID_name,annot$cellID_name),]
if(any(colData(sce)$cellID_name!=annot$cellID_name)) {
	stop("Attention: there are cells that have not been annotated properly")
}
colData(sce)$annotation <- annot$annotation

if (!is.null(opt$xlim)) {
	xlim <- unlist(strsplit(opt$xlim, split=","))
	cat("Subsetting the dataset for X>=", xlim[1], " and X<=", xlim[2], "\n")
	sce <- sce[,colData(sce)$LocX>=as.integer(xlim[1]) & colData(sce)$LocX<=as.integer(xlim[2])]
}

if (!is.null(opt$ylim)) {
	ylim <- unlist(strsplit(opt$ylim, split=","))
	cat("Subsetting the dataset for Y>=", ylim[1], " and Y<=", ylim[2], "\n")
	sce <- sce[,colData(sce)$LocY>=as.integer(ylim[1]) & colData(sce)$LocY<=as.integer(ylim[2])]
}


if (!is.null(opt$ids)) {
	ids <- read.table(opt$ids, header=F, stringsAsFactors=F, sep="\t")
	cat("Subsetting the dataset for the IDs in the ids file. \n")
	sce <- sce[,ids]
}

temp <- as.data.frame(table(colData(sce)$annotation)) %>% mutate(perc = Freq / sum(Freq)*100)
write.table(temp, file=file.path(opt$outdir_table, paste0(opt$label, "_cellfreqs.txt")), row.names=F, quote=F,sep="\t")



