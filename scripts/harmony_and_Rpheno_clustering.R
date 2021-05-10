#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the clustering with Rphenograph.
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
  	require(CATALYST, quietly=T)
  	require(SingleCellExperiment, quietly=T)
  	require(Rphenograph, quietly=T)
  	require(dplyr, quietly=T)
  	require(cowplot, quietly=T)
  	require(gridExtra, quietly=T)
  	require(ggplot2, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile_pref"), type="character", default=NULL, help="File that contains the sce object for which to do the clustering", metavar="character"),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--datatransf"), type="character", default="exprs",  help="Data transformation for running the clustering on. Options: counts=no transformation, exprs= arcsinh 
		transformation, scaled=scaled transformation of the exprs, scaledtrim=the scaled and trimmed values (q=0.01) of the exprs. Can include more than one option 
		in comma separated style, e.g. exprs,scaled. (default=exprs)"),
    make_option(c("--outputtable"), type="logical", action="store_true", help="Include flag to create the Zegami table."),
    make_option(c("--outdirZegami"), type="character", default=NULL,  help="output directory for the Zegami table (to be same as input from preprocessing)", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

source("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/github_scripts/scripts/plotting_functions.R")

###########################
# loading the files and running the analysis
cat("Starting the clustering with Rphenograph... \n")
datatransf <- unlist(strsplit(opt$datatransf, split=","))

for (i in datatransf) {
	cat("Running the Harmony and Rphenograph for ",i," \n")	
	cat("Loading",paste0(opt$infile_pref,"_",i,".RData"), " \n")	
	load(paste0(opt$infile_pref,"_",i,".RData"))
	sce_pheno <- data.frame(reducedDims(sce)[[paste0("TSNE_",i)]], reducedDims(sce)[[paste0("UMAP_",i)]])
	colnames(sce_pheno) <- paste0(c("TSNE1","TSNE2", "UMAP1","UMAP2"),"_",i)

	exp <- assay(sce, i)[rowData(sce)$marker_name[rowData(sce)$marker_class=="type"],]
	my_harmony_embeddings <- HarmonyMatrix(
  		data_mat  = exp,
  		meta_data = colData(sce)$sample,
		do_pca    = TRUE)

	R_pheno_out <- Rphenograph(my_harmony_embeddings, k=30)
	cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
	sce[[paste0("harmony_phenograph_cluster_",i,"_k30")]] <- factor(membership(R_pheno_out[[2]]))

	sce_pheno <- uwot::umap(my_harmony_embeddings,  n_components = 2)
	colnames(sce_pheno) <- c("UMAP1","UMAP2")
	sce_pheno <- cbind(sce_pheno, sce@colData)
	p2a <- ggplot(as.data.frame(sce_pheno), aes(UMAP1, UMAP2, col=annotation)) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c26)
	p2b <- ggplot(as.data.frame(sce_pheno), aes(UMAP1, UMAP2, col=annotation_all)) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c26)
	p2c <- ggplot(as.data.frame(sce_pheno), aes(UMAP1, UMAP2, col=sample_id)) + geom_point(size = 0.5)+theme_bw()


