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

if (any(is.null(opt$infile_pref),is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

source("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/github_scripts/scripts/plotting_functions.R")

###########################
# loading the files and running the analysis
cat("Starting the clustering with Rphenograph... \n")
datatransf <- unlist(strsplit(opt$datatransf, split=","))

ks <- c(15,30,50)
for (i in datatransf) {
	cat("Running the Rphenograph for ",i," \n")	
	cat("Loading",paste0(opt$infile_pref,"_",i,".RData"), " \n")	
	load(paste0(opt$infile_pref,"_",i,".RData"))
	sce_pheno <- data.frame(reducedDims(sce)[[paste0("TSNE_",i)]], reducedDims(sce)[[paste0("UMAP_",i)]])
	colnames(sce_pheno) <- paste0(c("TSNE1","TSNE2", "UMAP1","UMAP2"),"_",i)
	for (k in ks) { 

		R_pheno_out <- Rphenograph(t(assay(sce, i))[,rowData(sce)$marker_name[rowData(sce)$marker_class=="type"]], k=k)

		# Modularity is one measure of the structure of networks or graphs. It was designed to measure the strength of division of a network into modules (clusters/communities).
		# Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules. Modularity is often used in optimization methods
		# for detecting community structure in networks.
		# It has been shown that modularity suffers a resolution limit and, therefore, it is unable to detect small communities.
		cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")

		sce[[paste0("phenograph_cluster_",i,"_k", k)]] <- factor(membership(R_pheno_out[[2]]))

		one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id=paste0("phenograph_cluster_",i,"_k", k),
			features=rowData(sce)$marker_name[rowData(sce)$marker_class=="type"], row_anno = TRUE, bars=T) 
		cat("Plotting the heatmap plots of the  phenograph clusters... \n")
		pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_clustering_plots_heatmaps_",i,"_k", k,".pdf")), height=10, width=15)
		print(one_plot_heatmap)
		dev.off()

		one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id=paste0("phenograph_cluster_",i,"_k", k), features = rowData(sce)$marker_name[rowData(sce)$marker_class=="type"]) 
		cat("Plotting the expression density plots of the phenograph clusters... \n")
		pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_clustering_plots_exprsdens_",i,"_k", k,".pdf")), height=10, width=15)
		print(one_plot_exprs)
		dev.off()

		assign(paste0("R_pheno_",i,"_k", k), R_pheno_out)
	}

	pheno_ks <- paste0("phenograph_cluster_",i,"_k", ks)
	sce_pheno <- cbind(sce_pheno, sce@colData)
	one_plot_tsne <- function(k)  ggplot(sce_pheno, aes_string(x=paste0("TSNE1_",i), y=paste0("TSNE2_",i), col=k)) + geom_point(size = 1)+theme_bw()
	one_plot_umap <- function(k)  ggplot(sce_pheno, aes_string(x=paste0("UMAP1_",i), y=paste0("UMAP2_",i), col=k)) + geom_point(size = 1)+theme_bw()
	meta_plots_t <- plyr::llply(pheno_ks, one_plot_tsne)
	meta_plots_u <- plyr::llply(pheno_ks, one_plot_umap) 

	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_TSNEUMAP_", opt$label, "_",i,".pdf")), height=15, width=25)
	    print(meta_plots_t)
	    print(meta_plots_u)
	dev.off()
	cat("Updating the sce object in ", file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")),"\n")
	save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, "_", i, ".RData")))

}

###########################
## Creating table for Zegami

if (opt$outputtable) {
	source("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/github_scripts/scripts/zegami_function.R")
	cat("\nMaking the table for Zegami...\n")
	makeZegami(sce, output_dir=opt$outdirZegami, label=opt$label, zegami_suffix="")
}

###########################
# saving
rm(ks, my_mat, other_vars, my_clusters2, Zeg_table)
session <- sessionInfo()
save.image(file=file.path(opt$outdir, "RData", paste0("Rpheno_", opt$label, ".RData")))


