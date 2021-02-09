#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the clustering with Rphenograph.
##
###########################
# libraries 
stopifnot(
  require(optparse, quietly=T), #lib.loc = "/home/staff/erepapi/R/x86_64-pc-linux-gnu-library/3.6"),
  require(CATALYST, quietly=T),
  require(SingleCellExperiment, quietly=T),
  require(Rphenograph, quietly=T),
  require(dplyr, quietly=T),
  require(cowplot, quietly=T),
  require(gridExtra, quietly=T),
  require(ggplot2, quietly=T)
)

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the clustering", metavar="character"),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
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

###########################
# loading the files and running the analysis
cat("Starting the clustering with Rphenograph... \n")
load(opt$infile)

sce_pheno <- data.frame(reducedDims(sce)$TSNE, reducedDims(sce)$UMAP)
colnames(sce_pheno) <- c("TSNE1","TSNE2", "UMAP1","UMAP2")

ks <- c(15,30,50)
for (k in ks) { 
	R_pheno_out_all <- Rphenograph(t(assay(sce)), k=k)
	R_pheno_out_noDNA <- Rphenograph(t(assay(sce))[,setdiff(rownames(sce),c("DNA1","DNA3"))], k=k)
	R_pheno_out_noDNAH3 <- Rphenograph(t(assay(sce))[,setdiff(rownames(sce),c("DNA1","DNA3","H3"))], k=k)

# Modularity is one measure of the structure of networks or graphs. It was designed to measure the strength of division of a network into modules (clusters/communities).
# Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules. Modularity is often used in optimization methods
# for detecting community structure in networks.
# It has been shown that modularity suffers a resolution limit and, therefore, it is unable to detect small communities.
	cat("\n The modularity of the graph is with all markers is ", modularity(R_pheno_out_all[[2]]), "\n")
	cat("The modularity of the graph is without DNA1/3 markers is ", modularity(R_pheno_out_noDNA[[2]]), "\n")
	cat("The modularity of the graph is without DNA1/3,H3 markers is ", modularity(R_pheno_out_noDNAH3[[2]]), "\n")

	sce[[paste0("phenograph_cluster_all_k", k)]] <- sce_pheno[[paste0("phenograph_cluster_all_k", k)]] <- factor(membership(R_pheno_out_all[[2]]))
	sce[[paste0("phenograph_cluster_noDNA_k", k)]] <- sce_pheno[[paste0("phenograph_cluster_noDNA_k", k)]] <- factor(membership(R_pheno_out_noDNA[[2]]))
	sce[[paste0("phenograph_cluster_noDNAH3_k", k)]] <- sce_pheno[[paste0("phenograph_cluster_noDNAH3_k", k)]] <- factor(membership(R_pheno_out_noDNAH3[[2]]))

	assign(paste0("R_pheno_k_", k), R_pheno_out_noDNAH3)
}

###########################
# plots
pheno_ks <- c(paste0("phenograph_cluster_all_k", ks), paste0("phenograph_cluster_noDNA_k", ks),paste0("phenograph_cluster_noDNAH3_k", ks))
one_plot_tsne <- function(k)  ggplot(sce_pheno, aes_string(x="TSNE1", y="TSNE2", col=k)) + geom_point(size = 1)+theme_bw()
one_plot_umap <- function(k)  ggplot(sce_pheno, aes_string(x="UMAP1", y="UMAP2", col=k)) + geom_point(size = 1)+theme_bw()
meta_plots_t <- plyr::llply(pheno_ks, one_plot_tsne)
meta_plots_u <- plyr::llply(pheno_ks, one_plot_umap) 

pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_TSNEUMAP_", opt$label, ".pdf")), height=15, width=25)
    grid.arrange(grobs=meta_plots_t)
    grid.arrange(grobs=meta_plots_u)
dev.off()

###########################
## Creating table for Zegami
paste_matrix <- function(a,mat,sep = "",collapse = NULL){
	matnames <- colnames(mat)
	n <- nrow(mat)
	p <- ncol(mat)
	mat <- matrix(paste(a,mat,sep = sep,collapse = collapse),n,p)
	colnames(mat) <- matnames
	return(mat)
}

if (opt$outputtable) {
	cat("\nMaking the table for Zegami...\n")

	my_mat <- t(assay(sce))
	my_dims <- cbind(reducedDims(sce)$TSNE, reducedDims(sce)$UMAP)
	colnames(my_dims ) <- c("TSNE1","TSNE2","UMAP1","UMAP2")
	other_vars <- sce@colData
	other_vars[,grepl("phenograph", names(other_vars))] <- paste_matrix("cl", as.matrix(other_vars[,grepl("phenograph", names(other_vars))]))
	my_clusters2 <- sapply(sce@metadata$cluster_codes, function(x) x[cluster_ids(sce)])
	my_clusters2 <- paste_matrix("cl", my_clusters2)

	Zeg_table <- cbind(other_vars[!is.na(my_dims[,1]) | !is.na(my_dims[,3]),],
	    my_clusters2[!is.na(my_dims[,1]) | !is.na(my_dims[,3]),], 
	    my_dims[!is.na(my_dims[,1]) | !is.na(my_dims[,3]),], 
	    my_mat[!is.na(my_dims[,1]) | !is.na(my_dims[,3]),])

	cat("The dimensions of the table for Zegami are ", dim(Zeg_table), "\n")
	if(length(sce@metadata$experiment_info$sample_id)==1) {
		roi <- tolower(sce@metadata$experiment_info$ROI)
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, sample, roi, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	} else if(length(unique(sce@metadata$experiment_info$sample_name)==1)) {
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, sample, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")

	} else {
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	}

#	write.table(Zeg_table, file=file.path(opt$outdir, "output_tables", paste0(opt$label, "_Zegamitable.txt")), quote=F, row.names=F, sep="\t")

}

###########################
# saving
rm(R_pheno_out_all, R_pheno_out_noDNA, ks, my_mat, other_vars, my_clusters2, Zeg_table)
session <- sessionInfo()
save.image(file=file.path(opt$outdir, "RData", paste0("Rpheno_", opt$label, ".RData")))

cat("Updating the sce object in ", file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")),"\n")
save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")))

