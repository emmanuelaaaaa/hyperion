#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the integration with harmony and then clustering with Rphenograph. It will load a file and run harmony on the various datatransformations specified with the relevant flag and then
## run Phenograph on the results. At the moment it 
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
#  	require(CATALYST, quietly=T) # this is causing harmony to quit with weird exit message, but sourcing the plotting functions script that loads CATALYST too seems to be fine (!?!?!?)
  	require(SingleCellExperiment, quietly=T)
  	require(Rphenograph, quietly=T)
  	require(harmony, quietly=T)
  	require(cowplot, quietly=T)
  	require(gridExtra, quietly=T)
  	require(ggplot2, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the integration and clustering", metavar="character"),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--outdir_fig"), type="character", default=NULL,  help="output directory for figures (in datashare)", metavar="character"),
    make_option(c("--datatransf"), type="character", default="exprs",  help="Data transformation for running the clustering on. Options: counts=no transformation, exprs= arcsinh 
		transformation, scaled=scaled transformation of the exprs, scaledtrim=the scaled and trimmed values (q=0.01) of the exprs. Can include more than one option 
		in comma separated style, e.g. exprs,scaled. (default=exprs)"),
    make_option(c("--var"), type="character", default=NULL, help="The var from the data that will be removed. It needs to be a column in colData."),
    make_option(c("--annot"), type="character", default=NULL,  help="File with annotations to be loaded for the plots (optional).", metavar="character"),
    make_option(c("--outputtable"), type="logical", action="store_true", help="Include flag to create the Zegami table."),
    make_option(c("--outdirZegami"), type="character", default=NULL,  help="output directory for the Zegami table (to be same as input from preprocessing)", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$var),is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

source("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/github_scripts/scripts/plotting_functions.R")

###########################
c30 <- c("dodgerblue2", "#E31A1C", # red
                "green4",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black", "gold1",
                "skyblue2", "#FB9A99", # lt pink
                "palegreen2",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "darkgrey", "khaki2",
                "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
                "darkturquoise", "green1", "yellow4", "yellow3",
                "darkorange4", "brown","grey",
		"cyan","firebrick","hotpink","aquamarine"
              )

###########################
# loading the files and running the analysis

cat("Loading ",opt$infile, " \n")	
load(opt$infile)

if (opt$var=="sample") {
	colData(sce)$sample <- sapply(colData(sce)$sample_id, function(x) substring(x, first=1, last=nchar(x)-6))
}

datatransf <- unlist(strsplit(opt$datatransf, split=","))

if(!is.null(opt$annot)) {
	annot <- read.table(opt$annot, header=T, stringsAsFactors=F, sep="\t")
	annot <- annot[match(colData(sce)$cellID_name,annot$cellID_name),]
	colData(sce)$annotation <- annot$annotation
}

for (i in datatransf) {

	for (pca in c(T,F)) {
		cat("Running Harmony and Rphenograph for ",i, " PCA ",pca," \n")	

		exp <- assay(sce, i)[rowData(sce)$marker_name[rowData(sce)$marker_class=="type"],]
		my_harmony_embeddings <- HarmonyMatrix(
	  		data_mat  = exp,
	  		meta_data = colData(sce)[,opt$var],
			do_pca    = pca)
		
		if (pca==T) {
			R_pheno_out <- Rphenograph(my_harmony_embeddings, k=30)
			cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
			sce[[paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30")]] <- factor(membership(R_pheno_out[[2]]))
			sce_umap <- uwot::umap(my_harmony_embeddings,  n_components = 2)
			} else {
			R_pheno_out <- Rphenograph(t(my_harmony_embeddings), k=30)
			cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
			sce[[paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30")]] <- factor(membership(R_pheno_out[[2]]))
			sce_umap <- uwot::umap(t(my_harmony_embeddings),  n_components = 2)
		}
		reducedDims(sce)[[paste0("HarmIntegr_pca",pca,"_UMAP_",i)]] <- sce_umap

	}

	sce_pheno <- data.frame(reducedDims(sce)[["UMAP_scaledtrim"]],reducedDims(sce)[[paste0("HarmIntegr_pcaTRUE_UMAP_",i)]], reducedDims(sce)[[paste0("HarmIntegr_pcaFALSE_UMAP_",i)]])
	colnames(sce_pheno ) <- c("UMAP_scaledtrim_1","UMAP_scaledtrim_2","HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2","HarmIntegr_pcaFALSE_UMAP_1","HarmIntegr_pcaFALSE_UMAP_2")
	sce_pheno <- cbind(sce_pheno, sce@colData)

	p1a <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
	p2a <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
	p3a <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaFALSE_UMAP_1, HarmIntegr_pcaFALSE_UMAP_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))

	if(is.null(opt$annot)) {
		p1b <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=phenograph_cluster_scaledtrim_k30))+ scale_color_manual(values=c30) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=phenograph_cluster_scaledtrim_k30)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
		p3b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaFALSE_UMAP_1, HarmIntegr_pcaFALSE_UMAP_2, col=phenograph_cluster_scaledtrim_k30)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
	} else {
		p1b <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=annotation))+ scale_color_manual(values=c30) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=annotation)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
		p3b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaFALSE_UMAP_1, HarmIntegr_pcaFALSE_UMAP_2, col=annotation)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
	}

	if (length(unique(sce_pheno[,paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30")]))<=30) {
		p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
		p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
		p3c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaFALSE_UMAP_1", "HarmIntegr_pcaFALSE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
	} else {
		p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
		p3c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaFALSE_UMAP_1", "HarmIntegr_pcaFALSE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
	}

	if (length(unique(sce_pheno[,paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30")]))<=30) {
		p1d <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
		p2d <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
		p3d <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaFALSE_UMAP_1", "HarmIntegr_pcaFALSE_UMAP_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
	} else {
		p1d <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2d <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
		p3d <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaFALSE_UMAP_1", "HarmIntegr_pcaFALSE_UMAP_2", col=paste0("harmony_pcaFALSE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
	}

	pall <- plot_grid(p1a,p1b,p1c,p1d, p2a,p2b,p2c,p2d, p3a,p3b,p3c,p3d, ncol=4)
	ggsave(pall, filename=file.path(opt$outdir_fig, "figures", opt$label, paste0("harmonyInt_Rpheno_UMAP_", opt$label, "_",i,".png")), width = 30, height = 15, dpi = 300, type="cairo-png")

	}
}

save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_harmonycorr_", opt$label, ".RData")))

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

	my_mat <- t(assay(sce, "exprs"))

	my_dims <- data.frame(reducedDims(sce)[[1]])
	colnames(my_dims) <- paste(names(reducedDims(sce))[1], c(1,2), sep="_")
	for (i in 2:length(names(reducedDims(sce))) {
		my_dims_temp <- data.frame(reducedDims(sce)[[i]])
		colnames(my_dims_temp) <- paste(names(reducedDims(sce))[i], c(1,2), sep="_")
		my_dims <- cbind(my_dims, my_dims_temp)
	}

	other_vars <- sce@colData
	other_vars[,grepl("phenograph", names(other_vars))] <- paste_matrix("cl", as.matrix(other_vars[,grepl("phenograph", names(other_vars))]))

	Zeg_table <- cbind(other_vars, my_dims, my_mat)

	cat("The dimensions of the table for Zegami are ", dim(Zeg_table), "\n")
	if(length(sce@metadata$experiment_info$sample_id)==1) {
		roi <- tolower(sce@metadata$experiment_info$ROI)
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, sample, roi, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	} else if(length(unique(sce@metadata$experiment_info$sample_name))==1) {
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, sample, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")

	} else if(length(unique(sce@metadata$experiment_info$condition))==1) {
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(opt$outdirZegami, cond, "cellDataWithClustering.csv")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	} else {
		cat("The samples coming from different conditions, there is no relevant folder for saving the zegami. \n")
		filename <- file.path("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/output_tables", paste0(opt$label, "_Zegami.txt"))
		cat("Writing the table in ", filename," instead. \n")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	}

}
