#!/usr/bin/env Rscript
##################
## Description 
##
## This script runs the clustering with FlowSOM and makes the relevant plots.
##
###########################
## libraries 

suppressPackageStartupMessages({
  require(optparse, quietly=T)
  require(CATALYST, quietly=T)
  require(flowCore, quietly=T)
  require(dplyr, quietly=T)
  require(SingleCellExperiment, quietly=T)
  require(ConsensusClusterPlus, quietly=T)
  require(FlowSOM, quietly=T)
  require(ggplot2, quietly=T)
  require(cowplot, quietly=T)
  require(gridExtra, quietly=T)
  require(purrr, quietly=T)
  require(tidyr, quietly=T)
})

source("/t1-data/user/erepapi/Fellowship/COVID19_CyTOF/scripts/other_plots.R")

# parsing the arguments
option_list = list(
    make_option(c("--infile"), type="character", default=NULL, help="File that contains the sce object for which to do the clustering", metavar="character"),
    make_option(c("--datatransf"), type="character", default="exprs",  help="Data transformation for running the clustering on. Options: counts=no transformation, exprs= arcsinh 
		transformation, scaled=scaled transformation of the exprs, scaledtrim=the scaled and trimmed values (q=0.01) of the exprs. Can include more than one option 
		in comma separated style, e.g. exprs,scaled. (default=exprs)"),
    make_option(c("--xdim"), type="integer", default=10,  help="xdim for FlowSOM"),
    make_option(c("--ydim"), type="integer", default=10,  help="ydim for FlowSOM"),
    make_option(c("--maxmeta"), type="integer", default=NULL,  help="maximum metaclusters for FlowSOM"),
    make_option(c("--outputtable"), type="logical", action="store_true", help="Include flag to create the Zegami table."),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile),is.null(opt$xdim),is.null(opt$ydim),is.null(opt$outdir),is.null(opt$maxmeta),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

dir.create(file.path(opt$outdir, "figures", opt$label), showWarnings = FALSE)

########### 
## loading the sce object from the initial QC analysis - no need to preprocess again
cat("Starting the clustering with FlowSOM... \n")
load(opt$infile)

##changing the names of the location x,y because it's interfering with the UMAP/TSNE
names(colData(sce))[names(colData(sce))=="x"] <- "LocX"
names(colData(sce))[names(colData(sce))=="y"] <- "LocY"
colData(sce)$sample_name <- sapply(colData(sce)$sample_id, function(x) strsplit(x, split="_")[[1]][1])

datatransf <- unlist(strsplit(opt$datatransf, split=","))
###########################
# running dimentionality reduction

for (i in datatransf) {
	cat("Running dimentionality reduction for ", i, " ... \n")
	cat("starting with TSNE... \n")
	sce <- runDR(sce, dr = "TSNE",  features = "type", assay=i)
	reducedDim(sce, paste0("TSNE_",i)) <- reducedDim(sce, "TSNE")

	cat("UMAP... \n")
	sce <- runDR(sce, dr = "UMAP", features = "type", assay=i)
	reducedDim(sce, paste0("UMAP_",i)) <- reducedDim(sce, "UMAP")

	if(length(sce@metadata$experiment_info$sample_id)>1) {
		cat("Plotting the UMAP and TSNE plots by sample_name and sample_id... \n")
		plot_tsne_c <-plotDR(sce, "TSNE", color_by = "sample_id")
		plot_umap_c <-plotDR(sce, "UMAP", color_by = "sample_id")

		if(length(unique(sce@metadata$experiment_info$sample_name))>1) {
			plot_tsne_b <-plotDR(sce, "TSNE", color_by = "sample_name")
			plot_umap_b <-plotDR(sce, "UMAP", color_by = "sample_name")
			pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("TSNEUMAPpersample_", opt$label, ".pdf")), height=10, width=10)
			    print(plot_tsne_b)
			    print(plot_tsne_c)
			    print(plot_umap_b)
			    print(plot_umap_c)
			dev.off()
		} else {
		pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("TSNEUMAPpersample_", opt$label, ".pdf")), height=10, width=10)
		    print(plot_tsne_c)
		    print(plot_umap_c)
		dev.off()
		}
	}
}

###########################
## FlowSOM clustering (with code from CATALYST, but not calling CATALYST:::cluster so that the MST can be build as well)

meta_values <- paste0("meta",seq(from=(opt$maxmeta-10), to=opt$maxmeta, by=1))

for (i in datatransf) {
	fSOM <- ReadInput(flowFrame(t(assay(sce, i))))
	set.seed(123)
	fSOM <- BuildSOM(fSOM, colsToUse=rowData(sce)$marker_name[rowData(sce)$marker_class=="type"], xdim=opt$xdim, ydim=opt$ydim)
	fSOM <- BuildMST(fSOM, tSNE=TRUE)

	# FlowSOM consensus clustering
	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("consensus_plots_",i,".pdf")), height=15, width=15)
		mc <- suppressMessages(ConsensusClusterPlus(t(fSOM$map$codes),maxK = opt$maxmeta, reps = 100,
	        	distance = "euclidean", seed = 123, plot = NULL, clusterAlg = "hc", pItem = 0.9, pFeature = 1))
	dev.off()

	# store the clustering in the SCE object
	k <- opt$xdim * opt$ydim
	mcs <- seq_len(opt$maxmeta)[-1]
	codes <- data.frame(seq_len(k), map(mc[-1], "consensusClass"))
	codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
	colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", mcs))
	sce$cluster_id <- factor(fSOM$map$mapping[, 1])
	metadata(sce)$cluster_codes <- codes
	metadata(sce)$SOM_codes <- fSOM$map$codes
	metadata(sce)$delta_area <- CATALYST:::.plot_delta_area(mc)

	## plotting the UMAP/TSNE

	cat("Plotting the UMAP and TSNE plots with the meta-clusters... \n")
	one_plot_tsne <- function(meta)  plotDR(sce, paste0("TSNE_",i), color_by = meta)
	one_plot_umap <- function(meta)  plotDR(sce, paste0("UMAP_",i), color_by = meta)
	meta_plots_t <- plyr::llply(meta_values, one_plot_tsne)
	meta_plots_u <- plyr::llply(meta_values, one_plot_umap) 

	png(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_TSNE_", i, ".png")), height=1000, width=1500)
	    grid.arrange(grobs=meta_plots_t)
	dev.off()

	png(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_UMAP_", i, ".png")), height=1000, width=1500)
	    grid.arrange(grobs=meta_plots_u)
	dev.off()

	my_clusters2 <- as.data.frame(sapply(sce@metadata$cluster_codes, function(x) x[cluster_ids(sce)]))
	my_clusters2$cellID <- colData(sce)$cellID
	write.table(my_clusters2, file=file.path(opt$outdir, "output_tables/otherclusterings", paste0(opt$label, "_", i, "_clusters.txt")), quote=F, row.names=F, sep="\t")	

	## plotting the flowsom plotstars
	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_plotstars", "_", i, ".pdf")), height=8, width=12)
	for (j in (opt$maxmeta-10):opt$maxmeta ) {
	    PlotStars(fSOM, backgroundValues = as.factor(metadata(sce)$cluster_codes[,j]))
	}
	dev.off()

	## plotting the heatmaps, density and abundance plots per sample

	cat("Plotting the heatmap, density and abundance plots with the meta-clusters... \n")
	one_plot_heatmap <- function(meta)  plotExprHeatmap(sce, by="cluster_id", k = meta, m = NULL, features="type",
		         row_anno = TRUE, bars=T) #, draw_freqs = TRUE, scale=T)
	meta_plots_heat <- plyr::llply(meta_values, one_plot_heatmap)
	cat("Plotting the heatmap plots of the meta-clusters... \n")
	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_heatmaps", "_", i, ".pdf")), height=10, width=15)
	print(meta_plots_heat)
	dev.off()

	one_plot_exprs <- function(meta) plotClusterExprs(sce, k = meta, features = "type") 
	meta_plots_exprs <- plyr::llply(meta_values, one_plot_exprs)
	cat("Plotting the expression density plots of the meta-clusters... \n")
	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_exprsdens", "_", i,".pdf")), height=10, width=15)
	print(meta_plots_exprs)
	dev.off()

	if(length(sce@metadata$experiment_info$sample_id)>1) {
		one_plot_abundance <- function(meta) plotAbundances(sce, k = meta, by = "cluster_id", group_by = "sample_id") 
		meta_plots_abund <- plyr::llply(meta_values, one_plot_abundance)
		cat("Plotting the abundance plots of the meta-clusters... \n")
		pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_abundance_per_sampleid", "_", i, ".pdf")), height=10, width=15)
		print(meta_plots_abund)
		dev.off()
		if(length(unique(sce@metadata$experiment_info$sample_name))>1) {
			one_plot_abundance_id <- function(meta) plotAbundances(sce, k = meta, by = "cluster_id", group_by = "sample_name") 
			meta_plots_abund_id <- plyr::llply(meta_values, one_plot_abundance_id)
			cat("Plotting the abundance plots of the meta-clusters per sample_id... \n")
			pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("clustering_plots_abundance_per_sample", "_", i, ".pdf")), height=10, width=15)
			print(meta_plots_abund_id)
			dev.off()
		}
	}


	cat("Saving the sce object in ", file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, "_", i, ".RData")),"\n")
	save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, "_", i, ".RData")))

}

###########################
# saving the rest
session <- sessionInfo()
cat("Saving in", file.path(opt$outdir, "RData", paste0("clustering_analysis_", opt$label, ".RData")),"\n")
save(opt, session, fSOM, meta_plots_heat, 
     meta_plots_t, meta_plots_u,
    file=file.path(opt$outdir, "RData", paste0("clustering_analysis_", opt$label, ".RData")))
## the size of RData is getting too big and will make it very long to load, better to only keep the new stuff (the sce is stored on its own)

