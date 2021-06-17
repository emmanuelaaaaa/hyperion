#!/usr/bin/env Rscript
##################
## Description 
##
## This script performs the integration with harmony and then clustering with Rphenograph. It will load a file and run harmony on the various datatransformations specified with the relevant flag and then
## run Phenograph on the results. This the script I run after the exploratory harmony_and_Rpheno_clustering.R which is similar to this but has more options for the integration.
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
    make_option(c("--var"), type="character", default=NULL, help="The var from the data that will be removed. It needs to be a column in colData. (e.g. sample_name or sample_id)"),
    make_option(c("--annot"), type="character", default=NULL,  help="File with annotations to be loaded for the plots.Can include more than one option with comma separated files (will be added as annotation2, etc). (optional)", metavar="character"),
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
### setting up a better palette
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
### Create necessary folders

if (!dir.exists(file.path(opt$outdir, "figures"))) dir.create(file.path(opt$outdir, "figures"), showWarnings = FALSE)
if (!dir.exists(file.path(opt$outdir, "RData"))) dir.create(file.path(opt$outdir, "RData"), showWarnings = FALSE)

### loading the files and running the analysis

cat("Loading ",opt$infile, " \n")	
load(opt$infile)

datatransf <- unlist(strsplit(opt$datatransf, split=","))

if(!is.null(opt$annot)) {
	annotations <- unlist(strsplit(opt$annot, split=","))
	for (i in annotations) {
		annot <- read.table(i, header=T, stringsAsFactors=F, sep="\t")
		annot <- annot[match(colData(sce)$cellID_name,annot$cellID_name),]
		if (i==annotations[1]) {
			name <- "annotation"
		} else {
			name <- paste0("annotation", which(i==annotations))
		}
		colData(sce)[,name] <- annot$annotation
	}
}

if(is.null(opt$outdir_fig)) opt$outdir_fig <- opt$outdir

for (i in datatransf) {

	pca <- T
	cat("Running Harmony and Rphenograph for ",i, " PCA ",pca," \n")	

	exp <- assay(sce, i)[rowData(sce)$marker_name[rowData(sce)$marker_class=="type"],]
	my_harmony_embeddings <- HarmonyMatrix(
  		data_mat  = exp,
  		meta_data = colData(sce)[,opt$var],
		do_pca    = pca)

	R_pheno_out <- Rphenograph(my_harmony_embeddings, k=30)
	cat("\n The modularity of the graph is ", modularity(R_pheno_out[[2]]), "\n")
	sce[[paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30")]] <- factor(membership(R_pheno_out[[2]]))
	sce_umap <- uwot::umap(my_harmony_embeddings,  n_components = 2)
	reducedDims(sce)[[paste0("HarmIntegr_pcaTRUE_UMAP_",i)]] <- sce_umap

	sce_pheno <- data.frame(reducedDims(sce)[["UMAP_scaledtrim"]],reducedDims(sce)[[paste0("HarmIntegr_pcaTRUE_UMAP_",i)]])
	colnames(sce_pheno ) <- c("UMAP_scaledtrim_1","UMAP_scaledtrim_2","HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2")
	sce_pheno <- cbind(sce_pheno, sce@colData)

	p1a <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
	p2a <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=sample_id)) + geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))

	if(is.null(opt$annot)) {
		p1b <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=phenograph_cluster_scaledtrim_k30))+ scale_color_manual(values=c30) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=phenograph_cluster_scaledtrim_k30)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
	} else {
		p1b <- ggplot(as.data.frame(sce_pheno), aes(UMAP_scaledtrim_1, UMAP_scaledtrim_2, col=annotation))+ scale_color_manual(values=c30) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2b <- ggplot(as.data.frame(sce_pheno), aes(HarmIntegr_pcaTRUE_UMAP_1, HarmIntegr_pcaTRUE_UMAP_2, col=annotation)) + scale_color_manual(values=c30)+ geom_point(size = 0.5)+theme_bw()+ guides(colour = guide_legend(override.aes = list(size=4)))
	}

	if (length(unique(sce_pheno[,paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30")]))<=30) {
		p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
		p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+ scale_color_manual(values=c30)+ guides(colour = guide_legend(override.aes = list(size=4)))
	} else {
		p1c <- ggplot(as.data.frame(sce_pheno), aes_string("UMAP_scaledtrim_1", "UMAP_scaledtrim_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw() + guides(colour = guide_legend(override.aes = list(size=4)))
		p2c <- ggplot(as.data.frame(sce_pheno), aes_string("HarmIntegr_pcaTRUE_UMAP_1", "HarmIntegr_pcaTRUE_UMAP_2", col=paste0("harmony_pcaTRUE_phenograph_cluster_",i,"_k30"))) + geom_point(size = 0.5)+theme_bw()+guides(colour = guide_legend(override.aes = list(size=4)))
	}

	pall <- plot_grid(p1a,p1b,p1c, p2a,p2b,p2c, ncol=3)
	ggsave(pall, filename=file.path(opt$outdir_fig, "figures", opt$label, paste0("harmonyInt_Rpheno_UMAP_", opt$label, "_",i,".png")), width = 30, height = 15, dpi = 300, type="cairo-png")

	one_plot_heatmap <- plotExprHeatmap_updated(sce, cluster_id=paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30"),
		features=rowData(sce)$marker_name[rowData(sce)$marker_class=="type"], row_anno = TRUE, bars=T) 
	cat("Plotting the heatmap plots of the  phenograph clusters... \n")
	pdf(file=file.path(opt$outdir_fig, "figures", opt$label, paste0("Rpheno_harmonyInt_heatmaps_",i,".pdf")), height=10, width=15)
	print(one_plot_heatmap)
	dev.off()

	one_plot_exprs <- plotClusterExprs_updated(sce, cluster_id=paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30"), features = rowData(sce)$marker_name[rowData(sce)$marker_class=="type"]) 
	cat("Plotting the expression density plots of the phenograph clusters... \n")
	pdf(file=file.path(opt$outdir_fig, "figures", opt$label, paste0("Rpheno_harmonyInt_exprsdens_",i,".pdf")), height=10, width=15)
	print(one_plot_exprs)
	dev.off()

	one_plot_abund_box <- plotAbundances_updated(sce, cluster_id=paste0("harmony_pca", pca,"_phenograph_cluster_",i,"_k30"), group_by="sample_name", by="cluster_id" ) 
	one_plot_abund_stripe <- plotAbundances_updated(sce, cluster_id=paste0("phenograph_cluster_",i,"_k30"), group_by="sample_name", by="sample_id") 
	cat("Plotting the cluster abundances of the phenograph clusters per sample... \n")
	pdf(file=file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_harmonyInt_abund_",i,"_k30.pdf")), height=10, width=15)
	print(one_plot_abund_box)
	print(one_plot_abund_stripe)
	dev.off()

	if(opt$outdir_fig!=opt$outdir_fig) {
		file.symlink(file.path(opt$outdir_fig, "figures", opt$label, paste0("harmonyInt_Rpheno_UMAP_", opt$label, "_",i,".png")),file.path(opt$outdir, "figures", opt$label, paste0("harmonyInt_Rpheno_UMAP_", opt$label, "_",i,".png")))
		file.symlink(file.path(opt$outdir_fig, "figures", opt$label, paste0("Rpheno_harmonyInt_heatmaps_",i,".pdf")),file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_harmonyInt_heatmaps_",i,".pdf")))
		file.symlink(file.path(opt$outdir_fig, "figures", opt$label, paste0("Rpheno_harmonyInt_exprsdens_",i,".pdf")),file.path(opt$outdir, "figures", opt$label, paste0("Rpheno_harmonyInt_exprsdens_",i,".pdf")))
	}

}

save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_harmonycorr_", opt$label, "_final.RData")))

### Creating table for Zegami

if (opt$outputtable) {
	source("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/github_scripts/scripts/zegami_function.R")
	cat("\nMaking the table for Zegami...\n")
	makeZegami(sce, output_dir=opt$outdirZegami, label=opt$label, zegami_suffix="_integration_final")
}
