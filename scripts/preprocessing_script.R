#!/usr/bin/env Rscript
##################
## Description 
##
## This script runs the initial processing and some exploratory plots.
##
###########################
## libraries 

suppressPackageStartupMessages({
  require(optparse, quietly=T)
  require(CATALYST, quietly=T)
  require(flowCore, quietly=T)
  require(dplyr, quietly=T)
  require(SingleCellExperiment, quietly=T)
  require(ggplot2, quietly=T)
  require(cowplot, quietly=T)
  require(gridExtra, quietly=T)
  require(purrr, quietly=T)
  require(plyr, quietly=T)
  require(tidyr, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--input_dir"), type="character", default=NULL, help="Directory with the input files to process", metavar="character"),
    make_option(c("--panel_file"), type="character", default=NULL, help="Panel text file with the markers to process (needs to have channel_name,marker_name,marker_class).", metavar="character"),
    make_option(c("--metadata_file"), type="character", default=NULL, help="Text file with metadata information (sample_id, sample_name, condition, ROI)", metavar="character"),
    make_option(c("--transform"), type="character", default=FALSE, help="Arcsinh transform the data? (default=F)", metavar="character"),
    make_option(c("--only_include"), type="character", default=NULL, help="Name of specific samples (sample_id or sample_name) to process (comma separated), if NULL it will process everything in the metadata file.", metavar="character"),
    make_option(c("--samples_to_exclude"), type="character", default=NULL, help="Sample_id to exclude from the analysis(comma separated).", metavar="character"),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$input_dir),is.null(opt$panel_file),is.null(opt$metadata_file),is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

### Create necessary folders

if (!dir.exists(file.path(opt$outdir, "figures/"))) dir.create(file.path(opt$outdir, "figures/"), showWarnings = FALSE)
if (!dir.exists(file.path(opt$outdir, "figures/exploratory_plots"))) dir.create(file.path(opt$outdir, "figures/exploratory_plots"), showWarnings = FALSE)
if (!dir.exists(file.path(opt$outdir, "RData"))) dir.create(file.path(opt$outdir, "RData"), showWarnings = FALSE)
if (!dir.exists(file.path(opt$outdir, "figures", opt$label))) dir.create(file.path(opt$outdir, "figures", opt$label), showWarnings = FALSE)

### import and process the data

cat("Starting the preprocessing script... \n")
panel <- read.table(opt$panel_file, header=T, stringsAsFactors=F)
md_file <- read.table(opt$metadata_file, header=T, stringsAsFactors=F, sep="\t")
md_file$n_cells <- NA
if (!is.null(opt$samples_to_exclude)) {
	cat("Excluding samples: ",opt$samples_to_exclude,"\n")
	samples_to_exclude <- unlist(strsplit(opt$samples_to_exclude, split=","))
	md_file <- filter(md_file, !(sample_id %in% samples_to_exclude))
}
if (!is.null(opt$only_include)) {
	cat("Only running samples: ",opt$only_include,"\n")
	samples_to_include <- unlist(strsplit(opt$only_include, split=","))
	if(grepl("ROI", samples_to_include[1])) {
		md_file <- filter(md_file, (sample_id %in% samples_to_include))
	} else {
		md_file <- filter(md_file, (sample_name %in% samples_to_include))
	}
}

antib_data <- data.frame()
celldata <- data.frame()

for (i in md_file$sample_id) { 

	roi <- tolower(md_file[md_file$sample_id==i,]$ROI)
	sample <- paste0("sample_", strsplit(i, split="_")[[1]][3])
	cond <- tolower(strsplit(i, split="_")[[1]][1])
	filename <- file.path(opt$input_dir, cond, sample, roi, "Zegami/cellData.csv")

	exp1 <- read.csv(filename, header=T, stringsAsFactors=F)
	exp1$X <- NULL
	exp1$cellID_name <- gsub(".png", replacement="", exp1$Image.Name)
	names(exp1) <- gsub("_mean", replacement="",names(exp1))

	if (!file.exists(file.path(opt$outdir, "figures/exploratory_plots", paste0( i, "_varshistograms.pdf")))) {
		hist_plot <- exp1 %>%
  					keep(is.numeric) %>% 
  					gather() %>% 
  					ggplot(aes(value)) +
    					facet_wrap(~ key, scales = "free") +
    					geom_histogram()
		ggsave(hist_plot, filename=file.path(opt$outdir, "figures/exploratory_plots", paste0(i, "_varshistograms.pdf")), width = 12, height = 8, dpi = 300)
	}

	antib_data_temp <- select(exp1, panel$channel_name)
	celldata_temp <- select(exp1, setdiff(names(exp1),panel$channel_name))

	antib_data <- rbind(antib_data, antib_data_temp)
	celldata <- rbind.fill(celldata, celldata_temp)

	md_file[md_file$sample_id==i,"n_cells"] <-  nrow(exp1)
}

row.names(antib_data) <- celldata$cellID_name
row.names(celldata) <- celldata$cellID_name
row.names(panel) <- panel$channel_name
celldata$sample_id <- sapply(celldata$Image.Name, function(x) strsplit(x, split="__")[[1]][1])

if (opt$transform) {
	sce <- SingleCellExperiment(list(counts=t(as.matrix(antib_data))),
		colData=celldata,
		rowData=panel,
		metadata=list(experiment_info=md_file)
		)
	sce <- CATALYST:::.transform(sce, 5) # asinh(value/5)
	assay(sce, "scaled", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, 0)
	assay(sce, "scaledtrim", FALSE) <- CATALYST:::.scale_exprs(assay(sce, "exprs"), 1, 0.01)
} else {
	sce <- SingleCellExperiment(list(exprs=t(as.matrix(antib_data))),
    		colData=celldata,
    		rowData=panel,
    		metadata=list(experiment_info=md_file)
	)
}

### add the sample as column, removing the ROI 
colData(sce)$sample_name <- sapply(colData(sce)$sample_id, function(x) substring(x, first=1, last=nchar(x)-6))
colData(sce)$condition <- sapply(colData(sce)$sample_name, function(x) substring(x, first=1, last=nchar(x)-9))

### make expression plots, number of cells, heatmap and PCAs
expr_smoothhisto <- plotExprs(sce, color_by = "sample_id") 
expr_heat <- plotExprHeatmap(sce)
ncounts <- ggplot(metadata(sce)$experiment_info, aes(reorder(sample_id, -n_cells), n_cells, fill=condition)) + geom_bar(stat = "identity", width = 0.75) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size=12)) 
pdf(file=file.path(opt$outdir, "figures/exploratory_plots/", paste0(opt$label, "_transf", opt$transform, "_QCplots.pdf")), height=15, width=20)
    	ncounts
    	expr_smoothhisto
    	expr_heat
dev.off()

# running the scaled and non-scaled PCA:
if (length(unique(sce@metadata$experiment_info$sample_name))>1) {
	cs_by_s <- split(seq_len(ncol(sce)), sce$sample_id)
	es <- as.matrix(assay(sce, "exprs")[rowData(sce)$marker_name[rowData(sce)$marker_class=="type"],])
	ms <- vapply(cs_by_s, function(cs) rowMedians(es[, cs, drop = FALSE]), 
 	       numeric(length(which(rowData(sce)$marker_class=="type"))))
	rownames(ms) <- rowData(sce)$marker_name[rowData(sce)$marker_class=="type"]
	md1 <- metadata(sce)$experiment_info

	for (scaling in c(T,F)) {
		if (any(apply(t(ms), 2, min)==apply(t(ms), 2, max))) {
			df <- t(ms)[,apply(t(ms), 2, min)!=apply(t(ms), 2, max)]
		} else {
			df <- t(ms)
		}
		pca <- prcomp(df, scale=scaling)
		df_to_plot <- data.frame(pca$x[, 1:4])
		m <- match(rownames(df_to_plot), md1$sample_id)
		df_to_plot <- data.frame(df_to_plot, md1[m, ])
		names(df_to_plot)[1:4] <- paste0("PCA",1:4)
	    	p1a <- ggplot(data.frame(df_to_plot), aes(x=PCA1, y=PCA2)) +
	        	geom_point(aes(color=sample_name, shape=condition), size = 4) 
	    	p1b <- ggplot(data.frame(df_to_plot), aes(x=PCA1, y=PCA2)) +
		        geom_point(aes(color=sample_name, shape=condition), size = 2) 
	    	p2 <- ggplot(data.frame(df_to_plot), aes(x=PCA1, y=PCA3)) + theme(legend.position = "none") +
	        	geom_point(aes(color=sample_name, shape=condition), size = 2) 
	    	p3 <- ggplot(data.frame(df_to_plot), aes(x=PCA1, y=PCA4)) + theme(legend.position = "none") +
	        	geom_point(aes(color=sample_name, shape=condition), size = 2) 
	    	p4 <- ggplot(data.frame(df_to_plot), aes(x=PCA2, y=PCA3)) + theme(legend.position = "none") +
	        	geom_point(aes(color=sample_name, shape=condition), size = 2) 
	    	p5 <- ggplot(data.frame(df_to_plot), aes(x=PCA2, y=PCA4)) + theme(legend.position = "none") +
	        	geom_point(aes(color=sample_name, shape=condition), size = 2) 
	    	p6 <- ggplot(data.frame(df_to_plot), aes(x=PCA3, y=PCA4)) + theme(legend.position = "none") +
	        	geom_point(aes(color=sample_name, shape=condition), size = 2) 
		p7 <- plot_grid(p1b  + theme(legend.position = "none"), p2, p3, NULL, p4,p5, NULL, NULL, p6, ncol=3)
		pdf(file=file.path(opt$outdir, "figures/exploratory_plots/", paste0(opt$label,"_transf", opt$transform, "_PCA_1to4_scale", scaling,".pdf")), height=12, width=12)
			invisible(lapply(list(p1a, p7), print))
		dev.off()
		pdf(file=file.path(opt$outdir, "figures/exploratory_plots/", paste0(opt$label,"_transf", opt$transform, "_PCAbiplot_scale", scaling,".pdf")), height=12, width=12)
			biplot(pca)
		dev.off()

	}
}

### running dimentionality reduction
##changing the names of the location x,y because it's interfering with the UMAP/TSNE
names(colData(sce))[names(colData(sce))=="x"] <- "LocX"
names(colData(sce))[names(colData(sce))=="y"] <- "LocY"

if (opt$transform) {
	for (i in c("exprs","scaled","scaledtrim")) {
		cat("Running dimentionality reduction for ", i, " ... \n")
		cat("starting with TSNE... \n")
		sce <- runDR(sce, dr = "TSNE",  features = "type", assay=i)
		reducedDim(sce, paste0("TSNE_",i)) <- reducedDim(sce, "TSNE")

		cat("UMAP... \n")
		sce <- runDR(sce, dr = "UMAP", features = "type", assay=i)
		reducedDim(sce, paste0("UMAP_",i)) <- reducedDim(sce, "UMAP")

	}
} else {
	cat("Running dimentionality reduction for exprs ... \n")
	cat("starting with TSNE... \n")
	sce <- runDR(sce, dr = "TSNE",  features = "type", assay="exprs")

	cat("UMAP... \n")
	sce <- runDR(sce, dr = "UMAP", features = "type", assay="exprs")

}

save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")))

if(length(sce@metadata$experiment_info$sample_id)>1) {
	cat("Plotting the UMAP and TSNE plots by sample_name and sample_id... \n")
	plot_tsne_c <-plotDR(sce, "TSNE_exprs", color_by = "sample_id")
	plot_umap_c <-plotDR(sce, "UMAP_exprs", color_by = "sample_id")

	if(length(unique(sce@metadata$experiment_info$sample_name))>1) {
		plot_tsne_b <-plotDR(sce, "TSNE_exprs", color_by = "sample_name")
		plot_umap_b <-plotDR(sce, "UMAP_exprs", color_by = "sample_name")
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

session <- sessionInfo()

cat("Creating the sce object in ", file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")),"\n")
save(sce, file=file.path(opt$outdir, "RData", paste0("sceobj_", opt$label, ".RData")))

save(panel, md_file, session, file=file.path(opt$outdir, "RData", paste0("sce_initial_",opt$label, ".RData")))
