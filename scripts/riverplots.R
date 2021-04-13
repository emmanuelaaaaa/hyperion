#!/usr/bin/env Rscript
##################
## Description 
##
## This script creates the riverplots from one annotation to another.
##
###########################
# libraries 
suppressPackageStartupMessages({
  	require(optparse, quietly=T)
  	require(CATALYST, quietly=T)
  	require(SingleCellExperiment, quietly=T)
  	require(dplyr, quietly=T)
  	require(cowplot, quietly=T)
  	require(gridExtra, quietly=T)
  	require(ggplot2, quietly=T)
  	require(readxl, quietly=T)
  	require(networkD3, quietly=T)
  	require(htmltools, quietly=T)
})

# parsing the arguments
option_list = list(
    make_option(c("--infile_dir"), type="character", default=NULL, help="Directory with the input files", metavar="character"),
    make_option(c("--infile_combined"), type="character", default=NULL, help="File that contains the sce object for which to do the clustering", metavar="character"),
    make_option(c("--infile_individ"), type="character", default=NULL, help="File that contains the sce object for which we have the annotations", metavar="character"),
    make_option(c("--infile_annotations"), type="character", default=NULL, help="Excel with annotations per cluster", metavar="character"),
    make_option(c("--sheets_to_read"), type="character", default=NULL, help="Names of spreadsheets to import (need to have the same names as the sce objects that are being imported (case doesn't matter).", metavar="character"),
    make_option(c("--fromindtocombined"), type="logical", default=T, help="Include flag when the annnotations will be transfered from the large obj to the individuals (samples or ROIs)."),
    make_option(c("-o","--outdir"), type="character", default=NULL,  help="output directory", metavar="character"),
    make_option(c("--label"), type="character", default=NULL,  help="label for outputs")
); 
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$infile_dir),is.null(opt$infile_combined), is.null(opt$infile_individ), is.null(opt$infile_annotations), is.null(opt$sheets_to_read), is.null(opt$outdir),is.null(opt$label))) {
  print_help(opt_parser)
  stop("Arguments missing.n", call.=FALSE)
}

cat("Using the options:\n")
print(opt)

########### 
## loading the sce objects
cat("Loading all the individual files... \n")
infile_individ <- unlist(strsplit(opt$infile_individ, split=","))
ind_files <- c()
for (i in infile_individ) {
	load(file.path(opt$infile_dir,i))
	filename_ind <- gsub("sceobj_Jan21_", replacement="", i) %>% gsub("(_scaledtrim).RData", replacement="",.)
	assign(toupper(filename_ind),sce)
	ind_files <- c(ind_files,toupper(filename_ind))
}
cat("Loading the combined file... \n")
load(file.path(opt$infile_dir,opt$infile_combined))
filename <- gsub("sceobj_Jan21_", replacement="", opt$infile_combined) %>% gsub("(_scaledtrim).RData", replacement="",.) %>% toupper(.)
assign(filename,sce)

sheets_to_read <- unlist(strsplit(opt$sheets_to_read, split=","))
for (i in sheets_to_read) {
	annotated_clusters <- read_excel(opt$infile_annotations, sheet=i)
	sample <- toupper(i)
	sce_new <- get(sample)
	colData(sce_new)$annotation <- as.character(annotated_clusters$Annotation[colData(sce_new)$phenograph_cluster_scaledtrim_k30])
	assign(sample,sce_new)
}
rm(sce,sce_new)

## merging the relevant columns
cat("Merging the annotations... \n")
if(opt$fromindtocombined) {
	sce <- get(filename)
	samples_ann <- colData(sce)[,c("cellID_name","phenograph_cluster_scaledtrim_k30")]

	temp <- get(ind_files[1])
	temp <- colData(temp)[,c("cellID_name","annotation")]
	for (j in ind_files[2:length(ind_files)]) {
		temp <- rbind(temp, get(j)@colData[,c("cellID_name","annotation")])
	}

	samples_ann <- base::merge(samples_ann, temp, 
			by.x="cellID_name", by.y="cellID_name")
} else {
	sce <- get(filename)
	samples_ann <- colData(sce)[,c("cellID_name","annotation")]

	temp <- get(ind_files[1])
	temp <- colData(temp)[,c("cellID_name","phenograph_cluster_scaledtrim_k30")]
	for (j in ind_files[2:length(ind_files)]) {
		temp <- rbind(temp, get(j)@colData[,c("cellID_name","phenograph_cluster_scaledtrim_k30")])
	}

	samples_ann <- base::merge(samples_ann, temp, 
			by.x="cellID_name", by.y="cellID_name")

	samples_ann$ROI <- sapply(samples_ann$cellID_name, function(x) strsplit(x,"_")[[1]][5])
}

############
## making the plots
cat("Making the plots... \n")
if(opt$fromindtocombined) {
	edges_obj1 <- as.data.frame(table(samples_ann$annotation, samples_ann$phenograph_cluster_scaledtrim_k30))
	names(edges_obj1) <- c("N1","N2","Value")
	edges_obj1$N1 <- as.character(edges_obj1$N1)
	edges_obj1$N2 <- as.character(edges_obj1$N2)
	edges_obj1 <- edges_obj1 %>% subset(Value>0)

	myDF_obj1 <- list(
    	nodes=data.frame(name=unique(c(edges_obj1$N1, edges_obj1$N2))),
    		links=data.frame(source= match(edges_obj1$N2, unique(c(edges_obj1$N1, edges_obj1$N2)))-1,
                     target= match(edges_obj1$N1, unique(c(edges_obj1$N1, edges_obj1$N2)))-1,
                     value = edges_obj1$Value
    		)
	)

	sankeyNetwork(Links = myDF_obj1$links, Nodes = myDF_obj1$nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "cells", fontSize = 18, width=1500, height=1000,
              fontFamily = "sans-serif", iterations = 0) %>% 
		saveNetwork(file=file.path(opt$outdir, paste0(opt$label, "_riverplot.html")))
} else {

	for(k in unique(samples_ann$ROI)) {
		samples_temp <- filter(samples_ann, ROI==k)
		edges_obj1 <- as.data.frame(table(samples_temp$annotation, samples_temp$phenograph_cluster_scaledtrim_k30))
		names(edges_obj1) <- c("N1","N2","Value")
		edges_obj1$N1 <- as.character(edges_obj1$N1)
		edges_obj1$N2 <- as.character(edges_obj1$N2)
		edges_obj1 <- edges_obj1 %>% subset(Value>0)

		myDF_obj1 <- list(
    		nodes=data.frame(name=unique(c(edges_obj1$N1, edges_obj1$N2))),
    			links=data.frame(source= match(edges_obj1$N2, unique(c(edges_obj1$N1, edges_obj1$N2)))-1,
         	            target= match(edges_obj1$N1, unique(c(edges_obj1$N1, edges_obj1$N2)))-1,
          	           value = edges_obj1$Value
    			)
		)

		sankeyNetwork(Links = myDF_obj1$links, Nodes = myDF_obj1$nodes, Source = "source",
         	     Target = "target", Value = "value", NodeID = "name",
          	    units = "cells", fontSize = 18, width=1500, height=1000,
            	  fontFamily = "sans-serif", iterations = 0) %>% 
			saveNetwork(file=file.path(opt$outdir, paste0(opt$label, "_ROI",k,"_riverplot.html")))


	}

}
