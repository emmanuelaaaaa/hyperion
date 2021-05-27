makeZegami <- function(sce, output_dir, label, zegami_suffix="") {
	my_mat <- t(assay(sce, "exprs"))

	my_dims <- data.frame(reducedDims(sce)[[1]])
	colnames(my_dims) <- paste(names(reducedDims(sce))[1], c(1,2), sep="_")
	for (i in 2:length(names(reducedDims(sce)))) {
		my_dims_temp <- data.frame(reducedDims(sce)[[i]])
		colnames(my_dims_temp) <- paste(names(reducedDims(sce))[i], c(1,2), sep="_")
		my_dims <- cbind(my_dims, my_dims_temp)
	}

	other_vars <- sce@colData
	pheno_cols <- other_vars[,grepl("phenograph", names(other_vars))] 
	pheno_cols <- sapply(pheno_cols, function(x) sprintf("cl%02d", x))
	other_vars[,grepl("phenograph", names(other_vars))] <- pheno_cols

	Zeg_table <- cbind(other_vars, my_dims, my_mat)

	cat("The dimensions of the table for Zegami are ", dim(Zeg_table), "\n")
	if(length(sce@metadata$experiment_info$sample_id)==1) {
		roi <- tolower(sce@metadata$experiment_info$ROI)
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(output_dir, cond, sample, roi,  paste0("cellDataWithClustering", zegami_suffix, ".tab"))
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	} else if(length(unique(sce@metadata$experiment_info$sample_name))==1) {
		sample <- paste0("sample_", strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][3])
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(output_dir, cond, sample,  paste0("cellDataWithClustering", zegami_suffix, ".tab"))
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")

	} else if(length(unique(sce@metadata$experiment_info$condition))==1) {
		cond <- tolower(strsplit(sce@metadata$experiment_info$sample_name, split="_")[[1]][1])
		filename <- file.path(output_dir, cond,  paste0("cellDataWithClustering", zegami_suffix, ".tab"))
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	} else {
		cat("The samples coming from different conditions, there is no relevant folder for saving the zegami. \n")
		filename <- file.path("/t1-data/user/erepapi/Fellowship/Hyperion/COVID19/output_tables", paste0(label, "_Zegami", zegami_suffix, ".tab"))
		cat("Writing the table in ", filename," instead. \n")
		write.table(Zeg_table, file=filename, quote=F, row.names=F, sep="\t")
	}
}