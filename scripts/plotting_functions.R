library(CATALYST)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

plotExprHeatmap_updated <- function(sce, features = NULL, cluster_id = "cluster_id", 
    k = "meta20", q = 0.01, 
    row_anno = TRUE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, 
    bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
    k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, distance = c("euclidean", 
        "maximum", "manhattan", "canberra", "binary", "minkowski"), 
    linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
        "median", "centroid", "ward.D2")) 
 {
  sce <- sce[unique(CATALYST:::.get_features(sce, features)), ]
  z <- assay(sce, "exprs")
  z <- CATALYST:::.scale_exprs(z, 1, 0.01)
  assay(sce, "exprs", FALSE) <- z
  z <- CATALYST:::.agg(sce, by=cluster_id, "median", "exprs")
  z <- t(z)

  qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
  lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  cat_cols <- CATALYST:::.cluster_cols
  sce$cluster_id <- factor(sce[[cluster_id]])

  myColors <- c("#A6CEE3",  "#B2DF8A", "#33A02C", "#FB9A99", "#B15928",  "#FF7F00", "#CAB2D6", "#6A3D9A","#E31A1C", "#FFD700","#808000","#1F78B4","#FDBF6F","#8B7D6B","#D9D9D9", "#000000")#, "#000000")
  left_anno <- CATALYST:::.anno_clusters(sce, 
            k="meta20", m=NULL, k_pal=myColors, m_pal = myColors)
   if (bars) {
        right_anno <- CATALYST:::.anno_counts(sce[["cluster_id"]], perc)
    }
    else right_anno <- NULL


 Heatmap(matrix = z,show_column_dend = F, name ="median scaled expression", col = colorRamp2(seq(min(z), 
        max(z), l = n <- 100), colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(n)), heatmap_legend_param = lgd_aes,
 		left_annotation = left_anno, right_annotation = right_anno, 
        rect_gp = gpar(col = "white"),)
}

plotClusterExprs_updated <- function (x, cluster_id = "cluster_id", features = "type") 
{
    x$cluster_id <- x[[cluster_id]]
    features <- CATALYST:::.get_features(x, features)
    ms <- t(CATALYST:::.agg(x[features, ], "cluster_id", "median"))
    d <- dist(ms, method = "euclidean")
    o <- hclust(d, method = "average")$order
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.table::data.table(data.frame(t(es), cd, check.names = FALSE))
    df <- data.table::melt(df, id.vars = names(cd), variable.name = "antigen", 
        value.name = "expression")
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    df$cluster_id <- factor(df$cluster_id, levels = rev(c("avg", 
        levels(x$cluster_id)[o])), labels = rev(c("average", 
        paste0(names(fq), " (", fq, "%)")[o])))
    ggplot(df, aes_string(x = "expression", y = "cluster_id", 
        col = "avg", fill = "avg")) + facet_wrap(~antigen, scales = "free_x", 
        nrow = 2) + ggridges::geom_density_ridges(alpha = 0.2) + ggridges::theme_ridges() + 
        theme(legend.position = "none", strip.background = element_blank(), 
            strip.text = element_text(face = "bold"))
}

plotAbundances_updated <- function (x, cluster_id, by = c("sample_id", "cluster_id"), 
    group_by = "condition", shape_by = NULL, col_clust = TRUE, 
    distance = c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski"), linkage = c("average", "ward.D", 
        "single", "complete", "mcquitty", "median", "centroid", 
        "ward.D2"), k_pal = CATALYST:::.cluster_cols) 
{
    x$cluster_id <- x[[cluster_id]]
    by <- match.arg(by)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)
    stopifnot(is.logical(col_clust), length(col_clust) == 1)
    shapes <- CATALYST:::.get_shapes(x, shape_by)
    if (is.null(shapes)) 
        shape_by <- NULL
    if (by == "sample_id") {
        nk <- nlevels(x[[cluster_id]])
        if (length(k_pal) < nk) 
            k_pal <- colorRampPalette(k_pal)(nk)
    }
    ns <- table(cluster_id = x[[cluster_id]], sample_id = sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
    if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 
        1) {
        d <- dist(t(fq), distance)
        h <- hclust(d, linkage)
        o <- colnames(fq)[h$order]
        df$sample_id <- factor(df$sample_id, o)
    }
    p <- ggplot(df, aes_string(y = "Freq")) + labs(x = NULL, 
        y = "Proportion [%]") + theme_bw() + theme(panel.grid = element_blank(), 
        strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, 
            color = NA), axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.key.height = unit(0.8, "lines"))
    switch(by, sample_id = p + (if (!is.null(group_by)) facet_wrap(group_by, 
        scales = "free_x")) + geom_bar(aes_string(x = "sample_id", 
        fill = "cluster_id"), position = "fill", stat = "identity") + 
        scale_fill_manual("cluster_id", values = k_pal) + scale_x_discrete(expand = c(0, 
        0)) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 
        100, 25)) + theme(panel.border = element_blank(), panel.spacing.x = unit(1, 
        "lines")), cluster_id = {
        p <- p + scale_shape_manual(values = shapes) + guides(col = guide_legend(order = 1, 
            override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 3)))
        if (is.null(group_by)) {
            p + geom_boxplot(aes_string(x = "cluster_id"), alpha = 0.2, 
                position = position_dodge(), outlier.color = NA) + 
                geom_point(aes_string("cluster_id", shape = shape_by), 
                  position = position_jitter(width = 0.2))
        } else {
            p + facet_wrap("cluster_id", scales = "free_y", ncol = 4) + 
                geom_boxplot(aes_string(x = group_by, color = group_by, 
                  fill = group_by), position = position_dodge(), 
                  alpha = 0.2, outlier.color = NA, show.legend = FALSE) + 
                geom_point(aes_string(x = group_by, col = group_by, 
                  shape = shape_by), position = position_jitter(width = 0.2))
        }
    })
}

