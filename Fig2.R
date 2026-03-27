library(tidyverse)
library(reshape2)
library(data.table)
library(rhdf5)
library(Seurat)
library(scigenex)
library(patchwork)
library(scCustomize)
library(limma)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(circlize)
library(broom)
library(dplyr)
library(tidyr)
library(GSVA)
library(scalop)
library(ggsankey)
library(ggthemes)
library(msigdbr)
library(scCustomize)
library(SPATA2)
library(spatstat)
library(ComplexHeatmap)
library(viridis)
library(ggrepel)
library(clusterProfiler)
library(smplot2)
library(concaveman)  # For smoother outlines
library(dbscan)
library(sf)
library(ggalluvial)
library(monocle)


#color_v <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
outdir <- "../9.spaAnaly/NewCombined.Analy"
dir.create(outdir)

twogrp <- c("#F5A889", "#ACD6EC")
threegrp <- c("#91C299", "#F5A889", "#ACD6EC")

grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28", "#C85D4D")

heatgradiant <- c("#619DB8", "#AECDD7", "#E3EEEF", "#FAE7D9", "#F0B79A", "#C85D4D")

grad0 <- c("#000000", "#330033", "#660066", "#CC3300", "#FF9900", "#FFFF00")
grad1 <- c("#E5E5E5", "#BBD3BF", "#91C299", "#67B173", "#3DA04D", "#148F28")
grad2 <- c("#E5E5E5", "#E6CDDA", "#E7B6CF", "#E89FC4", "#E988B9", "#EA71AE")
grad3 <- c("#E5E5E5", "#EAD5B7", "#EFC689", "#F4B75B", "#F9A82D", "#FF9900")
grad4 <- c("#E5E5E5", "#C6D9F0", "#9BBCE1", "#709FD2", "#4582C3", "#1A65B4")
grad5 <- c("#E5E5E5", "#D9C2E5", "#C299E5", "#AB70E5", "#9447E5", "#7D1EE5")

grad6 <- c("#440154", "#822681", "#B63679", "#DE4968", "#F1605D", "#FBAE17", "#FDE725")

grad7 <- c("#67000D", "#CB181D", "#FEE0D2", "#DEEBF7", "#2171B5", "#08306B")

resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"

sub_colors <- c("#A8DAB5","#5C9E74","#F1C7DA","#E59CBF","#C973A0")

rescale_to_neg1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1]) * 2 - 1
}

rescale_to_0_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

#######################CNV class for each sample & basal-luminal score, CNV score, monocle3, cytotrace, basal, luminal#######################
##function to make HE image gray
gray_image <- apply(spadat@images$slice1@image[,,1:3], c(1,2), mean)
spadat@images$slice1@image[,,1] <- gray_image
spadat@images$slice1@image[,,2] <- gray_image
spadat@images$slice1@image[,,3] <- gray_image

dir.create("../9.spaAnaly/NewCombined.Analy/03.CNVRes")
cnvres <- "../9.spaAnaly/NewCombined.Analy/03.CNVRes"
#######################allfeatures, spatial, seurat umap#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
bulkcnv <- fread("../BulkRawData/04.GISTIC/all_thresholded.by_genes.txt", header = T, stringsAsFactors = F, data.table = F)

for (sam in usefulsams) {
  sample = allsams_id[allsams_id$ID == sam, ]$MyID
  cnvdir <- paste("../9.spaAnaly/NewCombined.Analy/02.InferCNV/", sam, "/InferCNV02", sep = "")
  seudat <- readRDS(file.path(cnvdir, "visual.score.seurat.rds"))
  
  tmpcnvmtr <- bulkcnv[, paste(sam, ".sorted_markdup", sep = "")]
  names(tmpcnvmtr) <- bulkcnv$`Gene Symbol`
  
  newcnvmtr <- readRDS(file.path(cnvdir, "22_denoise.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
  newcnvmtr <- as.data.frame(t(newcnvmtr@expr.data))
  
  geneorder <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/02.refers/gencode.v25.annotation.uniq.bed",
                     header = FALSE, stringsAsFactors = FALSE, data.table = FALSE)
  geneorder <- geneorder[geneorder$V1 %in% colnames(newcnvmtr), ]
  geneorder <- geneorder[geneorder$V1 %in% names(tmpcnvmtr), ]
  geneorder$V2 <- factor(geneorder$V2, levels = paste0("chr", 1:22))
  
  ordered_mat <- matrix(tmpcnvmtr[geneorder$V1], ncol = 1)
  rownames(ordered_mat) <- geneorder$V1
  colnames(ordered_mat) <- sample
  
  colannot <- data.frame(Chr = geneorder$V2)
  rownames(colannot) <- geneorder$V1
  
  chr_colors <- rep(c("lightgrey", "darkgrey"), length.out = length(levels(geneorder$V2)))
  names(chr_colors) <- levels(geneorder$V2)
  
  top_ann <- HeatmapAnnotation(
    chrs = anno_block(gp = gpar(fill = chr_colors)),
    show_legend = FALSE
  )
  
  png(file.path(cnvres, paste0(sam, ".CNV_Annot.WES.png")), 
      width = 1800, height = 900, res = 300)
  print(Heatmap(t(ordered_mat),
                name = "inferCNV",
                col = structure(
                  c("#377EB8", "#9ecae1", "#FFFFFF", "#fcae91", "#E41A1C"),
                  names = c("-2", "-1", "0", "1", "2")
                ),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                row_names_side = "left",
                show_column_names = FALSE,
                column_split = colannot$Chr,
                top_annotation = top_ann,
                column_title = NULL,
                row_gap = unit(0, "mm"),
                column_gap = unit(0, "mm"),
                border = "gray85",
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 12),
                  legend_height = unit(5, "cm"),
                  title_position = 'topleft'
                )
  ))
  dev.off()
  
  
}





for (sam in usefulsams) {
  cnvdir <- paste("../9.spaAnaly/NewCombined.Analy/02.InferCNV/", sam, "/InferCNV02", sep = "")
  seudat <- readRDS(file.path(cnvdir, "visual.score.seurat.rds"))
  gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
  seudat@images$slice1@image[,,1] <- gray_image
  seudat@images$slice1@image[,,2] <- gray_image
  seudat@images$slice1@image[,,3] <- gray_image
  
  seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
  
  newsize = 8
  p2 <- SpatialDimPlot(seudat, group.by = "seurat_clusters", image.alpha=0.2, stroke = NA, pt.size.factor = newsize) + 
    scale_fill_manual(values = c(grpcolor2, grpcolors))+
    ggtitle("SeuratClusters") +
    theme(
      plot.title = element_text(size = 20),
      legend.position = "none"
    )
  
  p3 <- SpatialDimPlot(seudat, group.by = "Annot_V2", image.alpha=0.2, stroke = NA, pt.size.factor = newsize) + 
    scale_fill_manual(values = grpcolors)+
    ggtitle("CNV") +
    theme(
      plot.title = element_text(size = 20),
      legend.position = "none"
    )
  
  p4 <- SpatialFeaturePlot(seudat, features = "CNVscore", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = c("#E5E5E5", "#E5E5E5", "#D9C2E5", "#C299E5", "#AB70E5", "#AB70E5", "#9447E5", "#9447E5", "#7D1EE5"))+
    ggtitle("CNVscore") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p5 <- SpatialFeaturePlot(seudat, features = "pseudotime", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad4)+
    ggtitle("pseudotime") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p6 <- SpatialFeaturePlot(seudat, features = "cytotrace", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad3)+
    ggtitle("cytotrace") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p7 <- SpatialFeaturePlot(seudat, features = "basal", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = c(grad2))+
    ggtitle("basal") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p8 <- SpatialFeaturePlot(seudat, features = "luminal", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad1)+
    ggtitle("luminal") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  seuptsize = 2
  # p2: UMAP plot for SeuratClusters
  p12 <- DimPlot(seudat, group.by = "seurat_clusters", pt.size = seuptsize) +
    scale_color_manual(values = c(grpcolor2, grpcolors))+
    theme(legend.position = "none")
  
  # p3: UMAP plot for Annot_V2
  p13 <- DimPlot(seudat, group.by = "Annot_V2", pt.size = seuptsize) +
    scale_color_manual(values = grpcolors)+
    theme(legend.position = "none")
  
  # p4: Feature plot for CNVscore
  p14 <- FeaturePlot(seudat, features = "CNVscore", pt.size = seuptsize) +
    scale_color_gradientn(colors = c("#E5E5E5", "#E5E5E5", "#D9C2E5", "#C299E5", "#AB70E5", "#AB70E5", "#9447E5", "#9447E5", "#7D1EE5"))+
    theme(legend.position = "none")
  
  # p5: Feature plot for pseudotime
  p15 <- FeaturePlot(seudat, features = "pseudotime", pt.size = seuptsize) +
    scale_color_gradientn(colors = grad4)+
    theme(legend.position = "none")
  
  # p6: Feature plot for cytotrace
  p16 <- FeaturePlot(seudat, features = "cytotrace", pt.size = seuptsize) +
    scale_color_gradientn(colors = grad3)+
    theme(legend.position = "none")
  
  # p7: Feature plot for basal
  p17 <- FeaturePlot(seudat, features = "basal", pt.size = seuptsize) +
    scale_color_gradientn(colors = c(grad2))+
    theme(legend.position = "none")
  
  # p8: Feature plot for luminal
  p18 <- FeaturePlot(seudat, features = "luminal", pt.size = seuptsize) +
    scale_color_gradientn(colors = grad1)+
    theme(legend.position = "none")
  
  pdf(file.path(cnvres, paste(sam, ".visual.score.seurat.pdf", sep = "")), width = 30, height = 8)
  print((p2|p3|p4|p5|p6|p7|p8)/(p12|p13|p14|p15|p16|p17|p18))
  dev.off()
  
  
  spot_class <- fread(file.path(cnvdir, "CNV_Annot.clusters.tsv"), header = T, stringsAsFactors = F, data.table = F)
  allidsort <- spot_class$ID
  #seudat <- readRDS(file.path(outdir, "seurat.V1.rds"))
  cnvmtr <- readRDS(file.path(cnvdir, "22_denoise.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
  cnvmtr <- as.data.frame(t(cnvmtr@expr.data))
  cnvmtr <- cnvmtr[allidsort, ]
  cnvmtr$Class <- spot_class$Annot_V2
  cnvmtr <- aggregate(.~ Class, cnvmtr, mean)
  rownames(cnvmtr) <- cnvmtr$Class
  cnvmtr <- cnvmtr[, -1]
  #cnvmtr <- cnvmtr[1:200, 1:300]
  
  shif <- 0.05
  newcnvmtr <- cnvmtr
  newcnvmtr[newcnvmtr < (1+shif) & newcnvmtr > (1-shif)] <- 1
  ############heatmap
  geneorder <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/02.refers/gencode.v25.annotation.uniq.bed",
                     header = F, stringsAsFactors = F, data.table = F)
  geneorder <- geneorder[geneorder$V1 %in% colnames(newcnvmtr), ]
  chr_gene_counts <- as.numeric(table(geneorder$V2))
  #geneorder$V2 <- factor(geneorder$V2, c("chr1", "chr2", "chr3"))
  geneorder$V2 <- factor(geneorder$V2, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                         "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                                         "chr21", "chr22"))
  
  gorder <- as.data.frame(table(geneorder$V2))
  gorder$Freq2 <- cumsum(gorder$Freq)
  
  colannot <- data.frame(Chr = geneorder$V2)
  rownames(colannot) <- geneorder$V1
  mat_colors <- list(chr = c(rep(c("lightgrey", "darkgrey"), 11)))
  names(mat_colors$chr) <- unique(colannot$Chr)
  
  chr_colors <- c(rep(c("lightgrey", "darkgrey"), 11))
  
  # Prepare the top annotation
  top_ann <- HeatmapAnnotation(
    chrs = anno_block(
      gp = gpar(fill = chr_colors)
    ),
    show_legend = FALSE
  )
  
  
  pdf(file.path(cnvres, paste(sam, ".CNV_Annot.clusters.pdf", sep = "")), width = 6, height = 3)
  print(Heatmap(as.matrix(newcnvmtr[,geneorder$V1]), name = "inferCNV",
                col = colorRamp2(c(0.93,1,1.07), c("#377EB8","#FFFFFF","#E41A1C")),
                cluster_rows = F,
                cluster_columns = F,
                show_row_names=T,
                row_names_side = "left",
                show_column_names = F,
                #column_names_side = "bottom",
                show_heatmap_legend = T,
                top_annotation = top_ann,
                #left_annotation = left_anno,
                #right_annotation = rowAnnotation(cnv.score = anno_points(cluster.cnv.score[rownames(avg.cnv.df),'cnv.score'])),
                row_split = rownames(newcnvmtr),
                #row_title = NULL,
                column_split = colannot$Chr,
                column_title = NULL,
                row_gap = unit(0, "mm"),
                column_gap = unit(0, "mm"),
                border ="gray85",
                heatmap_legend_param = list(
                  #at = c(0.9, 0, 1.1),
                  #labels = c("low", "zero", "high"),
                  #title = "Some values",
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 12),
                  legend_height = unit(5, "cm"),
                  title_position ='topleft'
                )))
  dev.off()
  
  library(data.tree)
  library(phangorn)
  library(patchwork)  # For arranging multiple plots
  library(data.table)
  library(ggtree)
  
  methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  plots <- list()  # Store plots
  #length(mymtapinfo$ID)
  spot_class <- fread(file.path(cnvdir, "CNV_Annot.clusters.tsv"), header = T, stringsAsFactors = F, data.table = F)
  spot_class <- spot_class[order(spot_class$Annot_V2), ]
  allidsort <- spot_class$ID
  
  cnvmtr <- readRDS(file.path(cnvdir, "22_denoise.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
  cnvmtr <- as.data.frame(t(cnvmtr@expr.data))
  cnvmtr <- cnvmtr[allidsort, ]
  cnvmtr$Class <- spot_class$Annot_V2
  cnvmtr <- aggregate(.~ Class, cnvmtr, mean)
  rownames(cnvmtr) <- cnvmtr$Class
  cnvmtr <- cnvmtr[, -1]
  #cnvmtr <- cnvmtr[1:200, 1:300]
  
  tmpcnvmtr <- cnvmtr
  tmpcnvmtr[tmpcnvmtr < (1-shif)] <- (1-shif)
  tmpcnvmtr[tmpcnvmtr > (1+shif)] <- (1+shif)
  for (method in methods) {
    exp_cluster_treeNJ <- NJ(dist(tmpcnvmtr, method = method))
    
    meta.data <- data.frame(group = exp_cluster_treeNJ$tip.label, 
                            group2 = exp_cluster_treeNJ$tip.label)
    rownames(meta.data) <- exp_cluster_treeNJ$tip.label
    
    groupInfo <- split(rownames(meta.data), meta.data$group)
    treeB2 <- groupOTU(exp_cluster_treeNJ, groupInfo)
    
    gp <- ggtree(treeB2, aes(color = group), ladderize = TRUE, size = 4) +
      theme_tree() +
      geom_tippoint(alpha = 1.0, size = 10) +
      ggtitle(paste("Method:", method))+
      scale_color_manual(values = c("tum_c1" = "#C05050", "tum_c2" = "#80A0C5", "Urothelium" = "gray80"))
    
    plots[[method]] <- gp
  }
  # Combine all trees into one figure
  combined_plot <- wrap_plots(plots, ncol = 3)
  
  pdf(file.path(cnvres, paste(sam, ".CNV_Annot.clusters.tree.pdf", sep = "")), width = 15, height = 10)
  print(combined_plot)
  dev.off()
  
  seudat <- readRDS(file.path(cnvdir, "visual.score.seurat.rds"))
  seudat$Class <- ifelse(seudat$basal > seudat$luminal, "Basal", "Luminal")
  
  library(ggpubr)
  library(ggplot2)
  library(dplyr)
  
  # Create the table
  tbl <- as.data.frame(table(seudat$Class, seudat$Annot_V2))
  
  # Create pie plot for tum_c1
  tum_c1 <- tbl %>% 
    filter(Var2 == "tum_c1") %>% 
    ggpie(x = "Freq", label = "Var1", 
          fill = "Var1", 
          color = "white", 
          palette = c("#EA71AE", "#148F28")) +
    ggtitle("Tum_c1")+
    theme_void()
  
  # Create pie plot for tum_c2
  tum_c2 <- tbl %>% 
    filter(Var2 == "tum_c2") %>% 
    ggpie(x = "Freq", label = "Var1", 
          fill = "Var1", 
          color = "white", 
          palette = c("#EA71AE", "#148F28")) +
    ggtitle("Tum_c2")+
    theme_void()
  
  
  pdf(file.path(cnvres, paste(sam, ".pieplot.pdf", sep = "")), width = 5, height = 3)
  print(ggarrange(tum_c1, tum_c2, ncol = 2))
  dev.off()
  
  spot_class <- fread(file.path(cnvdir, "CNV_Annot.clusters.tsv"), header = T, stringsAsFactors = F, data.table = F)
  spot_class <- spot_class[spot_class$Annot_V1 == "Tumor", ]
  spot_class <- spot_class[order(spot_class$Annot_V2), ]
  allidsort <- spot_class$ID
  
  cnvmtr <- readRDS(file.path(cnvdir, "22_denoise.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
  cnvmtr <- as.data.frame(t(cnvmtr@expr.data))
  cnvmtr <- cnvmtr[allidsort, ]
  
  sample_map <- spot_class %>% select(Sample = ID, Class = Annot_V2)
  unique_classes <- unique(sample_map$Class)
  
  cnvmtr1 <- cnvmtr[sample_map[sample_map$Class == "tum_c1", ]$Sample, ]
  cnvmtr2 <- cnvmtr[sample_map[sample_map$Class == "tum_c2", ]$Sample, ]
  
  cnvmtr <- rbind(cnvmtr1, cnvmtr2)
  
  deval <- 0.1
  amp_threshold <- 1 + deval
  del_threshold <- 1 - deval
  allclass_cnv <- data.frame()
  for (class_name in unique_classes) {
    
    # Subset CNV matrix by class
    class_cells <- sample_map %>%
      filter(Class == class_name) %>%
      pull(Sample)
    
    class_cnv_matrix <- cnvmtr[rownames(cnvmtr) %in% class_cells, , drop = FALSE]
    
    # Compute CNV frequencies for each gene
    cnv_stats <- data.frame(Gene = colnames(class_cnv_matrix))
    
    compute_cnv_frequencies <- function(gene) {
      gene_values <- class_cnv_matrix[, gene, drop = FALSE]
      amp_freq <- sum(gene_values > amp_threshold) / nrow(gene_values)
      del_freq <- sum(gene_values < del_threshold) / nrow(gene_values)
      return(c(amp_freq, del_freq))
    }
    
    cnv_stats[, 2:3] <- t(sapply(colnames(class_cnv_matrix), compute_cnv_frequencies))
    colnames(cnv_stats) <- c("Gene", "Amp_Freq", "Del_Freq")
    
    # Merge CNV stats with geneorder to get chromosome positions
    merged_data <- merge(cnv_stats, geneorder, by.x = "Gene", by.y = "V1")
    merged_data$Chromosome <- as.numeric(gsub("chr", "", merged_data$V2))
    
    # Assign chromosome positions based on gene counts
    merged_data <- merged_data[order(merged_data$Chromosome), ]
    
    # Assign start positions using chromosome lengths from gene counts
    merged_data$StartPos <- match(merged_data$Gene, colnames(class_cnv_matrix))
    
    # Smooth CNV frequencies using spline
    spline_amp <- smooth.spline(merged_data$StartPos, merged_data$Amp_Freq, spar = 0.1)
    spline_del <- smooth.spline(merged_data$StartPos, merged_data$Del_Freq, spar = 0.1)
    
    merged_data$Amp_Freq_smooth <- predict(spline_amp, merged_data$StartPos)$y
    merged_data$Del_Freq_smooth <- predict(spline_del, merged_data$StartPos)$y
    
    # Add class information
    merged_data$Class <- class_name
    
    # Append to final dataframe
    allclass_cnv <- rbind(allclass_cnv, merged_data)
  }
  df <- data.frame(
    chromName = levels(geneorder$V2),
    chromlength = gorder$Freq,
    chromNum = 1:length(levels(geneorder$V2))
  )
  
  df$chromlengthCumsum <- cumsum(df$chromlength)
  df$chromStartPosFrom0 <- c(0, df$chromlengthCumsum[-nrow(df)])
  df$chromMidelePosFrom0 <- df$chromStartPosFrom0 + (df$chromlength / 2)
  df$ypos <- rep(c(0.4, 0.45), length.out = nrow(df))
  
  allclass_cnv$Amp_Freq_smooth <- ifelse(allclass_cnv$Amp_Freq_smooth < 0.05, 0, allclass_cnv$Amp_Freq_smooth)
  allclass_cnv$Del_Freq_smooth <- ifelse(allclass_cnv$Del_Freq_smooth < 0.05, 0, allclass_cnv$Del_Freq_smooth)
  p2 <- ggplot(allclass_cnv, aes(StartPos)) +
    geom_ribbon(aes(ymin = 0, ymax = Amp_Freq_smooth), fill = "red", alpha = 0.4) +
    geom_ribbon(aes(ymin = 0, ymax = -Del_Freq_smooth), fill = "blue", alpha = 0.4) +
    geom_vline(data = df, aes(xintercept = chromlengthCumsum), linetype = 2, color = "gray85") +
    
    # Set consistent x-axis limits based on gene counts
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(df$chromlengthCumsum))) + 
    
    ylim(-0.7, 0.7) +
    
    # Remove gray grid lines and add a black border
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),     # Remove major grid lines
      panel.grid.minor = element_blank(),     # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 0.5), # Add black border
      strip.background = element_blank(),     # Remove background of facet labels
      axis.line = element_blank()             # Remove axis lines
    ) +
    
    facet_wrap(~ Class, ncol = 1, scales = "free_y")
  
  # Save to PDF
  pdf(file.path(cnvres, paste(sam, ".CNV_Annot.clusters.freq.pdf", sep = "")), width = 6, height = 4)
  print(p2)
  dev.off()
  
}

#######################score changes across all samples#######################

allres <- data.frame()
for (sam in usefulsams) {
  cnvdir <- paste("../9.spaAnaly/NewCombined.Analy/02.InferCNV/", sam, "/InferCNV02", sep = "")
  seudat <- readRDS(file.path(cnvdir, "visual.score.seurat.rds"))
  gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
  seudat@images$slice1@image[,,1] <- gray_image
  seudat@images$slice1@image[,,2] <- gray_image
  seudat@images$slice1@image[,,3] <- gray_image
  
  seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
  
  tmpres <- seudat@meta.data[, c("Annot_V2", "CNVscore", "pseudotime", "cytotrace", "basal", "luminal")]
  normalized_data <- tmpres %>%
    group_by(Annot_V2) %>%
    summarise(Mean_CNV = median(CNVscore, na.rm = TRUE),
              Mean_Trajec = median(pseudotime, na.rm = TRUE),
              Mean_Cyto = median(cytotrace, na.rm = TRUE),
              Mean_Basal = median(basal, na.rm = TRUE),
              Mean_Lum = median(luminal, na.rm = TRUE)) %>%
    ungroup()
  normalized_data$Sample = sam
  allres <- rbind(allres, normalized_data)
}


allsams_id <- fread(file.path(cnvres, "cnv.clone.lineage.info"), header = T, data.table = F, stringsAsFactors = F)
allsams_id$Merged_ID <- paste(allsams_id$Sample, allsams_id$Group, sep = ":")
allres$Merged_ID <- paste(allres$Sample, allres$Annot_V2, sep = ":")

allsaminfo <- merge(allres, allsams_id, by.x = "Merged_ID", by.y = "Merged_ID")

res <- allsaminfo[, c("MyID", "Lineage", "Mean_CNV", "Mean_Trajec", "Mean_Cyto", "Mean_Basal", "Mean_Lum")]

res <- melt(res, id.vars = c("MyID", "Lineage"))

my_comp <- list(c("Basal", "Luminal"))
p <- ggplot(res, aes(x = Lineage, y = value, group = MyID)) +
  geom_point(aes(color = Lineage), size = 3) +  # Add points
  geom_line(aes(group = MyID), size = 1) +  # Connect dots for each sample
  #stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  stat_compare_means(comparisons = my_comp, method = "wilcox", paired = T)+
  scale_color_manual(values = c("Basal" = "#EA71AE", "Luminal" = "#148F28"))+
  facet_wrap( ~ variable, scales = "free_y", ncol = 5)+
  theme_classic2()


ggsave(file.path(cnvres, "V2503.figure2.infercnv.combined.scores.pdf"), p, width = 10, height = 3)



#######################DEG and pathway differences in basal-luminal clones#######################
resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
combined_seurat <- readRDS(file.path(resdir, "allsam.merged.seurat.umap.rds"))

allsaminfo <- data.frame()
for (j in 1:length(usefulsams)) {
  print(j)
  sample = usefulsams[j]
  
  cnvdir <- paste("../9.spaAnaly/NewCombined.Analy/02.InferCNV/", sample, "/InferCNV02", sep = "")
  seudat <- readRDS(file.path(cnvdir, "visual.score.seurat.rds"))
  seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
  allsaminfo <- rbind(allsaminfo, seudat@meta.data)
  
}

allsams_id <- fread(file.path(cnvres, "cnv.clone.lineage.info"), header = T, data.table = F, stringsAsFactors = F)
allsams_id$Merged_ID <- paste(allsams_id$MyID, allsams_id$Group, sep = ":")
allsaminfo$Merged_ID <- paste(allsaminfo$BCID, allsaminfo$Annot_V2, sep = ":")

allsaminfo <- merge(allsaminfo, allsams_id, by.x = "Merged_ID", by.y = "Merged_ID")

combined_seurat <- subset(combined_seurat, subset = CellID %in% allsaminfo$CellID)

newinfo <- allsaminfo[, c("CellID", "Lineage")]
rownames(newinfo) <- newinfo$CellID

combined_seurat <- AddMetaData(combined_seurat, newinfo)

Idents(combined_seurat) <- combined_seurat$Lineage
markers <- FindAllMarkers(combined_seurat, only.pos = FALSE, min.pct = 0.0, logfc.threshold = 0.0)
markers <- markers[markers$cluster == "Basal", ]

library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(geneList)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
hall_kegg <- rbind(kegg, hall)
hall_kegg <- hall_kegg[gsub(".*_", "", hall_kegg$gs_name) != "DISEASE", ]


gene_list <- markers$avg_log2FC
names(gene_list) <- markers$gene
gene_list <- sort(gene_list, decreasing = TRUE)

path_dat <- GSEA(gene_list, TERM2GENE = hall_kegg, pvalueCutoff = 1)


p <- enrichmentNetwork(path_dat@result, drawEllipses = TRUE, fontSize = 2.5)

p

pdf(file.path(cnvres, "V2502.figure2.infercnv.gene.fishertest.basal2luminal.pathways.pdf"), width = 12, height = 6)
print(p)
dev.off()





