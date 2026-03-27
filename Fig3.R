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

#############################################################figure3 part#############################################################
dir.create("../9.spaAnaly/NewCombined.Analy/04.subtype")
subres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"

#######################spatial epithelial umap, markers, signatures, trajectory, genes with trajectory#######################
spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
Idents(spadat) <- spadat$BLClass
markers <- FindAllMarkers(spadat, only.pos = TRUE, min.pct = 0.0, logfc.threshold = 0.1)
top50 <- top_n(group_by(markers,cluster), n = 50, wt = avg_log2FC)

#############umap
ptsize = 4
pdf(file.path(subres, "V2503.figure3.spatial.class.umap.pdf"), width = 4.7, height = 4)
print(DimPlot(spadat, reduction = "umap.ischia12", label = TRUE, group.by = "BLClass", pt.size = ptsize, raster = T)+
        scale_color_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0")))
dev.off()
pdf(file.path(subres, "V2503.figure3.spatial.class.pseudotime.pdf"), width = 4.5, height = 4)
print(FeaturePlot(spadat, features = "pseudotime", pt.size = ptsize, raster = T) +
        scale_color_gradientn(colors = grad6))
dev.off()

pdf(file.path(subres, "V2503.figure3.spatial.class.tmp.pdf"), width = 4.5, height = 4)
print(FeaturePlot(spadat, features = "EGFR", pt.size = ptsize, raster = T) +
        scale_color_gradientn(colors = grad6))
print(FeaturePlot(spadat, features = "FGFR3", pt.size = ptsize, raster = T) +
        scale_color_gradientn(colors = grad6))
dev.off()


#############markers
Idents(spadat) <- spadat$BLClass
selgenes <- c("UPK3A", "GATA3", "CLDN4", "ERBB2", "FGFR3",
              "CCND1", "KRT20", "PPARG", "EPHX3",
              "FABP4", "VEGFA", "EGFR", "TGFB1", 
              "KRT6A", "TP63", "KRT5", "CXCL10", "SPARC", "WNT5A", "MMP11",
              "HIF1A", "TGFB2", "CD274", "TIGIT", "IL10", "MMP9", "PIK3CA",
              "IFNG", "IFNA1")
Idents(spadat) <- spadat$BLClass
avg_exp <- AverageExpression(spadat, features = selgenes, return.seurat = FALSE)$RNA
# Optional: scale expression per gene (row)
scaled_matrix <- t(scale(t(avg_exp)))  # z-score scaling by row
scaled_matrix <- scaled_matrix[selgenes, c("Lum", "IM", "B.EMT", "B.IS")]
# Set output file
pdf(file.path(subres, "V2503.figure3.spatial.class.markergene.pdf"), width = 8, height = 2)
# Draw heatmap
ht <- Heatmap(t(scaled_matrix),
              name = "Scaled\nExpression",
              col = colorRamp2(c(-2, 0.1, 1), c("#619DB8", "#E3EEEF", "#C85D4D")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10),
              show_column_dend = FALSE,
              show_row_dend = FALSE)
draw(ht)
decorate_heatmap_body("Scaled\nExpression", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1.5))
})
dev.off()


#############fraction along trajectory
data <- data.frame(
  pseudotime = spadat$pseudotime,
  BLClass = spadat$BLClass
) %>% filter(!is.na(BLClass))

data <- data %>% arrange(pseudotime)
data$group <- ceiling(seq_along(data$pseudotime) / 10)

proportion_data <- data %>%
  group_by(group, BLClass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  complete(group, BLClass, fill = list(proportion = 0))

proportion_data$relative_dist <- proportion_data$group / max(proportion_data$group)

proportion_data <- proportion_data %>%
  group_by(BLClass) %>%
  arrange(relative_dist) %>%
  mutate(proportion_smooth = predict(loess(proportion ~ relative_dist, span = 0.3))) %>%
  ungroup() %>%
  mutate(proportion_smooth = pmax(0, proportion_smooth))

proportion_data <- proportion_data %>%
  group_by(group) %>%
  mutate(
    sum_prop = sum(proportion_smooth),
    proportion_smooth = ifelse(sum_prop == 0, 0, proportion_smooth / sum_prop)
  ) %>%
  select(-sum_prop) %>%
  ungroup()

# Step 6: area
p_area <- ggplot(proportion_data, aes(x = relative_dist, y = proportion_smooth, fill = BLClass)) +
  geom_area(position = "stack", color = "black", size = 0.2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Relative distance", y = "Proportion") +
  scale_fill_manual(values = c(
    Lum = "#5C9E74",
    IM = "#A8DAB5",
    B.EMT = "#F1C7DA",
    B.IS = "#C973A0"
  )) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  )

pdf(file.path(subres, "V2503.figure3.spatial.class.proportion.trajectory.area.pdf"), width = 10, height = 4)
print(p_area)
dev.off()




#############gene changes along trajectory
library(monocle)
selgenes <- c("UPK3A", "GATA3", "CLDN4", "ERBB2", "FGFR3",
              "CCND1", "KRT20", "PPARG", "EPHX3",
              "PLK1", "FABP4", "VEGFA", "EGFR", "TGFB1", 
              "KRT6A", "KRT5", "TP63", "CXCL10", "SPARC", "WNT5A", "MMP11",
              "HIF1A", "CD274", "TIGIT", "IL10", "MMP9", "TGFB2", "PIK3CA",
              "IFNG", "IFNA1")

data <- data.frame(
  pseudotime = spadat$pseudotime,
  BLClass = spadat$BLClass
) %>% filter(!is.na(BLClass))

cds <- readRDS(file.path(subres, "V2503.figure3.allspots.basal.luminal.monocle2.rds"))
cds$Pseudotime <- data[colnames(cds), "pseudotime"]
cds$BLClass <- data[colnames(cds), "BLClass"]


cds_subset <- cds[selgenes, ]
newdata <- data.frame(Pseudotime = seq(min(cds_subset$Pseudotime), max(cds_subset$Pseudotime),length.out = 100)) 

m <- genSmoothCurves(cds_subset, cores=8, trend_formula = '~sm.ns(Pseudotime, df=3)',  
                     relative_expr = T, new_data = newdata)

#remove genes with no expression in any condition
m=m[!apply(m,1,sum)==0,]
m = vstExprs(cds_subset, expr_matrix=m)
# Row-center the data.
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
m[is.nan(m)] = 0
m[m>3] = 3
m[m<-3] = -3

heatmap_matrix <- m

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

library(colorRamps)
bks <- seq(-2.1,2.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)

pdf(file.path(subres, "V2503.figure3.spatial.class.gene.trajectory.pdf"), width = 8, height = 8)
print(pheatmap(heatmap_matrix, 
               cluster_cols=FALSE, 
               cluster_rows=FALSE, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               silent=TRUE,
               breaks=bks,
               border_color = NA,
               color=hmcols))

dev.off()


#############density plot
spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
bclineage <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/bc.linerge.txt",
                   header = T, stringsAsFactors = F, data.table = F)
bclineage$Signature <- paste("TCGA.", bclineage$Signature, sep = "")
names(bclineage) <- c("Signature", "Gene")
allscores <- split(bclineage$Gene, bclineage$Signature)

spadat <- AddModuleScore(
  object = spadat,
  features = allscores,
  nbin = 12,
  name = "TCELLSIGNATURE"
)

colnames(spadat@meta.data)[grep("TCELLSIGNATURE", colnames(spadat@meta.data))] <- names(allscores)

spadat$pseudotime <- rescale_to_neg1_1(spadat$pseudotime)
spadat$TCGA.EMT <- rescale_to_neg1_1(spadat$TCGA.EMT)

tmpresres <- spadat@meta.data[, c("BLClass", "tmp", "luminal")]
names(tmpresres) <- c("BLClass", "basal", "luminal")

res <- tmpresres[, c("BLClass", "basal", "luminal")]
res$basal <- rescale_to_neg1_1(res$basal)
res$luminal <- rescale_to_neg1_1(res$luminal)

names(res) <- c("BLClass", "basal", "luminal")
res$BLClass <- factor(res$BLClass, c("Lum", "IM", "B.EMT", "B.IS"))

a <- quantile(res$basal,probs = seq(0, 1, by=0.01), na.rm=T)
res$basal[res$basal>a[[length(a)-1]]] <- a[[length(a)-1]];
res$basal[res$basal< a[[2]]] <- a[[2]];

a <- quantile(res$luminal,probs = seq(0, 1, by=0.01), na.rm=T)
res$luminal[res$luminal>a[[length(a)-1]]] <- a[[length(a)-1]];
res$luminal[res$luminal< a[[2]]] <- a[[2]];


res$basal <- rescale_to_neg1_1(res$basal)
res$luminal <- rescale_to_neg1_1(res$luminal)

minx <- -1
miny <- -1
maxx <- 1
maxy <- 1

tmpres1 <- res[res$BLClass == "Lum", ]
tmpres2 <- res[res$BLClass == "IM", ]
p1 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#5C9E74") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "IM", ]
tmpres2 <- res[res$BLClass == "IM", ]
p2 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#A8DAB5") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "B.EMT", ]
tmpres2 <- res[res$BLClass == "IM", ]
p3 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#F1C7DA") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "B.IS", ]
tmpres2 <- res[res$BLClass == "IM", ]
p4 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#C973A0") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(subres, "V2503.figure3.spatial.class.density.pdf"), width = 6, height = 4)
print((p1|p2)/(p3|p4))
dev.off()



#######################single cell epithelial umap, trajectory, genes with trajectory, density plot#######################
spadat <- readRDS(file.path(subres, "V2503.figure3.scdata.basal.luminal.classification.final.rds"))
#############umap
ptsize = 6
pdf(file.path(subres, "V2503.figure3.scdata.class.umap.pdf"), width = 4.7, height = 4)
print(DimPlot(spadat, reduction = "umap.ischia12", label = TRUE, group.by = "BLClass", pt.size = ptsize, raster = T)+
        scale_color_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0")))
dev.off()
pdf(file.path(subres, "V2503.figure3.scdata.class.pseudotime.pdf"), width = 4.5, height = 4)
print(FeaturePlot(spadat, reduction = "umap.ischia12", features = "pseudotime", pt.size = ptsize, raster = T) +
        scale_color_gradientn(colors = grad6))
dev.off()

#############genes

selgenes <- c("UPK3A", "GATA3", "CLDN4", "ERBB2", "FGFR3",
              "CCND1", "KRT20", "PPARG", "EPHX3",
              "VEGFA", "FABP4", "EGFR", 
              "KRT6A", "TP63", "KRT5", "SPARC", "WNT5A","IFNA1", "MMP11",
              "CD274", "MMP9", "PIK3CA", "SNAI1", "VIM", "ITGA5", "IL10",
              "IFNG", "TGFB1", "CXCL10", "TIGIT")
Idents(spadat) <- spadat$BLClass
avg_exp <- AverageExpression(spadat, features = selgenes, return.seurat = FALSE)$RNA
# Optional: scale expression per gene (row)
scaled_matrix <- t(scale(t(avg_exp)))  # z-score scaling by row
scaled_matrix <- scaled_matrix[selgenes, c("Lum", "IM", "B.EMT", "B.IS")]
# Set output file
pdf(file.path(subres, "V2503.figure3.scdata.class.markergene.pdf"), width = 8, height = 2)
# Draw heatmap
ht <- Heatmap(t(scaled_matrix),
              name = "Scaled\nExpression",
              col = colorRamp2(c(-2, 0.4, 1), c("#619DB8", "#E3EEEF", "#C85D4D")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10),
              show_column_dend = FALSE,
              show_row_dend = FALSE)
draw(ht)
decorate_heatmap_body("Scaled\nExpression", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1.5))
})
dev.off()



#############density plot
bclineage <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/bc.linerge.txt",
                   header = T, stringsAsFactors = F, data.table = F)
bclineage$Signature <- paste("TCGA.", bclineage$Signature, sep = "")
names(bclineage) <- c("Signature", "Gene")

allscores <- split(bclineage$Gene, bclineage$Signature)
spadat <- AddModuleScore(
  object = spadat,
  features = allscores,
  nbin = 12,
  name = "TCELLSIGNATURE"
)
# Rename columns to proper signature names
colnames(spadat@meta.data)[grep("TCELLSIGNATURE", colnames(spadat@meta.data))] <- names(allscores)

spadat$pseudotime <- rescale_to_neg1_1(spadat$pseudotime)
spadat$TCGA.basal <- rescale_to_neg1_1(spadat$TCGA.basal)

res <- spadat@meta.data[, c("BLClass", "tmp", "luminal")]
names(res) <- c("BLClass", "basal", "luminal")
res$BLClass <- factor(res$BLClass, c("Lum", "IM", "B.EMT", "B.IS"))

a <- quantile(res$basal,probs = seq(0, 1, by=0.005))
res$basal[res$basal>a[[length(a)-1]]] <- a[[length(a)-1]];
res$basal[res$basal< a[[2]]] <- a[[2]];

a <- quantile(res$luminal,probs = seq(0, 1, by=0.005))
res$luminal[res$luminal>a[[length(a)-1]]] <- a[[length(a)-1]];
res$luminal[res$luminal< a[[2]]] <- a[[2]];


rescale_to_neg1_1 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1]) * 2 - 1
}
res$basal <- rescale_to_neg1_1(res$basal)
res$luminal <- rescale_to_neg1_1(res$luminal)

minx <- -1
miny <- -1
maxx <- 1
maxy <- 1

tmpres1 <- res[res$BLClass == "Lum", ]
tmpres2 <- res[res$BLClass == "IM", ]
p1 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#5C9E74") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "IM", ]
tmpres2 <- res[res$BLClass == "IM", ]
p2 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#A8DAB5") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "B.EMT", ]
tmpres2 <- res[res$BLClass == "IM", ]
p3 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#F1C7DA") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

tmpres1 <- res[res$BLClass == "B.IS", ]
tmpres2 <- res[res$BLClass == "IM", ]
p4 <- ggplot(tmpres1, aes(x = luminal, y = basal)) +
  geom_density_2d(data = tmpres2, aes(x = luminal, y = basal), color = "gray80")+
  geom_density_2d(color = "#C973A0") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  xlim(miny, maxy) +
  ylim(minx, maxx) +
  theme_minimal() +
  coord_fixed() +
  labs(x = "Luminal Z-score", y = "Basal Z-score") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(subres, "V2503.figure3.scdata.class.density.pdf"), width = 6, height = 4)
print((p1|p2)/(p3|p4))
dev.off()






#######################pathway enriched between different groups#######################
spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))

Idents(spadat) <- spadat$BLClass
markers <- FindAllMarkers(spadat, only.pos = FALSE, min.pct = 0.0, logfc.threshold = 0.0)

hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
names(hall) <- c("Signature", "Genes")
names(c2) <- c("Signature", "Genes")
hall_c2 <- rbind(hall, c2)


library(clusterProfiler)
library(tibble)

gsea_results <- list()
markers$cluster <- factor(markers$cluster, c("Lum", "IM", "B.EMT", "B.IS"))

for (cl in unique(markers$cluster)) {
  cat("Running GSEA for cluster:", cl, "\n")
  
  # Prepare gene list for current cluster
  gene_list <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC)) %>%
    distinct(gene, .keep_all = TRUE)
  
  # Create named vector
  gene_vector <- setNames(gene_list$avg_log2FC, gene_list$gene)
  
  # Run GSEA
  gsea <- GSEA(gene_vector, TERM2GENE = hall_c2, verbose = FALSE, pvalueCutoff = 1)
  
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    gsea_results[[cl]] <- as.data.frame(gsea) %>%
      mutate(cluster = cl)
  }
}

library(tidyr)

# Combine all results
gsea_all <- bind_rows(gsea_results)

# Filter for adjusted p-value < 0.05 and reshape
sig_gsea <- gsea_all %>% 
  filter(NES > 0 & p.adjust < 0.05)

gsea_heat <- gsea_all[gsea_all$ID %in% sig_gsea$ID, ] %>% 
  select(cluster, ID, NES) %>%
  pivot_wider(names_from = cluster, values_from = NES, values_fill = 0)

write.table(gsea_all, file.path(subres, "V2503.figure3.spadat.pathway.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)
write.table(gsea_heat, file.path(subres, "V2503.figure3.spadat.pathway.heatmap.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)

res4heat <- fread(file.path(subres, "V2503.figure3.spadat.pathway.heatmap.selected.txt"),
                  header = T, stringsAsFactors = F, data.table = F)
selgsea_all <- gsea_all[gsea_all$ID %in% res4heat$ID, ]

plotdata <- selgsea_all %>%
  mutate(sig = ifelse(p.adjust < 0.05 & NES > 0, "*", "")) %>%
  mutate(ID = factor(ID, levels = rev(res4heat$ID)))

plotdata$cluster <- factor(plotdata$cluster, c("Lum", "IM", "B.EMT", "B.IS"))

plotdata$NES[plotdata$NES == 0] <- 0.2
plotdata$NES[plotdata$NES > 2] <- 2
plotdata$NES[plotdata$NES < -2] <- -2

plotdata <- plotdata %>%
  complete(cluster, ID, fill = list(NES = -0.3, sig = ""))
p <- ggplot(plotdata, aes(x = cluster, y = ID, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), color = "black", size = 5, vjust = 0.8) +
  scale_fill_gradient2(low = "#619DB8", mid = "#E3EEEF", high = "#C85D4D", 
                       midpoint = 0, name = "NES") +
  theme_minimal() +
  theme(#axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank())

pdf(file.path(subres, "V2503.figure3.spadat.pathway.heatmap.pdf"), width = 7, height = 6)
print(p)
dev.off()


#######################pick overlapped mid high genes#######################
spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
Idents(spadat) <- spadat$BLClass
markers1 <- FindAllMarkers(spadat, only.pos = TRUE, min.pct = 0.0, logfc.threshold = 0.0)


seudat <- readRDS(file.path(subres, "V2503.figure3.scdata.basal.luminal.classification.final.rds"))
Idents(seudat) <- seudat$BLClass
markers2 <- FindAllMarkers(seudat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.1)

cutoff <- 0.7
IM_ovep <- intersect(markers1[markers1$cluster == "Lum" & markers1$p_val_adj < 0.05 & markers1$avg_log2FC > cutoff, ]$gene,
                     markers2[markers2$cluster == "Lum" & markers2$p_val_adj < 0.05 & markers2$avg_log2FC > cutoff, ]$gene)
EMT_ovep <- intersect(markers1[markers1$cluster == "IM" & markers1$p_val_adj < 0.05 & markers1$avg_log2FC > cutoff, ]$gene,
                      markers2[markers2$cluster == "IM" & markers2$p_val_adj < 0.05 & markers2$avg_log2FC > cutoff, ]$gene)

fdaapprove <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/06.HPA/protein_class_FDA.tsv",
                    header = T, stringsAsFactors = F, data.table = F)
fdagenes <- unique(c(fdaapprove$Gene, "NECTIN4"))


length(markers1[markers1$cluster == "IM" & markers1$p_val_adj < 0.05 & markers1$avg_log2FC > cutoff2, ]$gene)
length(markers2[markers2$cluster == "IM" & markers2$p_val_adj < 0.05 & markers2$avg_log2FC > cutoff2, ]$gene)
length(intersect(fdagenes, EMT_ovep))
intersect(fdagenes, EMT_ovep)

a <- markers2[markers2$cluster == "IM" & markers2$p_val_adj < 0.05 & markers2$avg_log2FC > cutoff2, ]$gene
b <- markers1[markers1$cluster == "IM" & markers1$p_val_adj < 0.05 & markers1$avg_log2FC > cutoff2, ]$gene
c <- fdagenes

length(setdiff(a, c(b, c)))
length(setdiff(intersect(a, b), c))
length(setdiff(b, c(a, c)))
length(setdiff(intersect(a, c), b))
length(intersect(intersect(a, c), b))
length(setdiff(intersect(b, c), a))
length(setdiff(c, c(a, b)))


#######################EGFR and cell dependency#######################
cellinfo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/04.CCLE/sample_info.csv",
                  header = T, stringsAsFactors = F, data.table = F, sep = ",")
rownames(cellinfo) <- cellinfo$DepMap_ID
cellinfo <- cellinfo[cellinfo$primary_disease == "Bladder Cancer", ]

ccleexp <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/04.CCLE/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv",
                 header = T, stringsAsFactors = F, data.table = F, sep = ",")
rownames(ccleexp) <- ccleexp$V1
ccleexp <- ccleexp[, -1]
ccleexp <- t(ccleexp)
rownames(ccleexp) <- gsub(" .*", "", rownames(ccleexp))
ccleexp[1:5, 1:5]


ccledat <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/04.CCLE/CRISPRGeneEffect.csv",
                 header = T, stringsAsFactors = F, data.table = F, sep = ",")
rownames(ccledat) <- ccledat$V1
ccledat <- ccledat[, -1]
ccledat <- t(ccledat)
rownames(ccledat) <- gsub(" .*", "", rownames(ccledat))
ccledat[1:5, 1:5]

comsams <- intersect(cellinfo$DepMap_ID, intersect(colnames(ccleexp), colnames(ccledat)))

cellinfo <- cellinfo[comsams, ]
ccleexp <- ccleexp[, comsams]
ccledat <- ccledat[, comsams]

hall <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/bc.linerge.txt", 
              header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
allscores <- split(hall$Genes, hall$Signature)
names(allscores) <- paste(names(allscores), "_Signature", sep = "")

sig_tfs <- gsva(as.matrix(ccleexp), allscores, kcdf="Gaussian",method = "gsva",parallel.sz=10)


refdat <- data.frame(Sample = colnames(exp_data),
                     Class = classification,
                     basal_score = as.numeric(sig_tfs["basal_Signature", ]),
                     EGFR = as.numeric(ccledat["EGFR", ]),
                     FGFR3 = as.numeric(ccledat["FGFR3", ]),
                     Celltype = cellinfo$primary_disease,
                     Type = cellinfo$primary_or_metastasis)


p <- ggplot(data = refdat, mapping = aes(x = basal_score, y = EGFR)) +
  geom_point(shape = 21, fill = "gray50", color = 'gray50', size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E59CBF", fill = "#E59CBF", alpha = 0.2) +
  sm_statCorr(color = "#E59CBF")+
  ylab("EGFR CERES score")+
  xlab("Basal activity")

ggsave(file.path(subres, "V2503.figure3.ccle.gene.effects.pdf"), p, width = 4, height = 3.5)

#######################spatial distribution in each sample#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data
usefulsams <- unique(spadat$BCID)

pdf(file.path(subres, "V2503.figure3.spatial.examples.pdf"), width = 16, height = 4)
for (sam in usefulsams) {
  
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
  seudat@images$slice1@image[,,1] <- gray_image
  seudat@images$slice1@image[,,2] <- gray_image
  seudat@images$slice1@image[,,3] <- gray_image
  
  seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
  metainfo <- seudat@meta.data
  rownames(metainfo) <- paste(metainfo$BCID, rownames(metainfo), sep = ".")
  tmpallinfo <- allinfo[rownames(metainfo), ]
  
  seudat$BLClass <- tmpallinfo$BLClass
  seudat$pseudotime <- tmpallinfo$pseudotime
  seudat$EGFR <- seudat@assays$Spatial$data["EGFR", ]
  seudat$FGFR3 <- seudat@assays$Spatial$data["FGFR3", ]
  seudat$basal <- tmpallinfo$basal
  seudat$luminal <- tmpallinfo$luminal
  
  seudat$pseudotime <- rescale_to_neg1_1(seudat$pseudotime)
  seudat$EGFR <- rescale_to_neg1_1(seudat$EGFR)
  seudat$FGFR3 <- rescale_to_neg1_1(seudat$FGFR3)
  seudat$basal <- rescale_to_neg1_1(seudat$basal)
  seudat$luminal <- rescale_to_neg1_1(seudat$luminal)
  
  newsize = allsams[allsams$Sample == sample, ]$Size
  p1 <- SpatialDimPlot(seudat, group.by = "BLClass", image.alpha=0.2, stroke = NA, pt.size.factor = newsize) + 
    scale_fill_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0"))+
    ggtitle(sam) +
    theme(
      plot.title = element_text(size = 20),
      legend.position = "none"
    )
  
  p2 <- SpatialFeaturePlot(seudat, features = "basal", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = c("#148F28", "#E5E5E5", "#EA71AE"))+
    ggtitle("basal") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p3 <- SpatialFeaturePlot(seudat, features = "pseudotime", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad4)+
    ggtitle("pseudotime") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p4 <- SpatialFeaturePlot(seudat, features = "EGFR", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad3)+
    ggtitle("EGFR") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  p5 <- SpatialFeaturePlot(seudat, features = "FGFR3", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
    scale_fill_gradientn(colors = grad3)+
    ggtitle("FGFR3") +
    theme(
      plot.title = element_text(size = 20)
    )
  
  print(p1|p2|p3|p4|p5)
  
}

dev.off()

#######################spatial distribution in each sample (enlarged)#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

spadat <- readRDS(file.path(subres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data

sam <- "BC05"
sample <- allsams_id[allsams_id$MyID == sam, ]$ID
resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
seudat@images$slice1@image[,,1] <- gray_image
seudat@images$slice1@image[,,2] <- gray_image
seudat@images$slice1@image[,,3] <- gray_image

seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
metainfo <- seudat@meta.data
rownames(metainfo) <- paste(metainfo$BCID, rownames(metainfo), sep = ".")
tmpallinfo <- allinfo[rownames(metainfo), ]

seudat$BLClass <- tmpallinfo$BLClass
seudat$pseudotime <- tmpallinfo$pseudotime
seudat$EGFR <- seudat@assays$Spatial$data["EGFR", ]
seudat$FGFR3 <- seudat@assays$Spatial$data["FGFR3", ]
seudat$basal <- tmpallinfo$basal
seudat$luminal <- tmpallinfo$luminal

seudat$pseudotime <- rescale_to_neg1_1(seudat$pseudotime)
seudat$EGFR <- rescale_to_neg1_1(seudat$EGFR)
seudat$FGFR3 <- rescale_to_neg1_1(seudat$FGFR3)
seudat$basal <- rescale_to_neg1_1(seudat$basal)
seudat$luminal <- rescale_to_neg1_1(seudat$luminal)

coords <- GetTissueCoordinates(seudat)
p <- ggplot(coords, aes(x = x, y = y)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = "Spatial Spot Coordinates",
       x = "x (imagecol)", y = "y (imagerow)")
pdf(file.path(subres, "V2503.figure3.spatial.BC05.examples.test.pdf"), width = 4, height = 4)
print(p)
dev.off()

coords <- GetTissueCoordinates(seudat)
subset_spots <- rownames(coords[coords$x > 8000 & coords$x < 10000 &
                                  coords$y > 15000 & coords$y < 23000, ])


tmpseudat <- subset(seudat, subset = ID %in% subset_spots)
newsize = 14
p1 <- SpatialDimPlot(tmpseudat, group.by = "BLClass", image.alpha=0.2, stroke = NA, pt.size.factor = newsize) + 
  scale_fill_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0"))+
  ggtitle(sam) +
  theme(
    plot.title = element_text(size = 20),
    legend.position = "none"
  )

p4 <- SpatialFeaturePlot(tmpseudat, features = "EGFR", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
  scale_fill_gradientn(colors = c("#E5E5E5", grad2, "#EA71AE"))+
  ggtitle("EGFR") +
  theme(
    plot.title = element_text(size = 20)
  )

p5 <- SpatialFeaturePlot(tmpseudat, features = "FGFR3", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
  scale_fill_gradientn(colors = grad1)+
  ggtitle("FGFR3") +
  theme(
    plot.title = element_text(size = 20)
  )

pdf(file.path(subres, "V2503.figure3.spatial.BC05.examples.enlarged.pdf"), width = 4, height = 9)
print(p1/p4/p5)
dev.off()



sam <- "BC13"
sample <- allsams_id[allsams_id$MyID == sam, ]$ID
resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
seudat@images$slice1@image[,,1] <- gray_image
seudat@images$slice1@image[,,2] <- gray_image
seudat@images$slice1@image[,,3] <- gray_image

seudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
metainfo <- seudat@meta.data
rownames(metainfo) <- paste(metainfo$BCID, rownames(metainfo), sep = ".")
tmpallinfo <- allinfo[rownames(metainfo), ]

seudat$BLClass <- tmpallinfo$BLClass
seudat$pseudotime <- tmpallinfo$pseudotime
seudat$EGFR <- seudat@assays$Spatial$data["EGFR", ]
seudat$FGFR3 <- seudat@assays$Spatial$data["FGFR3", ]
seudat$basal <- tmpallinfo$basal
seudat$luminal <- tmpallinfo$luminal

seudat$pseudotime <- rescale_to_neg1_1(seudat$pseudotime)
seudat$EGFR <- rescale_to_neg1_1(seudat$EGFR)
seudat$FGFR3 <- rescale_to_neg1_1(seudat$FGFR3)
seudat$basal <- rescale_to_neg1_1(seudat$basal)
seudat$luminal <- rescale_to_neg1_1(seudat$luminal)

coords <- GetTissueCoordinates(seudat)
p <- ggplot(coords, aes(x = x, y = y)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = "Spatial Spot Coordinates",
       x = "x (imagecol)", y = "y (imagerow)")
pdf(file.path(subres, "V2503.figure3.spatial.BC13.examples.test.pdf"), width = 4, height = 4)
print(p)
dev.off()

coords <- GetTissueCoordinates(seudat)
subset_spots <- rownames(coords[coords$x > 18000 & coords$x < 23000 &
                                  coords$y > 27000 & coords$y < 32000, ])


tmpseudat <- subset(seudat, subset = ID %in% subset_spots)
newsize = 22
p1 <- SpatialDimPlot(tmpseudat, group.by = "BLClass", image.alpha=0.2, stroke = NA, pt.size.factor = newsize) + 
  scale_fill_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0"))+
  ggtitle(sam) +
  theme(
    plot.title = element_text(size = 20),
    legend.position = "none"
  )

tmpseudat$EGFR[tmpseudat$EGFR > 0.7] <- 0.7
tmpseudat$EGFR[tmpseudat$EGFR < 0.1] <- 0.1
p4 <- SpatialFeaturePlot(tmpseudat, features = "EGFR", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
  scale_fill_gradientn(colors = grad2)+
  ggtitle("EGFR") +
  theme(
    plot.title = element_text(size = 20)
  )

tmpseudat$FGFR3[tmpseudat$FGFR3 > 0.0] <- 0.0
tmpseudat$FGFR3[tmpseudat$FGFR3 < -1] <- -1
p5 <- SpatialFeaturePlot(tmpseudat, features = "FGFR3", image.alpha = 0.2, stroke = NA, pt.size.factor = newsize) +
  scale_fill_gradientn(colors = grad1)+
  ggtitle("FGFR3") +
  theme(
    plot.title = element_text(size = 20)
  )

pdf(file.path(subres, "V2503.figure3.spatial.BC13.examples.enlarged.pdf"), width = 4, height = 9)
print(p1/p4/p5)
dev.off()


