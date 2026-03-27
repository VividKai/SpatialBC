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

#######################clinical information#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
rownames(allsams_id) <- allsams_id$ID
allsams_id <- allsams_id[order(allsams_id$Responder, decreasing = T), ]
allsams_id$MyID <- factor(allsams_id$MyID, allsams_id$MyID)

nature_palette <- c(
  "MVAC" = "#D9CC77", "GC" = "#B3A2A5", "GTA" = "#A0BCD9",  # Softer Treatment colors
  "P" = "#377EB8", "N" = "#FF7F00", "NS" = "gray70",        # MTAP colors (unchanged P/N)
  "T2" = "#F9EAE1", "T3" = "#D5A6BD", "T4" = "#C85A7A",     # Softer Stage colors
  "Yes" = "#4DAF4A", "No" = "#E41A1C",                      # Responder colors (unchanged)
  "PT" = "#6B97C7", "IT" = "#E99D50"                        # Adjusted colors for better contrast
)

allsamsdis <- data.frame()
for (i in 1:length(allsams$Sample)) {
  sample = allsams$Sample[i]
  print(i)
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  rctddat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  nspots <- nrow(rctddat@meta.data)
  if ("Tumor Cell" %in% rctddat$Level2_Annot) {
    rctddat <- subset(rctddat, subset = Level2_Annot == "Tumor Cell")
    DefaultAssay(rctddat) <- "Spatial"
    rctddat$MTAP <- as.numeric(GetAssayData(rctddat, slot = "data")["MTAP", ])
    rctddat$CDKN2A <- as.numeric(GetAssayData(rctddat, slot = "data")["CDKN2A", ])
    tmpres <- data.frame(ID = unique(rctddat$BCID),
                         nSpots = nspots,
                         TumorSpots = nrow(rctddat@meta.data))
  }else{
    tmpres <- data.frame(ID = unique(rctddat$BCID),
                         nSpots = nspots,
                         TumorSpots = 0)
  }
  
  allsamsdis <- rbind(allsamsdis, tmpres)
  
}

allsams_id <- merge(allsams_id, allsamsdis, by.x = "MyID", by.y = "ID")
allsams_id <- allsams_id[order(allsams_id$MyID), ]


# Prepare data for heatmap (categorical variables)
allsams_melt <- allsams_id %>%
  select(MyID, Responder, Treatment, Stage, TLS) %>%
  melt(id.vars = "MyID", variable.name = "Annotation", value.name = "Value")

# Heatmap (Top panel)
allsams_melt$Type <- ifelse(allsams_melt$Annotation %in% c("Responder", "Treatment", "Stage"), "tile", "circle")

# Heatmap
heatmap_plot <- ggplot(allsams_melt, aes(x = MyID, y = Annotation, fill = Value)) +
  geom_tile(data = subset(allsams_melt, Type == "tile"), color = "black") +
  geom_point(data = subset(allsams_melt, Type == "circle"), shape = 21, size = 6, color = "black") +
  scale_fill_manual(values = nature_palette) +
  labs(title = "Annotation Heatmap", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )





# Prepare data for bar plots (numeric variables)
allsams_melt_numeric <- allsams_id %>%
  select(MyID, nSpots, TumorSpots) %>%
  melt(id.vars = "MyID", variable.name = "Variable", value.name = "Value")

# Bar plot (Bottom panel)
bar_plot <- ggplot(allsams_melt_numeric, aes(x = MyID, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("nSpots" = "#C9A9A6", "TumorSpots" = "#92A8D1")) +
  labs(x = "Sample ID", y = "Expression / Tumor Spots") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

# Combine plots using patchwork
final_plot <- heatmap_plot / bar_plot + plot_layout(heights = c(1, 2))  # Heatmap smaller than bar plot

pdf(file.path(outdir, "V2503.figure1.sample.info.pdf"), width = 8, height = 6)
print(final_plot)
dev.off()


#######################clinical annotation distribution#######################
cliant <- fread("../color.map", header = T, stringsAsFactors = F, data.table = F)
cliant <- cliant[cliant$Level1_Annot != "Unlabeled", ]

desired_order <- c('Tumor Region', 'Normal Region', 'Other Region')
cliant$Level1_Annot <- factor(cliant$Level1_Annot, levels = desired_order)
order2 <- c('Tumor Cell', 'Tumor Stromas', 'Immune Cells', 'Necrosis', 'TLS/TIL/LA', 
            'Vessels', 'Urotheliums', 'Normal Stromas', 'Muscles', 'Others')
cliant$Level2_Annot <- factor(cliant$Level2_Annot, levels = order2)

cliant <- cliant[order(cliant$Level1_Annot, cliant$Level2_Annot, cliant$Official_Annot, decreasing = T), ]

cliant$Official_Color <- cliant$Level2_Color

unique_annotations <- unique(cliant[, c("Level1_Annot", "Level1_Color")])
color_mapping1 <- setNames(unique_annotations$Level1_Color, unique_annotations$Level1_Annot)
unique_annotations <- unique(cliant[, c("Level2_Annot", "Level2_Color")])
color_mapping2 <- setNames(unique_annotations$Level2_Color, unique_annotations$Level2_Annot)
unique_annotations <- unique(cliant[, c("Official_Annot", "Official_Color")])
color_mapping3 <- setNames(unique_annotations$Official_Color, unique_annotations$Official_Annot)
color_mapping <- c(color_mapping1, color_mapping2, color_mapping3)

# Step 1
df <- cliant %>%
  make_long(Level1_Annot, Level2_Annot, Official_Annot)
df

df$node <- factor(df$node, c(as.character(unique(cliant$Level1_Annot)), as.character(unique(cliant$Level2_Annot)), unique(cliant$Official_Annot)))
df$next_node <- factor(df$next_node, c(as.character(unique(cliant$Level1_Annot)), as.character(unique(cliant$Level2_Annot)), unique(cliant$Official_Annot)))

# Chart 1
pl <- ggplot(df, aes(x = x
                     , next_x = next_x
                     , node = node
                     , next_node = next_node
                     , fill = factor(node)
                     , label = node)
)
pl <- pl +geom_sankey(flow.alpha = 0.5
                      , node.color = "black"
                      ,show.legend = FALSE)
pl <- pl +geom_sankey_text(size = 3, color = "black", hjust = -0.5)
pl <- pl +  theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(fill = 'Nodes')
pl <- pl + scale_fill_manual(values = color_mapping )

ggsave(file.path(outdir, "V2503.figure1.clinic.annot.sankey.pdf"), pl, width = 6, height = 5)


alldats <- data.frame()
for (sams in 1:length(allsams$Sample)) {
  ###deal with spatial spot data
  sample = allsams$Sample[sams]
  #sample = "S-04-021004"
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  
  spadat <- subset(seudat, subset = Level1_Annot %in% c("Tumor Region", "Normal Region", "Other Region"))
  spadat$Sample <- sample
  
  ###run nmf
  #spadat <- runNMF(spadat, assay="SCT", slot="scale.data", k=5)
  alldats <- rbind(alldats, spadat@meta.data)
}

write.table(alldats, file.path(outdir, "V2503.figure1.metainfo.summary.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

finalres <- fread(file.path(outdir, "V2503.figure1.metainfo.summary.txt"), header = T, stringsAsFactors = F, data.table = F)

cellratio <- prop.table(table(finalres$Level2_Annot, as.character(finalres$Sample)), margin = 1)
cellratio <- as.data.frame(cellratio)

order2 <- rev(c('Tumor Cell', 'Tumor Stromas', 'Immune Cells', 'Necrosis', 'TLS/TIL/LA', 
                'Vessels', 'Urotheliums', 'Normal Stromas', 'Muscles', 'Others'))
colors <- rev(c('#F59696', '#5f88c9', '#FB40B0', '#F7A7D6', '#F7C29B', 
                '#FFEBAE', '#5fc9ad', '#5fadc9', '#FFF400', '#CCCCCC'))

cellratio$Var1 <- factor(cellratio$Var1, order2)

p <- ggplot(cellratio) +
  geom_bar(aes(x =Var1, y= Freq, fill = Var2), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  #coord_flip()+
  theme_classic()+
  scale_y_continuous(expand  = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
p

pdf(file.path(outdir, "V2503.figure1.metainfo.level22sample.summary.bar.pdf"), width = 6, height = 3)
print(p)
dev.off()


cellnum <- table(finalres$Level2_Annot)
cellnum <- as.data.frame(cellnum)
cellnum$Var1 <- factor(cellnum$Var1, order2)

# Ensure that 'colors' is defined correctly before this
p <- ggplot(cellnum) +
  geom_col(aes(x = Var1, y = Freq, fill = Var1), width = 0.7) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Sample', y = 'Ratio') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_flip()

# Output to PDF
pdf(file.path(outdir, "V2503.figure1.metainfo.level22sample.spotsnumber.pdf"), width = 6, height = 3)
print(p)
dev.off()

#######################show pathologist annotation examples#######################
n_rows <- 2
n_cols <- 11
total_plots <- n_rows * n_cols
plot_list <- vector("list", total_plots)
for (i in 1:length(allsams$Sample)) {
  sample = allsams$Sample[i]
  print(i)
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  spadat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  gray_image <- apply(spadat@images$slice1@image[,,1:3], c(1,2), mean)
  spadat@images$slice1@image[,,1] <- gray_image
  spadat@images$slice1@image[,,2] <- gray_image
  spadat@images$slice1@image[,,3] <- gray_image
  
  spadat$Level2_Annot <- factor(spadat$Level2_Annot, unique(spadat$Level2_Annot))
  color_v=unique(spadat@meta.data[order(spadat$Level2_Annot), ]$Level2_Color)
  
  p <- SpatialDimPlot(spadat, group.by = "Level2_Annot", image.alpha=0.6, stroke = NA, pt.size.factor = allsams[allsams$Sample == sample, ]$Size) + 
    scale_fill_manual(values = color_v)+
    ggtitle(sample) +
    theme(
      plot.title = element_text(size = 40),
      legend.position = "none"
    )
  plot_list[[i]] <- p 
}

colormap <- fread("../color.map", header = T, stringsAsFactors = F, data.table = F)
colormap <- colormap[colormap$Level2_Annot != "Unlabel", ]
colormap$Level2_Annot <- factor(colormap$Level2_Annot, unique(colormap$Level2_Annot))
color_v=unique(colormap[order(colormap$Level2_Annot), ]$Level2_Color)

p2 <- ggplot(colormap, aes(x = Level2_Annot, y = 1, color = Level2_Annot)) +
  geom_point(size = 5) +  # Adjust the size of the points
  scale_color_manual(values = color_v) +
  theme_void() +  # Remove axes, gridlines, and background
  theme(
    legend.position = "right",  # Adjust legend position
    legend.title = element_blank()  # Remove the legend title
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, "V2503.figure1.epi.clinical.info.pdf"), width = 48, height = 12)
gridExtra::grid.arrange(grobs = plot_list, nrow = n_rows, ncol = n_cols)
print(p2)
dev.off()


#######################plot RCTD proportions for each patients#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
allcellfrac <- data.frame()
for (i in 1:length(allsams$Sample)) {
  sample = allsams$Sample[i]
  print(i)
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  spadat1 <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  spadat1$Level2_Annot <- as.character(spadat1$Level2_Annot)
  
  if ("Tumor Cell" %in% unique(spadat1$Level2_Annot)) {
    tmpspadat <- subset(spadat1, subset = Level2_Annot %in% c("Tumor Cell"))
    Idents(tmpspadat) <- tmpspadat$Level2_Annot
    
    rctdrds <- paste("../9.spaAnaly/NewCombined.Analy/01.RCTD/", sample, sep = "")
    rctddat <- fread(file.path(rctdrds, "cell.fraction.full.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F)
    rownames(rctddat) <- rctddat$V1
    rctddat <- rctddat[, -which(names(rctddat) %in% c("V1"))]
    rctddat <- rctddat / rowSums(rctddat)
    
    ########rctd
    tmpspadat <- AddMetaData(tmpspadat, rctddat)
    
    metainfo <- tmpspadat@meta.data[, c(names(rctddat))]
    tmpres <- data.frame(Sample = sample,
                         Cells = names(metainfo),
                         Mean_frac = as.numeric(colMeans(metainfo)) / sum(as.numeric(colMeans(metainfo))))
    
    allcellfrac <- rbind(allcellfrac, tmpres)
    
  }
  
}

nature_palette <- c("#F4D7A3", "#7EBCD4", "#89D2B3", "#F7E6A8", "#6E9ECF", "#E6A987", "#D9A5C9", "#BBBBBB", 
                    "#D1C197", "#50A6A5", "#E29D68", "#B9A3E3", "#8CBF75", "#C78DA1", "#F59696")
allcellfrac$Cells <- factor(allcellfrac$Cells, c("T", "NK", "B", "Plasma", "LAMP3.DC", "cDC1", "cDC2", 
                                                 "pDC", "Macro", "Mast", "Mono", 
                                                 "Endo", "mCAF", "iCAF", "Epithelial"))
allcellfrac <- merge(allcellfrac, allsams_id, by.x = "Sample", by.y = "ID")

ordered_samples <- allcellfrac %>%
  filter(Cells == "Epithelial") %>%
  arrange(desc(Mean_frac)) %>%
  pull(MyID)

# Convert Sample to factor for ordering
allcellfrac$MyID <- factor(allcellfrac$MyID, levels = ordered_samples)

# Plot stacked bar plot
p <- ggplot(allcellfrac, aes(x = MyID, y = Mean_frac, fill = Cells)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = nature_palette) +  # Use the predefined nature palette
  labs(x = "MyID", y = "Mean Fraction", fill = "Cell Type") +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic2() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "right"
  )

pdf(file.path(outdir, "V2503.figure1.rctd.cellfrac.sample.pdf"), width = 5, height = 2.5)
print(p)
dev.off()


#######################RCTD proportions correlated with signature score#######################
siggenes <- fread("../../11.JJ_Landscape_BC/4.Results/combined.tme.epi.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)

allres <- data.frame()
for (j in 1:length(allsams$Sample)) {
  print(j)
  sample = allsams$Sample[j]
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  spadat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  rctdrds <- paste("../9.spaAnaly/NewCombined.Analy/01.RCTD/", sample, sep = "")
  rctddat <- fread(file.path(rctdrds, "cell.fraction.full.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F)
  rownames(rctddat) <- rctddat$V1
  
  tmpspadat1 <- spadat
  siggenes <- fread("../../11.JJ_Landscape_BC/4.Results/combined.tme.epi.markers.top50.txt", header = T, stringsAsFactors = F, data.table = F)
  tmppath_lst <- split(siggenes$gene, siggenes$cluster)
  tmpspadat1 <- AddModuleScore(tmpspadat1,
                               features = tmppath_lst,
                               name="TCELLSIGNATURE",
                               nbin = 10)
  colnames(tmpspadat1@meta.data)[grep("TCELLSIGNATURE", colnames(tmpspadat1@meta.data))] <- names(tmppath_lst)
  
  ########rctd
  tmpspadat <- spadat
  tmpspadat <- AddMetaData(tmpspadat, rctddat)
  for (cells in gsub("\\+", ".", unique(siggenes$cluster))) {
    if (cells %in% names(rctddat)) {
      res <- data.frame(Sig = as.numeric(tmpspadat1@meta.data[, gsub("\\.", "\\+", cells)]),
                        Deconv = as.numeric(tmpspadat@meta.data[, cells]))
      corre <- cor.test(res$Deconv,res$Sig,method="pearson")
      
      tmpres <- data.frame(ID1 = cells,
                           Cor = corre$estimate[[1]],
                           Sample = sample)
      allres <- rbind(allres, tmpres)
    }
  }
}


normalized_data <- allres %>%
  group_by(ID1) %>%
  summarise(Mean_Cor = median(Cor, na.rm = TRUE)) %>%
  ungroup()


nature_palette <- c("#F4D7A3", "#7EBCD4", "#89D2B3", "#F7E6A8", "#6E9ECF", "#E6A987", "#D9A5C9", "#BBBBBB", 
                    "#D1C197", "#50A6A5", "#E29D68", "#B9A3E3", "#8CBF75", "#C78DA1", "#F59696")
normalized_data$ID1 <- factor(normalized_data$ID1, c("T", "NK", "B", "Plasma", "LAMP3.DC", "cDC1", "cDC2", 
                                                     "pDC", "Macro", "Mast", "Mono", 
                                                     "Endo", "mCAF", "iCAF", "Epithelial"))

p <- ggplot(normalized_data, aes(x = ID1, y = Mean_Cor, fill = ID1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = nature_palette) +  # Use the predefined nature palette
  geom_hline(yintercept = mean(normalized_data$Mean_Cor), linetype = "dashed", color = "black", size = 0.8) +
  labs(x = "MyID", y = "Mean Fraction", fill = "Cell Type") +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic2() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = "right"
  )

pdf(file.path(outdir, "V2503.figure1.rctd.cellfrac.corr.signature.pdf"), width = 5, height = 3)
print(p)
dev.off()

#######################compare some things#######################
tlsgenes <- fread("/rsrch5/scratch/genomic_med/kyu3/softwares/istartls/data/TLS_signature_david2/tls/gene-names.txt",
                  header = F, stringsAsFactors = F, data.table = F)
tlssignature <- data.frame(Genes = tlsgenes$V1,
                           Signature = "TLS")
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
names(hall) <- c("Signature", "Genes")
hall$Signature <- gsub("HALLMARK_", "", hall$Signature)

hall_tls <- rbind(hall, tlssignature)


allscores <- split(hall_tls$Genes, hall_tls$Signature)

allcellfrac <- data.frame()
for (j in 1:length(allsams$Sample)) {
  print(j)
  sample = allsams$Sample[j]
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  spadat1 <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  spadat1$Level2_Annot <- as.character(spadat1$Level2_Annot)
  DefaultAssay(spadat1) <- "Spatial"
  spadat1 <- AddModuleScore(spadat1,
                            features = allscores,
                            nbin = 12,
                            name="TCELLSIGNATURE")
  colnames(spadat1@meta.data)[grep("TCELLSIGNATURE", colnames(spadat1@meta.data))] <- names(allscores)
  spadat1$EPCAM <- FetchData(spadat1, vars = "EPCAM")[, 1]
  spadat1$UPK3A <- FetchData(spadat1, vars = "UPK3A")[, 1]
  spadat1$GATA3 <- FetchData(spadat1, vars = "GATA3")[, 1]
  spadat1$FOXA1 <- FetchData(spadat1, vars = "FOXA1")[, 1]
  
  tmpspadat <- spadat1
  Idents(tmpspadat) <- tmpspadat$Level2_Annot
  
  rctdrds <- paste("../9.spaAnaly/NewCombined.Analy/01.RCTD/", sample, sep = "")
  rctddat <- fread(file.path(rctdrds, "cell.fraction.full.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F)
  rownames(rctddat) <- rctddat$V1
  rctddat <- rctddat[, -which(names(rctddat) %in% c("V1"))]
  #rctddat <- rctddat[, -which(names(rctddat) %in% c("V1", "Epithelial"))]
  rctddat <- rctddat / rowSums(rctddat)
  
  ########rctd
  tmpspadat <- AddMetaData(tmpspadat, rctddat)
  
  metainfo <- tmpspadat@meta.data[, c("Level2_Annot", "EPCAM", "UPK3A", "GATA3", "FOXA1", names(rctddat), names(allscores))]
  
  allcellfrac <- rbind(allcellfrac, metainfo)
  
  
}


allcellfrac <- allcellfrac[allcellfrac$Level2_Annot != "Others", ]


tmpfrac <- allcellfrac
tmpfrac$Class <- ifelse(tmpfrac$Level2_Annot == "Tumor Cell", "Tumor Cell", "Non-T regions")
tmpfrac$Class <- factor(tmpfrac$Class, c("Tumor Cell", "Non-T regions"))
colorv = c("#F59696", "#94FFD8")
my_comp <-  list(c("Tumor Cell", "Non-T regions"))
p1 <- ggplot(tmpfrac, aes(x = Class, y = Epithelial, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Epithelial proportion") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)

p2 <- ggplot(tmpfrac, aes(x = Class, y = GATA3, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "normalized GATA3 RNA") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")


tmpfrac <- allcellfrac
tmpfrac$Class <- ifelse(tmpfrac$Level2_Annot == "TLS/TIL/LA", "TLS", "Non-TLS regions")
tmpfrac$Class <- factor(tmpfrac$Class, c("TLS", "Non-TLS regions"))
colorv = c("#F59696", "#94FFD8")
my_comp <-  list(c("TLS", "Non-TLS regions"))
p3 <- ggplot(tmpfrac, aes(x = Class, y = TLS, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "TLS signature") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)

p4 <- ggplot(tmpfrac, aes(x = Class, y = B, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B proportion") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

p5 <- ggplot(tmpfrac, aes(x = Class, y = T, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "T proportion") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")


tmpfrac <- allcellfrac
tmpfrac$Class <- ifelse(tmpfrac$Level2_Annot %in% c("Vessels"), "Vessels", "Non-Vessels")
tmpfrac$Class <- factor(tmpfrac$Class, c("Vessels", "Non-Vessels"))
colorv = c("#F59696", "#94FFD8")
my_comp <-  list(c("Vessels", "Non-Vessels"))
p6 <- ggplot(tmpfrac, aes(x = Class, y = Endo, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Endo proportion") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)

tmpfrac$CAF <- tmpfrac$iCAF + tmpfrac$mCAF
tmpfrac$Class <- ifelse(tmpfrac$Level2_Annot %in% c("Muscles"), "Muscles", "Non-Muscles")
tmpfrac$Class <- factor(tmpfrac$Class, c("Muscles", "Non-Muscles"))
colorv = c("#F59696", "#94FFD8")
my_comp <-  list(c("Muscles", "Non-Muscles"))
p7 <- ggplot(tmpfrac, aes(x = Class, y = CAF, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Fibro proportion") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")


tmpfrac <- allcellfrac
tmpfrac$Class <- ifelse(tmpfrac$Level2_Annot %in% c("Necrosis"), "Necrosis", "Non-Necrosis")
tmpfrac$Class <- factor(tmpfrac$Class, c("Necrosis", "Non-Necrosis"))
colorv = c("#F59696", "#94FFD8")
my_comp <-  list(c("Necrosis", "Non-Necrosis"))
p8 <- ggplot(tmpfrac, aes(x = Class, y = HYPOXIA, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "HYPOXIA Signature") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)


p9 <- ggplot(tmpfrac, aes(x = Class, y = APOPTOSIS, fill = Class)) +
  geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.1, outlier.size = 0, position = position_dodge(0.9)) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "APOPTOSIS Signature") +
  stat_compare_means(comparisons = my_comp)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)

pdf(file.path(outdir, "V2503.figure1.pathologist.annot.compare.pdf"), width = 20, height = 4)
print(p1|p2|p3|p4|p5|p6|p7|p8|p9)
dev.off()












#######################pathologist annotation markers#######################
colormap <- fread("../color.map", header = T, stringsAsFactors = F, data.table = F)
colormap <- colormap[colormap$Level2_Annot != "Unlabel", ]
colormap$Level2_Annot <- as.character(colormap$Level2_Annot)

combined_seurat <- readRDS(file.path(resdir, "allsam.merged.seurat.umap.rds"))
combined_seurat <- subset(combined_seurat, subset = Level2_Annot != "Others")

spadat$Level2_Annot <- factor(spadat$Level2_Annot, unique(colormap$Level2_Annot))

Idents(combined_seurat) <- combined_seurat$Level2_Annot
markers <- FindAllMarkers(combined_seurat, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.1)

top50 <- top_n(group_by(markers,cluster), n = 6, wt = avg_log2FC)
write.table(markers, file.path(outdir, "V2503.figure1.pathologist.annot.markers.txt"), row.names = F, col.names = T, sep = "\t", quote = F)

top50$cluster <- factor(top50$cluster, unique(colormap$Level2_Annot))
top50 <- top50[order(top50$cluster, decreasing = T), ]
genes <- top50$gene

p <- DotPlot(combined_seurat, features = unique(genes))
tmpres <- p$data

tmpres$features.plot <- factor(tmpres$features.plot, unique(genes))
tmpres$id <- factor(tmpres$id, unique(colormap$Level2_Annot))

plot <- ggplot(tmpres, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
  geom_tile() +
  scale_fill_gradient2(low = "lightblue", mid = "white", high = "#750E21", midpoint = 0.0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Feature", y = "Cell Type", fill = "Scaled Exp")


pdf(file.path(outdir, "V2503.figure1.pathologist.annot.markers.pdf"), width = 18, height = 3)
print(plot)
dev.off()



