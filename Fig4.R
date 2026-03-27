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

dir.create("../9.spaAnaly/NewCombined.Analy/05.distance")
subres <- "../9.spaAnaly/NewCombined.Analy/05.distance"

#######################different region class distribution#######################
tmpsubres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data


combined_info <- data.frame()
for (i in 1:nrow(allsams)) {
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  print(i)
  sample <- allsams$Sample[i]
  if (file.exists(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))) {
    resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
    seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))
    
    tmpres <- seudat@meta.data
    rownames(tmpres) <- paste(tmpres$BCID, rownames(tmpres), sep = ".")
    
    combined_info <- rbind(combined_info, tmpres)
    
  }
}

finalres <- merge(allinfo, combined_info[, c("CellID", "Regions")], by.x = "CellID", by.y = "CellID")


chisq_test <- chisq.test(table(finalres$BLClass, finalres$Regions))
finalres$Regions <- factor(finalres$Regions, c("Tumor_Core", "Tumor_Inter", "Tumor_Bound"))
finalres$BLClass <- factor(finalres$BLClass, c("Lum", "IM", "B.EMT", "B.IS"))
tmpres <- as.data.frame(prop.table(table(finalres$Regions, finalres$BLClass), margin = 1))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0"))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()


pdf(file.path(subres, "V2503.figure4.spatial.region.subtype.pdf"), width = 7, height = 6)
print(p)
dev.off()

#######################spatial coherence#######################
tmpsubres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"
library(patchwork)
library(parallel)
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data
usefulsams <- unique(spadat$BCID)

#!!! Run nearest spots for each spots
for (sam in usefulsams) {
  print(sam)
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  
  metainfo <- seudat@meta.data
  rownames(metainfo) <- paste(metainfo$BCID, rownames(metainfo), sep = ".")
  tmpallinfo <- allinfo[rownames(metainfo), ]
  
  seudat$AssignedCombo <- tmpallinfo$BLClass
  seudat$AssignedCombo[is.na(seudat$AssignedCombo)] <- "NonTRegion"
  
  bar2mp <- seudat@meta.data[, c("ID", "AssignedCombo")]
  posfile <- paste("../1.spaceranger/", sample, "/outs/spatial/tissue_positions.csv", sep = "")
  posdat <- fread(posfile, header = T, sep = ",", stringsAsFactors = F, data.table = F)
  rownames(posdat) <- posdat$barcode
  posdat <- posdat[bar2mp$ID, ]
  
  coords <- posdat[, c("barcode", "array_row", "array_col")]
  names(coords) <- c("barcodes", "x", "y")
  all_spots <- ppp(coords$x, coords$y, window=owin(range(coords$x), range(coords$y)))
  allres <- data.frame()
  
  for (selbarcodes in coords$barcodes) {
    region_spots <- coords[coords$barcodes == selbarcodes, ]
    # Convert coordinates to spatial points
    region_points <- ppp(region_spots$x, region_spots$y, window=owin(range(coords$x), range(coords$y)))
    # Calculate distances from each spot to the closest spot in the region
    distances <- nncross(all_spots, region_points, what="dist")
    newdis <- distances / sort(distances)[2]
    closet <- coords$barcodes[newdis < 1.5]
    if (length(closet) > 4) {
      selmps <- bar2mp[selbarcodes, ]
      closetmps <- bar2mp[closet, ]
      tmpres <- data.frame(ID = selbarcodes,
                           Coherence = sum(closetmps$AssignedCombo %in% selmps$AssignedCombo) / nrow(closetmps),
                           Distance = 1,
                           MPs = selmps$AssignedCombo)
      allres <- rbind(allres, tmpres)
    }
    newdis <- distances / sort(distances)[2]
    closet <- coords$barcodes[newdis < 2.5]
    if (length(closet) > 12) {
      selmps <- bar2mp[selbarcodes, ]
      closetmps <- bar2mp[closet, ]
      tmpres <- data.frame(ID = selbarcodes,
                           Coherence = sum(closetmps$AssignedCombo %in% selmps$AssignedCombo) / nrow(closetmps),
                           Distance = 2,
                           MPs = selmps$AssignedCombo)
      allres <- rbind(allres, tmpres)
    }
    newdis <- distances / sort(distances)[2]
    closet <- coords$barcodes[newdis < 3.5]
    if (length(closet) > 24) {
      selmps <- bar2mp[selbarcodes, ]
      closetmps <- bar2mp[closet, ]
      tmpres <- data.frame(ID = selbarcodes,
                           Coherence = sum(closetmps$AssignedCombo %in% selmps$AssignedCombo) / nrow(closetmps),
                           Distance = 3,
                           MPs = selmps$AssignedCombo)
      allres <- rbind(allres, tmpres)
    }
    
  }
  
  saveRDS(allres, file.path(resdir, paste(sample, ".seurat.umap.spatial.coherence.rds", sep = "")))
}


allres <- data.frame()
for (i in 1:length(usefulsams)) {
  print(i)
  sam = usefulsams[i]
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  tmpres <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.coherence.rds", sep = "")))
  tmpres$Sample <- sample
  tmpres <- tmpres[!is.nan(tmpres$Coherence), ]
  allres <- rbind(allres, tmpres)
  
  
}

allres <- merge(allres, allsams_id, by.x = "Sample", by.y = "ID")
#allres <- allres[allres$Coherence > 0, ]
allres <- allres[allres$Distance %in% c(1,2,3), ]
allres <- allres[allres$MPs != "NonTRegion", ]


result_tf <- allres %>%
  group_by(Sample, MPs) %>%
  summarise(mean_Coherence = mean(Coherence, na.rm = TRUE), .groups = "drop")

result_tf <- merge(result_tf, allsams_id, by.x = "Sample", by.y = "ID")
result_tf$MPs <- factor(result_tf$MPs, c("Lum", "IM", "B.EMT", "B.IS"))

sample_order <- result_tf %>%
  group_by(MyID) %>%
  summarise(mean_mean_Coherence = mean(mean_Coherence)) %>%
  arrange(desc(mean_mean_Coherence)) %>% # Sort in descending order
  pull(MyID)

# Update the factor levels for Sample based on the sorted order
result_tf$MyID <- factor(result_tf$MyID, levels = sample_order)

summary_df <- result_tf %>%
  group_by(MyID) %>%
  summarise(
    mean_value = mean(mean_Coherence, na.rm = TRUE),
    se_value = sd(mean_Coherence, na.rm = TRUE) / sqrt(n())
  )

p <- ggplot(result_tf, aes(x = MyID, y = mean_Coherence)) +
  geom_point(aes(color = MPs), size = 2) + # Add jitter for better separation
  geom_errorbar(data = summary_df, aes(x = MyID, y = mean_value, 
                                       ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.2)+
  scale_color_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0")) + # Add appropriate colors for MPs
  labs(x = "MyID", y = "Mean Spatial Coherence", color = "MetaProgram") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) # Rotate x-axis labels for readability

class_annotation <- allsams_id
class_annotation <- class_annotation[class_annotation$MyID %in% sample_order, ]
class_annotation$MyID <- factor(class_annotation$MyID, levels = sample_order)

p2 <- ggplot(class_annotation, aes(x = MyID, y = 1, fill = Responder)) +
  geom_tile(width = 1, height = 1) +  # Add annotation bar for each sample
  scale_fill_manual(values = c("Yes" = "#4DAF4A", "No" = "#E41A1C")) +  # Set colors for Class annotation
  theme_void() +  # Remove axes and background for cleaner annotation bar
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank())  # Remove axis labels


combined_plot <- p / p2 + plot_layout(heights = c(8, 1)) # Adjust heights (4:1 means p2 is 0.2 of the height)

# Save to PDF
pdf(file.path(subres, "V2503.figure4.spatial.coherence.pdf"), width = 6, height = 4)
print(combined_plot)
dev.off()


#######################immune cell colocalization (calculate score)#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data

usefulsams <- unique(spadat$BCID)
allres <- data.frame()
for (sam in usefulsams) {
  print(sam)
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
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
  
  spadat1 <- seudat
  DefaultAssay(spadat1) <- "Spatial"
  
  spadat1 <- subset(spadat1, subset = Level2_Annot %in% c("Tumor Cell"))
  rctdrds <- paste("../9.spaAnaly/", sample, "/15.RCTD/", sep = "")
  rctddat <- fread(file.path(rctdrds, "cell.fraction.full.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F)
  rownames(rctddat) <- rctddat$V1
  rctddat <- rctddat[, -which(names(rctddat) %in% c("V1", "Epithelial"))]
  
  spadat1 <- AddMetaData(spadat1, rctddat)
  
  # Get the expression values for TROP2
  trop2_positive_spots <- spadat1@meta.data[spadat1$BLClass %in% "Lum", ]$ID
  
  for (cells in names(rctddat)) {
    dc_expression <- FetchData(spadat1, vars = cells)
    dc_positive_spots <- dc_expression > 0.1
    dc_positive_spots <- rownames(dc_expression)[dc_positive_spots]
    
    observed_overlap <- length(intersect(trop2_positive_spots, dc_positive_spots))
    print(paste("Observed overlap:", observed_overlap))
    n_permutations <- 10000
    null_overlaps <- numeric(n_permutations)
    all_spots <- colnames(spadat1)  # Assuming spadat1 is your Seurat object
    set.seed(123)  # For reproducibility
    for (i in 1:n_permutations) {
      # Assign the first n_trop2 spots to the "trop2" group and the first n_cd274 spots to the "cd274" group
      n_trop2 <- length(trop2_positive_spots)
      n_cd274 <- length(dc_positive_spots)
      
      shuffled_trop2 <- sample(all_spots, size = n_trop2, replace = F)
      shuffled_cd274 <- sample(all_spots, size = n_cd274, replace = F)
      
      # Calculate the overlap between the shuffled groups
      null_overlaps[i] <- length(intersect(shuffled_trop2, shuffled_cd274))
    }
    null_mean <- mean(null_overlaps)
    null_sd <- sd(null_overlaps)
    # Calculate Z-score
    z_score <- (observed_overlap - null_mean) / null_sd
    # Calculate two-tailed p-value using normal distribution
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    
    print(paste(cells, ": P-value:", p_value))
    
    tmpres <- data.frame(ID1 = "Luminal",
                         ID2 = cells,
                         ES = log2(observed_overlap+1) - log2(mean(null_overlaps) + 1),
                         Pvalue = p_value,
                         Sample = sample)
    allres <- rbind(allres, tmpres) 
  }
}

write.table(allres, file.path(subres, "Part1.2.spatial.TROP2.colocalization.immune.rctd.luminal.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)

allres <- data.frame()
for (sam in usefulsams) {
  print(sam)
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
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
  
  spadat1 <- seudat
  DefaultAssay(spadat1) <- "Spatial"
  
  spadat1 <- subset(spadat1, subset = Level2_Annot %in% c("Tumor Cell"))
  rctdrds <- paste("../9.spaAnaly/", sample, "/15.RCTD/", sep = "")
  rctddat <- fread(file.path(rctdrds, "cell.fraction.full.txt"), sep = "\t", header = T, stringsAsFactors = F, data.table = F)
  rownames(rctddat) <- rctddat$V1
  rctddat <- rctddat[, -which(names(rctddat) %in% c("V1", "Epithelial"))]
  #rctddat <- rctddat / rowSums(rctddat)
  
  spadat1 <- AddMetaData(spadat1, rctddat)
  
  # Get the expression values for TROP2
  trop2_positive_spots <- spadat1@meta.data[spadat1$BLClass %in% c("B.EMT", "B.IS"), ]$ID
  
  for (cells in names(rctddat)) {
    dc_expression <- FetchData(spadat1, vars = cells)
    dc_positive_spots <- dc_expression > 0.1
    dc_positive_spots <- rownames(dc_expression)[dc_positive_spots]
    
    observed_overlap <- length(intersect(trop2_positive_spots, dc_positive_spots))
    print(paste("Observed overlap:", observed_overlap))
    n_permutations <- 10000
    null_overlaps <- numeric(n_permutations)
    all_spots <- colnames(spadat1)  # Assuming spadat1 is your Seurat object
    set.seed(123)  # For reproducibility
    for (i in 1:n_permutations) {
      # Assign the first n_trop2 spots to the "trop2" group and the first n_cd274 spots to the "cd274" group
      n_trop2 <- length(trop2_positive_spots)
      n_cd274 <- length(dc_positive_spots)
      
      shuffled_trop2 <- sample(all_spots, size = n_trop2, replace = F)
      shuffled_cd274 <- sample(all_spots, size = n_cd274, replace = F)
      
      # Calculate the overlap between the shuffled groups
      null_overlaps[i] <- length(intersect(shuffled_trop2, shuffled_cd274))
    }
    null_mean <- mean(null_overlaps)
    null_sd <- sd(null_overlaps)
    # Calculate Z-score
    z_score <- (observed_overlap - null_mean) / null_sd
    # Calculate two-tailed p-value using normal distribution
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    
    print(paste(cells, ": P-value:", p_value))
    
    tmpres <- data.frame(ID1 = "Basal",
                         ID2 = cells,
                         ES = log2(observed_overlap+1) - log2(mean(null_overlaps) + 1),
                         Pvalue = p_value,
                         Sample = sample)
    allres <- rbind(allres, tmpres) 
  }
}

write.table(allres, file.path(subres, "Part1.2.spatial.TROP2.colocalization.immune.rctd.basal.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)

#######################immune cell colocalization (plot figures)#######################
library(ggplot2)
library(dplyr)
library(gghalves)
library(ggbeeswarm)
library(patchwork)
library(ggthemes)

allres <- fread(file.path(subres, "Part1.2.spatial.TROP2.colocalization.immune.rctd.luminal.txt"), header = T, stringsAsFactors = F, data.table = F)
allres <- allres %>%
  group_by(ID2) %>%
  mutate(median_ES = median(ES)) %>%
  ungroup() %>%
  arrange(median_ES) %>%
  mutate(ID2 = factor(ID2, levels = unique(ID2)))

# Create the plot with black borders
allres <- allres %>%
  group_by(ID2) %>%
  mutate(median_ES = median(ES)) %>%
  ungroup() %>%
  arrange(median_ES) %>%
  mutate(ID2 = factor(ID2, levels = unique(ID2))) %>%
  mutate(Significance = case_when(
    Pvalue < 0.05 & ES > 0 ~ "sigup",
    Pvalue < 0.05 & ES < 0 ~ "sigdown",
    TRUE ~ "nonsig"
  ),
  negLogP = -log10(Pvalue),
  Significance = factor(Significance, levels = c("sigdown", "nonsig", "sigup")))

custom_colors <- c("sigup" = "#CB181D", "sigdown" = "#2171B5", "nonsig" = "gray80")
# Step 2: Count Significant Up and Down for Each Cell Type
bar_counts <- allres %>%
  group_by(ID2, Significance) %>%
  summarise(Count = n(), .groups = "drop")

bar_counts$Significance <- factor(bar_counts$Significance, levels = c("sigdown", "nonsig", "sigup"))


# Step 3: Create the Scatter Plot with Updated Colors
scatter_plot <- ggplot(allres, aes(x = ES, y = ID2)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = Significance, size = negLogP),
    groupOnX = FALSE,
    width = 0.3,
    alpha = 0.75,
    shape = 16
  ) +
  scale_color_manual(values = custom_colors) +
  scale_size_continuous(range = c(0.4, 2.2)) +
  labs(
    x = "ES",
    y = "Cell Type",
    color = "Significance",
    size = "-log10(P value)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )

# Step 4: Create the Bar Plot for Significant Counts
bar_plot <- ggplot(bar_counts, aes(x = Count, y = ID2, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_manual(values = custom_colors) +
  theme_classic2() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Significant Counts")

# Step 5: Combine the Scatter Plot and Bar Plot
combined_plot <- scatter_plot + bar_plot + plot_layout(widths = c(1, 1.5))

pdf(file.path(subres, "V2503.figure4.spatial.immunecell.colocalization.luminal.pdf"), width = 6, height = 4)
print(combined_plot)
dev.off()


allres <- fread(file.path(subres, "Part1.2.spatial.TROP2.colocalization.immune.rctd.basal.txt"), header = T, stringsAsFactors = F, data.table = F)
allres <- allres %>%
  group_by(ID2) %>%
  mutate(median_ES = median(ES)) %>%
  ungroup() %>%
  arrange(median_ES) %>%
  mutate(ID2 = factor(ID2, levels = unique(ID2)))

# Create the plot with black borders
allres <- allres %>%
  group_by(ID2) %>%
  mutate(median_ES = median(ES)) %>%
  ungroup() %>%
  arrange(median_ES) %>%
  mutate(ID2 = factor(ID2, levels = unique(ID2))) %>%
  mutate(Significance = case_when(
    Pvalue < 0.05 & ES > 0 ~ "sigup",
    Pvalue < 0.05 & ES < 0 ~ "sigdown",
    TRUE ~ "nonsig"
  ),
  negLogP = -log10(Pvalue),
  Significance = factor(Significance, levels = c("sigdown", "nonsig", "sigup")))

custom_colors <- c("sigup" = "#CB181D", "sigdown" = "#2171B5", "nonsig" = "gray80")
# Step 2: Count Significant Up and Down for Each Cell Type
bar_counts <- allres %>%
  group_by(ID2, Significance) %>%
  summarise(Count = n(), .groups = "drop")

bar_counts$Significance <- factor(bar_counts$Significance, levels = c("sigdown", "nonsig", "sigup"))


# Step 3: Create the Scatter Plot with Updated Colors
scatter_plot <- ggplot(allres, aes(x = ES, y = ID2)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(color = Significance, size = negLogP),
    groupOnX = FALSE,
    width = 0.3,
    alpha = 0.75,
    shape = 16
  ) +
  scale_color_manual(values = custom_colors) +
  scale_size_continuous(range = c(0.4, 2.2)) +
  labs(
    x = "ES",
    y = "Cell Type",
    color = "Significance",
    size = "-log10(P value)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )

# Step 4: Create the Bar Plot for Significant Counts
bar_plot <- ggplot(bar_counts, aes(x = Count, y = ID2, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_manual(values = custom_colors) +
  theme_classic2() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Significant Counts")

# Step 5: Combine the Scatter Plot and Bar Plot
combined_plot <- scatter_plot + bar_plot + plot_layout(widths = c(1, 1.5))

pdf(file.path(subres, "V2503.figure4.spatial.immunecell.colocalization.basal.pdf"), width = 6, height = 4)
print(combined_plot)
dev.off()


#######################pathway changes gradiantly#######################
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
names(hall) <- c("Signature", "Genes")

kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)
names(kegg) <- c("Signature", "Genes")

metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cancer.mps.txt", header = T, stringsAsFactors = F, data.table = F)
metapros$ID <- "ID"
metapros <- melt(metapros, id.vars = "ID")
names(metapros) <- c("Class", "Signature", "Genes")
metapros <- metapros[, -1]
metapros <- rbind(metapros, hall, kegg)
allscores <- split(metapros$Genes, metapros$Signature)


usefulsams <- unique(spadat$BCID)
allmps <- data.frame()
for (sam in usefulsams) {
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  if (file.exists(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))){
    spadat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))
    spadat <- AddModuleScore(spadat,
                             features = allscores,
                             name="TCELLSIGNATURE")
    colnames(spadat@meta.data)[grep("TCELLSIGNATURE", colnames(spadat@meta.data))] <- names(allscores)
    
    tmpres <- as.data.frame(spadat@meta.data)
    tmpres <- tmpres[, c("ID", "BCID", "Distance", names(allscores))]
    
    allmps <- rbind(allmps, tmpres)
  }
}


allmps$Distance[allmps$Distance > 9] <- 9
allmps$Distance[allmps$Distance < -9] <- -9

allmps_melt <- melt(allmps, id.vars = c("ID", "BCID", "Distance"))

finalres_median <- allmps_melt %>%
  group_by(BCID, Distance, variable) %>%
  summarise(mean_frac = mean(value, na.rm = TRUE), .groups = "drop")


selpath <- c("HALLMARK_COMPLEMENT", "MP17 Interferon/MHC-II (I)", "KEGG_ECM_RECEPTOR_INTERACTION")
tmpres <- finalres_median[finalres_median$variable %in%selpath, ]

p <- ggplot(tmpres, aes(x = Distance, y = mean_frac)) +
  geom_smooth(method = "loess", se = TRUE, fill = "gray80", color = "orange") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  scale_x_continuous(breaks = c(-9, -6, -3, 0, 3, 6, 9)) +
  labs(
    x = "Distance",
    y = "Mean Fraction",
    color = "Class"
  ) +
  facet_wrap(~variable, ncol = 8, scales = "free_y") +
  theme_classic2()

pdf(file.path(subres, "V2503.figure4.spatial.distance.pathway.selected.pdf"), width = 8, height = 3)
print(p)
dev.off()



#######################R/NR different proportion#######################
tmpsubres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data


combined_info <- data.frame()
for (i in 1:nrow(allsams)) {
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  print(i)
  sample <- allsams$Sample[i]
  if (file.exists(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))) {
    resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
    seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))
    
    tmpres <- seudat@meta.data
    rownames(tmpres) <- paste(tmpres$BCID, rownames(tmpres), sep = ".")
    tmpres$Responder = allsams_id[allsams_id$ID == sample, ]$Responder
    combined_info <- rbind(combined_info, tmpres)
    
  }
}

finalres <- merge(allinfo, combined_info[, c("CellID", "Responder")], by.x = "CellID", by.y = "CellID")


chisq_test <- chisq.test(table(finalres$BLClass, finalres$Responder))
finalres$Responder <- factor(finalres$Responder, c("Yes", "No"))
finalres$BLClass <- factor(finalres$BLClass, c("Lum", "IM", "B.EMT", "B.IS"))
tmpres <- as.data.frame(prop.table(table(finalres$Responder, finalres$BLClass), margin = 1))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0"))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()


pdf(file.path(subres, "V2503.figure4.spatial.region.responders.pdf"), width = 5, height = 6)
print(p)
dev.off()

tmpres <- finalres[finalres$BCID == c("BC01", "BC17"), ]
chisq_test <- chisq.test(table(tmpres$BLClass, tmpres$BCID))

tbl <- as.data.frame(table(tmpres$BLClass, tmpres$BCID))

# Create pie plot for tum_c1
tum_c1 <- tbl %>% 
  filter(Var2 == "BC01") %>% 
  ggpie(x = "Freq", label = "Var1", 
        fill = "Var1", 
        color = "white", 
        palette = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0")) +
  ggtitle("BC01")+
  theme_void()

# Create pie plot for tum_c2
tum_c2 <- tbl %>% 
  filter(Var2 == "BC17") %>% 
  ggpie(x = "Freq", label = "Var1", 
        fill = "Var1", 
        color = "white", 
        palette = c(Lum = "#5C9E74", IM = "#A8DAB5", B.EMT = "#F1C7DA", B.IS = "#C973A0")) +
  ggtitle("BC17")+
  theme_void()


pdf(file.path(subres, "V2503.figure4.spatial.region.responders.pie.pdf"), width = 5, height = 3)
print(ggarrange(tum_c1, tum_c2, ncol = 2))
dev.off()



#######################T/N boundary gene expression#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data


alldemarkers <- data.frame()
for (j in 1:nrow(allsams)) {
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  print(j)
  sample <- allsams$Sample[j]
  if (file.exists(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))) {
    
    spadat1 <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))
    
    if ("Tumor_Bound" %in% unique(spadat1$Regions) & "Normal_Bound" %in% unique(spadat1$Regions)) {
      tmpspadat <- subset(spadat1, subset = Regions %in% c("Tumor_Bound", "Normal_Bound"))
      Idents(tmpspadat) <- tmpspadat$Regions
      markers <- FindAllMarkers(tmpspadat, only.pos = FALSE, min.pct = 0.0, logfc.threshold = 0.0)
      markers <- markers[markers$cluster == "Tumor_Bound", ]
      markers$logp <- -log10(markers$p_val_adj)
      markers$Class <- "NS"
      markers$Class <- ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.5, "Up", markers$Class)
      markers$Class <- ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC < 0.5, "Down", markers$Class)
      markers$Sample <- sample
      
      alldemarkers <- rbind(alldemarkers, markers)
    }
  }
}

write.table(alldemarkers, file.path(subres, "V2503.figure4.spatial.tnboundary.diff.gene.pathway.txt"), 
            row.names = F, col.names = T, sep = "\t", quote = F)



upgenes <- alldemarkers[alldemarkers$Class == "Up", ]
sort(table(upgenes$gene), decreasing = T)[1:50]
dngenes <- alldemarkers[alldemarkers$Class == "Down", ]
sort(table(dngenes$gene), decreasing = T)[1:50]


alldemarkers$Class <- ifelse(alldemarkers$p_val_adj < 0.05 & alldemarkers$avg_log2FC > 0.5, "Up", alldemarkers$Class)
alldemarkers$Class <- ifelse(alldemarkers$p_val_adj < 0.05 & alldemarkers$avg_log2FC < 0.5, "Down", alldemarkers$Class)
upgenes <- alldemarkers[alldemarkers$Class == "Up", ]
genes <- names(sort(table(upgenes$gene), decreasing = T)[sort(table(upgenes$gene), decreasing = T) >= 9])
length(genes)
sort(table(upgenes$gene), decreasing = T)["NECTIN4"]

fdaapprove <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/06.HPA/protein_class_FDA.tsv",
                    header = T, stringsAsFactors = F, data.table = F)
fdagenes <- c(fdaapprove$Gene, "NECTIN4")

length(genes)
length(fdagenes)
length(intersect(genes, fdagenes))


#######################spatial examples, region, class, expression#######################
tmpsubres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
allinfo <- spadat@meta.data


usefulsams <- unique(spadat$BCID)

pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.examples.pdf"), width = 10, height = 4)
for (sam in usefulsams) {
  
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  if (file.exists(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))){
    seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.distance.rds", sep = "")))
    gray_image <- apply(seudat@images$slice1@image[,,1:3], c(1,2), mean)
    seudat@images$slice1@image[,,1] <- gray_image
    seudat@images$slice1@image[,,2] <- gray_image
    seudat@images$slice1@image[,,3] <- gray_image
    
    metainfo <- seudat@meta.data
    rownames(metainfo) <- paste(metainfo$BCID, rownames(metainfo), sep = ".")
    tmpallinfo <- allinfo[rownames(metainfo), ]
    
    seudat$AssignedCombo <- tmpallinfo$BLClass
    
    
    tmpinfo <- seudat@meta.data
    tmpinfo <- tmpinfo[, c("Barcodes", "Original_Annot", "Level2_Annot", "AssignedCombo")]
    write.table(tmpinfo, "../Results/BC01.annotation.txt",
                row.names = F, col.names = T, sep = "\t", quote = F)
    
    tmpseudat <- subset(seudat, subset = AssignedCombo %in% c("Lum", "IM", "B.EMT", "B.IS"))
    tmpseudat$AssignedCombo <- ifelse(tmpseudat$AssignedCombo %in% c("Lum", "IM"), "Lum", "Basal")
    newsize = allsams[allsams$Sample == sample, ]$Size
    
    p2 <- SpatialDimPlot(tmpseudat, group.by = "AssignedCombo", image.alpha=0.2, crop = F, stroke = NA, pt.size.factor = newsize) + 
      scale_fill_manual(values = c(Lum = "#5C9E74", Basal = "#C973A0"))+
      ggtitle(sam) +
      theme(
        plot.title = element_text(size = 20),
        legend.position = "none"
      )
    
    seudat$Regions <- factor(seudat$Regions, c("Tumor_Core", "Tumor_Inter", "Tumor_Bound", "Normal_Bound", "Normal_Jaxta", "Normal_Distal"))
    p1 <- SpatialDimPlot(seudat, group.by = "Regions", image.alpha=0.2, crop = F, stroke = NA, pt.size.factor = newsize) + 
      scale_fill_manual(values = c(grad7))+
      ggtitle(sam) +
      theme(
        plot.title = element_text(size = 20),
        legend.position = "none"
      )
    
    p3 <- SpatialFeaturePlot(seudat, features = "NECTIN4", image.alpha = 0.2, crop = F, stroke = NA, pt.size.factor = newsize) +
      scale_fill_gradientn(colors = c("#E5E5E5", "#E5E5E5" ,grad2, "#EA71AE", "#EA71AE"))+
      ggtitle("basal") +
      theme(
        plot.title = element_text(size = 20)
      )
    
    
    print(p1|p2|p3)
    
  }
}

dev.off()





#######################NECTIN4 expression in TCGA, IM210, ABACUS, 2800#######################
###############TCGA###############
gene <- "PVRL4"

rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/gene_exp/TCGA.BLCA.sampleMap%2FHiSeqV2",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$sample
rnamtr <- rnamtr[, -1]
rnamtr <- rnamtr[, as.numeric(substr(names(rnamtr), 14, 15)) < 10]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/phenotype/TCGA.BLCA.sampleMap%2FBLCA_clinicalMatrix",
                  header = T, stringsAsFactors = F, data.table = F)

##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "sampleID")

finalres$Stage_Class <- finalres$pathologic_stage
finalres <- finalres[finalres$Stage_Class %in% c("Stage II", "Stage III", "Stage IV"), ]
finalres$Stage_Class <- ifelse(finalres$Stage_Class == "Stage IV", "Stage IV", "Stage II-III")

finalres$N_Class <- finalres$pathologic_N
finalres$N_Class <- ifelse(finalres$N_Class == "NX", "NA", finalres$N_Class)
finalres$N_Class <- ifelse(finalres$N_Class == "", "NA", finalres$N_Class)
finalres <- finalres[finalres$N_Class %in% c("N0", "N1", "N2", "N3"), ]
finalres$M_Class <- finalres$pathologic_M
finalres <- finalres[finalres$M_Class %in% c("M0", "M1"), ]

finalres$N_Class <- ifelse(finalres$M_Class == "M1", "M1", finalres$N_Class)
finalres$N_Class <- ifelse(finalres$N_Class %in% c("M1", "N1", "N2", "N3"), "m1-n1-3", finalres$N_Class)

finalres <- finalres[, c("bcr_patient_barcode", "Gene", "Stage_Class", "N_Class", "M_Class")]
top_1_3_cutoff <- quantile(finalres$Gene, probs = 3/4, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/4, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)

finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]
cellinfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/TCGA.cell.info.txt",
                  header = T, stringsAsFactors = F, data.table = F)
cellinfo <- cellinfo[, c("Case ID", "mRNA cluster")]
names(cellinfo) <- c("CaseID", "TCGACluster")

finalres <- merge(finalres, cellinfo, by.x = "bcr_patient_barcode", by.y = "CaseID")

library(ggalluvial)
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$TCGACluster))

tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$TCGACluster), margin = 1))
tmpres$Var2 <- factor(tmpres$Var2, c("Luminal_papillary", "Luminal", "Luminal_infiltrated", "Basal_squamous", "Neuronal"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0","#DEEBF7"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.tcga.pdf"), width = 7, height = 3)
print(p)
dev.off()


###############IM210###############
gene <- "PVRL4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/IMvigor210.geneexp.mtx",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/IMvigor210.cliinfo.txt",
                  header = T, stringsAsFactors = F, data.table = F)

##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "SampleID")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]


library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.im210.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "Lund2"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("UroA", "UroB", "Infiltrated", "Basal/SCC-like", "Genomically unstable"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0","#DEEBF7"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()


###############ABACUS###############
gene <- "NECTIN4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/ABACUS.geneexp.mtx",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/ABACUS.cliinfo.txt",
                  header = T, stringsAsFactors = F, data.table = F)

##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "SampleID")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]

library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.abacus.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "IMMUNE_PHENO"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("Excluded", "Desert", "Inflamed"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()


###############IM010###############
gene <- "NECTIN4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/wo29636_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <-saminfo
clininfo$NMF <- paste("NMF", clininfo$NMF, sep = "")
##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "trunc_anonymized_sample_ids")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]


library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.IM010.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("NMF1", "NMF2", "NMF3", "NMF4"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()

###############IM130###############
gene <- "NECTIN4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/wo30070_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <-saminfo
clininfo$NMF <- paste("NMF", clininfo$NMF, sep = "")
##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "trunc_anonymized_sample_ids")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]


library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.IM130.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("NMF1", "NMF2", "NMF3", "NMF4"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()

###############IM210###############
gene <- "NECTIN4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <-saminfo
clininfo$NMF <- paste("NMF", clininfo$NMF, sep = "")
##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "trunc_anonymized_sample_ids")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]


library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.IM210.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("NMF1", "NMF2", "NMF3", "NMF4"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()

###############IM211###############
gene <- "NECTIN4"
library(survival)
library(survminer)
rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29294_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

res_tfs <- data.frame(ID = names(rnamtr),
                      Gene = as.numeric(rnamtr[gene, ]))

##deal survival file
clininfo <-saminfo
clininfo$NMF <- paste("NMF", clininfo$NMF, sep = "")
##prepare file for plot
finalres <- merge(res_tfs, clininfo, by.x = "ID", by.y = "trunc_anonymized_sample_ids")

top_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)
bottom_1_3_cutoff <- quantile(finalres$Gene, probs = 1/2, na.rm = TRUE)

finalres$FinalGroup <- "Others"
finalres$FinalGroup <- ifelse(finalres$Gene >= top_1_3_cutoff, "1.High", finalres$FinalGroup)
finalres$FinalGroup <- ifelse(finalres$Gene <= bottom_1_3_cutoff, "0.Low", finalres$FinalGroup)
finalres <- finalres[finalres$FinalGroup %in% c("1.High", "0.Low"), ]


library(ggalluvial)
pdf(file.path(subres, "V2503.figure4.spatial.show.nectin4.IM211.pdf"), width = 7, height = 3)
finalres$Stage_Class <- finalres[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(finalres$FinalGroup, finalres$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(finalres$FinalGroup, finalres$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("NMF1", "NMF2", "NMF3", "NMF4"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74", "#A8DAB5", "#F1C7DA","#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()

#######################EV resistence UCUM3 (2 groups)#######################
rnamtr <- fread("../EVResisRNA/OE_Res_Low.matrix",
                header = T, stringsAsFactors = F, data.table = F)
rnamtr <- aggregate(. ~ symbol, data = rnamtr, FUN = mean, na.rm = TRUE)

rownames(rnamtr) <- rnamtr$symbol
rnamtr <- rnamtr[, -1]

res_tfs <- data.frame(
  ID = colnames(rnamtr),
  Group = c(rep("OE", 3), rep("HRes", 3), rep("LRes", 3))
)
rownames(res_tfs) <- res_tfs$ID

rnamtr <- round(rnamtr)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rnamtr,
                              colData = res_tfs,
                              design = ~ Group)

dds <- DESeq(dds)
norm_mtx <- counts(dds, normalized = TRUE)

res <- results(dds, contrast = c("Group", "LRes", "OE"))


de.lowvsoe <- fread("../EVResisRNA/LowvsOE.deg.txt",
                    header = T, stringsAsFactors = F, data.table = F)

overgenes <- de.lowvsoe[de.lowvsoe$log2FoldChange > 0 & de.lowvsoe$padj < 0.05, ]$symbol
dngenes <- de.lowvsoe[de.lowvsoe$log2FoldChange < 0 & de.lowvsoe$padj < 0.05, ]$symbol

hall <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
bclineage <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/bc.linerge.txt",
                   header = T, stringsAsFactors = F, data.table = F)
metapros <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/cancer.mps.txt", header = T, stringsAsFactors = F, data.table = F)
metapros$ID <- "ID"
metapros <- melt(metapros, id.vars = "ID")
names(metapros) <- c("Class", "Signature", "Gene")


basallist <- c(bclineage[bclineage$Signature %in% c("basal", "squamous-differentiation"), ]$Genes, "EGFR")
lumialist <- bclineage[bclineage$Signature == "luminal", ]$Genes


intersect(basallist, overgenes)
intersect(basallist, dngenes)
intersect(lumialist, overgenes)
intersect(lumialist, dngenes)

mhc2list <- metapros[metapros$Signature == "MP17 Interferon/MHC-II (I)", ]$Gene
intersect(mhc2list, overgenes)
intersect(mhc2list, dngenes)

genes <- intersect(mhc2list, overgenes)
selected_genes <- grep("^(IFI|ISG|STA|HLA)", genes, value = TRUE)

allgenes <- c("FGFR3", "EGFR", sort(selected_genes))


norm_mtx <- norm_mtx[, c(1:3, 7:9)]
res_tfs <- data.frame(
  ID = colnames(norm_mtx),
  Group = c(rep("OE", 3), rep("LRes", 3))
)
rownames(res_tfs) <- res_tfs$ID
selexpmtr <- norm_mtx[allgenes, ]
selexpmtr <- t(scale(t(selexpmtr)))

library(pheatmap)
tmpann_colors = list(
  Group = c(OE="#F5A889", LRes="#ACD6EC"))

#selexpmtr[selexpmtr > 1.5] <- 1.5
#selexpmtr[selexpmtr < -1.5] <- -1.5
pdf(file.path(subres, "V2503.figure4.umuc3.gene.heatmap.2groups.pdf"), width = 5, height = 7)
pheatmap(as.matrix(selexpmtr),
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         annotation_col = res_tfs,
         annotation_colors = tmpann_colors,
         color =colorRampPalette(c(heatgradiant, "#C85D4D", "#C85D4D"))(50),
         border_color = NA,
         fontsize = 10
)
dev.off()




