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

############################################################figure 6 part#############################################################
dir.create("../9.spaAnaly/NewCombined.Analy/06.tls")
subres <- "../9.spaAnaly/NewCombined.Analy/06.tls"
tmpsubres <- "../9.spaAnaly/NewCombined.Analy/04.subtype"

#######################show all samples with TLS#######################
tlsseurat <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.idendity.rds"))
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)

for (sam in unique(tlsseurat$BCID)) {
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  spadat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  p1 <- SpatialPlot(spadat, image.alpha=1, pt.size.factor = 0, crop = F) + 
    ggtitle(sam) +
    theme(
      plot.title = element_text(size = 40),
      legend.position = "none"
    )
  gray_image <- apply(spadat@images$slice1@image[,,1:3], c(1,2), mean)
  spadat@images$slice1@image[,,1] <- gray_image
  spadat@images$slice1@image[,,2] <- gray_image
  spadat@images$slice1@image[,,3] <- gray_image
  
  spadat$NewClass <- "Stroma"
  spadat$NewClass <- ifelse(spadat$Level2_Annot %in% c("Tumor Cell"), "Tumor Cell", spadat$NewClass)
  spadat$NewClass <- ifelse(spadat$Level2_Annot %in% c("TLS/TIL/LA"), "TLS/TIL/LA", spadat$NewClass)
  
  spadat$NewClass <- factor(spadat$NewClass, c("Tumor Cell", "Stroma", "TLS/TIL/LA"))
  color_v=c("#F59696", "#CCCCCC", "#FB40B0")
  
  
  allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
  # Normalize Responder column to Y/N
  allsams_id <- allsams_id %>%
    mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N"))
  tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))
  tlsdis_min <- tlsdis %>%
    group_by(TLSID) %>%
    slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(DisClass = ifelse(Distance < 24, "PT-TLS", "DT-TLS")) %>%
    select(TLSID, BCID, Distance, DisClass)
  # Add responder info from metadata
  tlsdis_min <- tlsdis_min %>%
    left_join(allsams_id %>% select(MyID, ResponderYN), by = c("BCID" = "MyID"))
  
  tlsdis_min$FinalClass <- NA
  tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "DT-TLS", "DT-TLS", tlsdis_min$FinalClass)
  tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "Y", "PT-TLS(R)", tlsdis_min$FinalClass)
  tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "N", "PT-TLS(NR)", tlsdis_min$FinalClass)
  
  tlsseurat <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.idendity.rds"))
  meta <- tlsseurat@meta.data
  meta <- meta %>%
    left_join(tlsdis_min %>% select(TLSID, FinalClass), by = "TLSID")
  
  
  spadat$NewClass <- as.character(spadat$NewClass)
  spadat$NewClass <- ifelse(spadat$NewClass == "TLS/TIL/LA", "Others", spadat$NewClass)
  spadat$NewClass <- ifelse(spadat$CellID %in% meta[meta$FinalClass == "DT-TLS", ]$CellID, "DT-TLS", spadat$NewClass)
  spadat$NewClass <- ifelse(spadat$CellID %in% meta[meta$FinalClass %in% c("PT-TLS(NR)", "PT-TLS(R)"), ]$CellID, "PT-TLS", spadat$NewClass)
  
  spadat$NewClass <- factor(spadat$NewClass, c("Tumor Cell", "Stroma", "PT-TLS", "DT-TLS", "Others"))
  color_v=c("#F59696", "#CCCCCC", "#619DB8", "#91C299", "gray95")
  
  p2 <- SpatialDimPlot(spadat, group.by = "NewClass", image.alpha=0.6, stroke = NA, pt.size.factor = allsams[allsams$Sample == sample, ]$Size) + 
    scale_fill_manual(values = color_v)+
    ggtitle(sample) +
    theme(
      plot.title = element_text(size = 40),
      legend.position = "right"
    )
  
  pdf(file.path(subres, paste(sample, ".V2503.figureS9.spatial.tls.example.pdf", sep = "")), width = 6, height = 4)
  print(p1)
  print(p2)
  dev.off()
  
}




#######################Classification of TLS#######################
#############distance class
tlsseurat <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.idendity.rds"))
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
spadat <- subset(spadat, subset = BLClass %in% c("Lum", "IM", "B.EMT", "B.IS"))

Idents(spadat) <- spadat$BLClass
markers <- FindAllMarkers(
  object = spadat,
  group.by = "BLClass",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE)

marker_list <- split(top_markers$gene, top_markers$cluster)

allres <- data.frame()
for (sam in unique(tlsseurat$BCID)) {
  print(sam)
  sample <- allsams_id[allsams_id$MyID == sam, ]$ID
  
  tmptls <- subset(tlsseurat, subset = BCID == sam)
  
  resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
  seudat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))
  
  tmpseudat <- subset(seudat, subset = Level2_Annot == "Tumor Cell")
  if (sam == "BC11") {
    tmpseudat <- AddModuleScore(
      object = tmpseudat,
      features = marker_list,
      name = names(marker_list)
    )
    score_cols <- paste0(names(marker_list), c("1", "2", "3", "4"))
    tmpseudat$BLClass <- apply(tmpseudat@meta.data[, score_cols], 1, function(x) {
      names(marker_list)[which.max(x)]
    })
    tmptumor <- subset(tmpseudat, subset = BCID == sam)
    
  }else{
    tmptumor <- subset(spadat, subset = BCID == sam)
    
  }
  
  spata_object <-
    asSPATA2(
      object = seudat,
      sample_name = sample,
      platform = "VisiumLarge",
      img_scale_fct = "lowres",
      assay_name = "Spatial",
      assay_modality = "gene"
    )
  coords <- as.data.frame(getCoordsDf(spata_object))
  
  for (i in unique(tmptls$TLSID)) {
    subtmptls <- subset(tmptls,subset = TLSID == i)
    for (j in unique(tmptumor$BLClass)) {
      subtmptumor <- subset(tmptumor,subset = BLClass == j)
      region_spots <- coords[coords$barcodes %in% subtmptls$ID, ]
      all_spots <- coords[coords$barcodes %in% subtmptumor$ID, ]
      all_points <- ppp(all_spots$x, all_spots$y, window=owin(range(coords$x), range(coords$y)))
      region_points <- ppp(region_spots$x, region_spots$y, window=owin(range(coords$x), range(coords$y)))
      # Calculate distances from each spot to the closest spot in the region
      distances <- nncross(all_points, region_points, what="dist")
      
      tmpres <- data.frame(TLSID = i,
                           Subtype = j,
                           Distance = min(distances),
                           Sample = sample,
                           BCID = sam)
      allres <- rbind(allres, tmpres)
    }
  }
}

saveRDS(allres, file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))


#######################TLS and basal luminal#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
#############distance
tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))

tlsdis <- tlsdis[tlsdis$Subtype != "IM", ]
tlsdis$Subtype <- ifelse(tlsdis$Subtype == "Lum", "Luminal", "Basal")

tlsdis_min <- tlsdis %>%
  group_by(TLSID, Subtype) %>%
  slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
  ungroup()

tlsdis_min <- as.data.frame(tlsdis_min)
tlsdis_min$Distance <- log2(tlsdis_min$Distance+1)
tlsdis_min$Responder <- ifelse(tlsdis_min$BCID %in% c("BC01", "BC11", "BC12", "BC14", "BC15"), "Y", "N")
#tlsdis_min$Distance <- tlsdis_min$Distance + runif(nrow(tlsdis_min), min = 1e-5, max = 1e-4)

tls_N <- tlsdis_min %>% filter(Responder == "N")
tls_Y <- tlsdis_min %>% filter(Responder == "Y")

my_comp <- list(c("Basal", "Luminal"))
p1 <- ggplot(tls_N, aes(x = Subtype, y = Distance, group = TLSID)) +
  geom_point(aes(color = Subtype), size = 3) +
  geom_line(size = 1) +
  stat_compare_means(comparisons = my_comp, method = "wilcox", paired = TRUE) +
  scale_color_manual(values = c("Basal" = "#EA71AE", "Luminal" = "#148F28")) +
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme_classic2() +
  labs(title = "Responder = N")

p2 <- ggplot(tls_Y, aes(x = Subtype, y = Distance, group = TLSID)) +
  geom_point(aes(color = Subtype), size = 3) +
  geom_line(size = 1) +
  stat_compare_means(comparisons = my_comp, method = "wilcox", paired = TRUE) +
  scale_color_manual(values = c("Basal" = "#EA71AE", "Luminal" = "#148F28")) +
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme_classic2() +
  labs(title = "Responder = Y")

p <- ggarrange(p1, p2, ncol = 2)


ggsave(file.path(subres, "V2503.figure5.spatial.tls.basal.distance.pdf"), p, width = 7, height = 3)



#############basal-luminal TLS stimulating marker expression
spadat <- readRDS(file.path(tmpsubres, "V2503.figure3.basal.luminal.classification.ISCHIA.umap.rds"))
selgenes <- c("CXCL13", "CCL19", "CCL21", "CXCL9", "CXCL10", "LTB", "LTA", 
              "TNFSF13", "TNFSF14", "TNFSF13B", "IL6", "ICAM1", "VCAM1", "CXCL11", "CCL5")
Idents(spadat) <- spadat$BLClass
avg_exp <- AverageExpression(spadat, features = selgenes, return.seurat = FALSE)$RNA
# Optional: scale expression per gene (row)
scaled_matrix <- t(scale(t(avg_exp)))  # z-score scaling by row
scaled_matrix <- scaled_matrix[selgenes, c("Lum", "IM", "B.EMT", "B.IS")]
# Set output file
pdf(file.path(subres, "V2503.figure5.spatial.class.tlsformation.pdf"), width = 8, height = 2)
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

pdf (file.path(subres, "V2503.figure5.spatial.class.tlsformation.dotplot.pdf"), height = 4, width = 12)
print(DotPlot(spadat, features = unique(selgenes))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete("")+
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        scale_colour_gradient2(low = "white", mid = "#7C93C3", high = "#750E21", midpoint = 1.0) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()



#######################Distribution of TLS#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
# Normalize Responder column to Y/N
allsams_id <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N"))


tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))

tlsdis_min <- tlsdis %>%
  group_by(TLSID) %>%
  slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(DisClass = ifelse(Distance < 24, "PT-TLS", "DT-TLS")) %>%
  select(TLSID, BCID, Distance, DisClass)

# Add responder info from metadata
tlsdis_min <- tlsdis_min %>%
  left_join(allsams_id %>% select(MyID, ResponderYN), by = c("BCID" = "MyID"))

# Construct full combination of all patients and TLS classes
all_patients <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N")) %>%
  select(BCID = MyID, ResponderYN)

all_combinations <- expand.grid(
  BCID = all_patients$BCID,
  DisClass = c("PT-TLS", "DT-TLS"),
  stringsAsFactors = FALSE
)

# Compute true TLS summary
tls_summary_real <- tlsdis_min %>%
  count(BCID, DisClass)

# Join and fill missing values with zero count
tls_summary <- all_combinations %>%
  left_join(tls_summary_real, by = c("BCID", "DisClass")) %>%
  mutate(n = replace_na(n, 0)) %>%
  left_join(all_patients, by = "BCID")

# Sort samples: Responder Y first, then N
tls_summary <- tls_summary %>%
  mutate(ResponderYN = factor(ResponderYN, levels = c("Y", "N"))) %>%
  arrange(ResponderYN, BCID)

tls_summary$BCID <- factor(tls_summary$BCID, levels = unique(tls_summary$BCID))

# Define background color by responder group
responder_bg <- tls_summary %>%
  select(BCID, ResponderYN) %>%
  distinct() %>%
  mutate(x = as.numeric(BCID),
         xmin = x - 0.5,
         xmax = x + 0.5,
         bg_fill = ifelse(ResponderYN == "Y", "#d0ecff", "#ffe6e6"))

# Plot
p <- ggplot() +
  geom_rect(data = responder_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = ResponderYN),
            alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_bar(data = tls_summary,
           aes(x = BCID, y = n, fill = DisClass),
           stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c(
    "PT-TLS" = "#C973A0",
    "DT-TLS" = "#5C9E74",
    "Y" = "#d0ecff",
    "N" = "#ffe6e6"
  ))+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 7.5))+
  labs(title = "TLS Distribution per Patient",
       x = "Patient (BCID)", y = "Count",
       fill = "TLS Classification") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(subres, "V2503.figure5.spatial.tls.sample.distribu.pdf"), p, width = 6, height = 3)


#######################TLS to responder/Nonresponder (3 class)#######################
allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
# Normalize Responder column to Y/N
allsams_id <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N"))

# Read TLS minimum distance data
tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))

tlsdis_min <- tlsdis %>%
  group_by(TLSID) %>%
  slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(DisClass = ifelse(Distance < 24, "PT-TLS", "DT-TLS")) %>%
  select(TLSID, BCID, Distance, DisClass)

# Add responder info from metadata
tlsdis_min <- tlsdis_min %>%
  left_join(allsams_id %>% select(MyID, ResponderYN), by = c("BCID" = "MyID"))

tlsdis_min <- unique(tlsdis_min[, c("BCID", "DisClass", "ResponderYN")])

all_patients <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N")) %>%
  select(BCID = MyID, ResponderYN)

all_patients$DisClass <- "Non-TLS"
all_patients <- all_patients[!(all_patients$BCID %in% tlsdis_min$BCID), ]

alldata <- rbind(tlsdis_min, all_patients)

alldata$DisClass <- factor(alldata$DisClass, c("PT-TLS", "DT-TLS", "Non-TLS"))
alldata_filtered <- alldata %>%
  group_by(BCID) %>%
  arrange(DisClass) %>%  # PT-TLS < DT-TLS alphabetically
  slice(1) %>%
  ungroup()

newinfo <- as.data.frame(alldata_filtered)
her2_colors <- c(
  "PT-TLS" = "#C973A0",
  "DT-TLS" = "#5C9E74",
  "Non-TLS" = "#FAE7D9"   # Light Blue
)

prim_site_colors <- c(
  "Y" = "#d0ecff",
  "N" = "#ffe6e6"
)

newinfo$ResponderYN <- factor(newinfo$ResponderYN, c("Y", "N"))
newinfo$RandomID <- 1:1
# Create the Sankey plot with enhanced coloring
p1 <- ggplot(newinfo, aes(axis1 = DisClass, axis2 = ResponderYN, y = RandomID)) +
  geom_alluvium(aes(fill = DisClass), width = 1/10, alpha = 0.8) +  # Alluvium colored by HER2.IHC
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "black") +  # Strata colored
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, hjust = 0.5) +  # Text labels
  scale_x_discrete(limits = c("DisClass", "ResponderYN"), expand = c(0.2, 0.2)) + # More space
  scale_fill_manual(values = c(her2_colors, prim_site_colors)) +  # Custom colors
  labs(
    title = "Flow of DisClass to ResponderYN",
    x = "Categories", y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # Zero background
    panel.grid = element_blank(),                               # No gridlines
    axis.text.x = element_text(size = 12, face = "bold"),       # X-axis text
    axis.title = element_text(size = 14, face = "bold")         # Axis titles
  )


M <- table(newinfo$ResponderYN, newinfo$DisClass)
data_matrix <- M
p_value_matrix <- matrix(NA, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
rownames(p_value_matrix) <- rownames(data_matrix)
colnames(p_value_matrix) <- colnames(data_matrix)
for (i in 1:nrow(data_matrix)) {
  for (j in 1:ncol(data_matrix)) {
    # Observed value: the specific cell
    observed <- data_matrix[i, j]
    
    # Remaining values: total matrix excluding the observed cell
    remaining <- data_matrix
    remaining[i, j] <- 0  # Set the observed cell to 0
    
    # Create a 2x2 contingency table
    contingency_table <- matrix(c(observed,
                                  sum(data_matrix[i, ]) - observed,  # Rest of the row
                                  sum(data_matrix[, j]) - observed,  # Rest of the column
                                  sum(data_matrix) - sum(data_matrix[i, ]) - sum(data_matrix[, j]) + observed),  # Rest of the matrix
                                nrow = 2, byrow = TRUE)
    
    # Perform the chi-square test
    p_value_matrix[i, j] <- fisher.test(contingency_table)$p.value
  }
}

p_value_matrix# Output the results

M_matrix <- prop.table(table(newinfo$ResponderYN, newinfo$DisClass), margin = 2)
propM <- as.data.frame(M_matrix)
propP <- as.data.frame(as.table(p_value_matrix))

# Get the raw counts for each combination
count_matrix <- as.data.frame(table(newinfo$ResponderYN, newinfo$DisClass))
colnames(count_matrix) <- c("Var1", "Var2", "Count")

# Merge counts into proportion data
propM$Var1 <- factor(propM$Var1, sort(unique(propM$Var1)))
propM$Pvalue <- propP$Freq
propM <- left_join(propM, count_matrix, by = c("Var1", "Var2"))

# Calculate cumulative frequency and p-value labels
propM <- propM %>%
  group_by(Var2) %>%
  arrange(desc(Var1)) %>%
  mutate(
    Midpoint_Freq = cumsum(Freq) - Freq / 2,
    Pvalue_label = ifelse(Pvalue < 0.06, "*", "")
  )

# Add total counts for each HER2.IHC group
her2_counts <- as.data.frame(colSums(table(newinfo$ResponderYN, newinfo$DisClass)))
her2_counts$Var2 <- rownames(her2_counts)
colnames(her2_counts) <- c("Total_Count", "Var2")

# Merge total counts with plot data
propM <- left_join(propM, her2_counts, by = "Var2")

# Perform Chi-square test
pva <- chisq.test(table(newinfo$ResponderYN, newinfo$DisClass))

# Create the plot
p2 <- ggplot(propM, aes(fill = Var1, y = Freq, x = Var2)) +
  geom_flow(aes(alluvium = Var1, stratum = Var1), alpha = 0.8,
            width = 0.6, knot.pos = 0.5) +
  scale_fill_manual(values = prim_site_colors) +
  geom_col(aes(fill = Var1), width = 0.6, color = "black") +
  
  # Add p-value annotation
  geom_text(aes(y = Midpoint_Freq - 0.04, label = Pvalue_label),
            color = "red", size = 4, vjust = -0.2) +
  
  # Add total counts at the top of each bar
  geom_text(aes(y = 1.05, label = paste0("n=", Total_Count)), 
            color = "black", size = 4) +
  
  # Add individual counts inside each bar
  geom_text(aes(y = Midpoint_Freq, label = Count), 
            color = "black", size = 3.5) +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +  # Adjust y-axis to fit total count labels
  xlab(paste("Overall Chi-square p-value = ", round(pva$p.value, 4))) +
  theme_classic()

# Save the plot to a PDF
pdf(file.path(subres, "V2503.figure5.spatial.tls.sample.flow.3class.pdf"), width = 10, height = 4)
print(p1 | p2)
dev.off()


#######################TLS score#######################
tlsseurat <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.idendity.rds"))
tlsmap <- tlsseurat@meta.data


allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
# Normalize Responder column to Y/N
allsams_id <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N"))
tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))
tlsdis_min <- tlsdis %>%
  group_by(TLSID) %>%
  slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(DisClass = ifelse(Distance < 24, "PT-TLS", "DT-TLS")) %>%
  select(TLSID, BCID, Distance, DisClass)
# Add responder info from metadata
tlsdis_min <- tlsdis_min %>%
  left_join(allsams_id %>% select(MyID, ResponderYN), by = c("BCID" = "MyID"))

tlsdis_min$FinalClass <- NA
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "DT-TLS", "DT-TLS", tlsdis_min$FinalClass)
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "Y", "PT-TLS(R)", tlsdis_min$FinalClass)
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "N", "PT-TLS(NR)", tlsdis_min$FinalClass)

tlsmap_summary <- aggregate(TLS_Mean ~ TLSID, data = tlsmap, FUN = mean, na.rm = TRUE)

tlsdis_min <- merge(tlsdis_min, tlsmap_summary, by.x= "TLSID", by.y = "TLSID")

tlsdis_min$FinalClass <- factor(tlsdis_min$FinalClass, c("DT-TLS", "PT-TLS(R)", "PT-TLS(NR)"))

my_comp <- list(c("DT-TLS", "PT-TLS(R)"), c("PT-TLS(NR)", "PT-TLS(R)"), c("DT-TLS", "PT-TLS(NR)"))
p1 <- ggplot(tlsdis_min, aes(x = FinalClass, y = TLS_Mean, fill = FinalClass)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "TLS score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = c("#91C299", "#F5A889", "#ACD6EC"))+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")
#facet_wrap( ~ variable, scales = "free_y", ncol = 5)

ggsave(file.path(subres, "V2503.figure5.spatial.tls.score.pdf"), p1, width = 4, height = 3.5)


#######################Markers of TLS#######################
####get gene list
tlsmarker <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/09.tlsmarker/TLS_related_markers/all_genes.csv",
                   header = T, stringsAsFactors = F, data.table = F)

tlsmarker_unique <- tlsmarker[, c("gene", "celltype")] |> 
  dplyr::distinct(gene, .keep_all = TRUE)

tlsmarker_unique <- tlsmarker_unique[order(tlsmarker_unique$celltype), ]
tlsmarker_unique$gene <- gsub("\\.", "-", tlsmarker_unique$gene)

tlsmarker_unique$celltype <- factor(tlsmarker_unique$celltype, c("B", "NaiveB", "MemB", "AtypicalMemB", "Breg",
                                                                 "GCB", "PC", "T", "Tn", "Tcm",
                                                                 "Teff", "Tfh", "Th17", "Tisg", "Treg",
                                                                 "Macrophage", "DC", "FDC", "apCAF", "CAF",
                                                                 "myCAF", "Fibroblasts", "Proliferation", "TCR signaling", "Exhaustion",
                                                                 "Cytokines and Chemokines", "Stress"))

tlsmarker_unique <- tlsmarker_unique[order(tlsmarker_unique$celltype), ]

allsams_id <- fread("sample.clinic.info", header = T, data.table = F, stringsAsFactors = F)
# Normalize Responder column to Y/N
allsams_id <- allsams_id %>%
  mutate(ResponderYN = ifelse(Responder == "Yes", "Y", "N"))
tlsdis <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.toTumorDistance.rds"))
tlsdis_min <- tlsdis %>%
  group_by(TLSID) %>%
  slice_min(order_by = Distance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(DisClass = ifelse(Distance < 24, "PT-TLS", "DT-TLS")) %>%
  select(TLSID, BCID, Distance, DisClass)
# Add responder info from metadata
tlsdis_min <- tlsdis_min %>%
  left_join(allsams_id %>% select(MyID, ResponderYN), by = c("BCID" = "MyID"))

tlsdis_min$FinalClass <- NA
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "DT-TLS", "DT-TLS", tlsdis_min$FinalClass)
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "Y", "PT-TLS(R)", tlsdis_min$FinalClass)
tlsdis_min$FinalClass <- ifelse(tlsdis_min$DisClass == "PT-TLS" & tlsdis_min$ResponderYN == "N", "PT-TLS(NR)", tlsdis_min$FinalClass)


tlsseurat <- readRDS(file.path(subres, "V2503.figure5.spatial.tls.idendity.rds"))

Idents(tlsseurat) <- tlsseurat$TLSID
geneexp <- AverageExpression(tlsseurat, features = tlsmarker_unique$gene, return.seurat = FALSE)$RNA

library(Matrix)
library(dplyr)
tls_map <- tlsdis_min %>% select(TLSID, FinalClass)
intersect_tls <- intersect(colnames(geneexp), tls_map$TLSID)
geneexp_sub <- geneexp[, intersect_tls]
class_vec <- tls_map$FinalClass[match(colnames(geneexp_sub), tls_map$TLSID)]
names(class_vec) <- colnames(geneexp_sub)
grouped_means <- lapply(unique(class_vec), function(cls) {
  cols_in_class <- which(class_vec == cls)
  if (length(cols_in_class) > 1) {
    rowMeans(geneexp_sub[, cols_in_class])
  } else {
    as.vector(geneexp_sub[, cols_in_class])  # 只有一个列时
  }
})
mean_matrix <- do.call(cbind, grouped_means)
colnames(mean_matrix) <- unique(class_vec)
rownames(mean_matrix) <- rownames(geneexp)
mean_matrix <- Matrix(mean_matrix, sparse = TRUE)
mean_matrix[1:5, ]

tlsmarker_unique <- tlsmarker_unique[tlsmarker_unique$gene %in% rownames(mean_matrix), ]
mean_matrix <- mean_matrix[tlsmarker_unique$gene, ]

gene_celltype <- tlsmarker_unique$celltype[match(rownames(mean_matrix), tlsmarker_unique$gene)]
names(gene_celltype) <- rownames(mean_matrix)

scaled_matrix <- t(scale(t(mean_matrix)))

pdf(file.path(subres, "V2503.figure5.spatial.class.tlsmarkers.pdf"), width = 60, height = 2)

ht <- Heatmap(
  t(scaled_matrix),
  name = "Scaled\nExpression",
  col = colorRamp2(c(-2, 0, 1), c("#619DB8", "#E3EEEF", "#C85D4D")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_split = gene_celltype,
  column_names_rot = 45,
  column_title_gp = gpar(fill = "#A8DAB5", col = NA, fontsize = 10),
  column_gap = unit(2, "mm")
)

draw(ht)
decorate_heatmap_body("Scaled\nExpression", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1.5))
})
dev.off()


pdf(file.path(subres, "V2503.figure5.spatial.class.tlsmarkers.text.pdf"), width = 60, height = 2)

ht <- Heatmap(
  t(scaled_matrix),
  name = "Scaled\nExpression",
  col = colorRamp2(c(-2, 0, 1), c("#619DB8", "#E3EEEF", "#C85D4D")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_split = gene_celltype,
  column_names_rot = 45,
  column_gap = unit(2, "mm")
)

draw(ht)
decorate_heatmap_body("Scaled\nExpression", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1.5))
})
dev.off()


genes_to_plot <- rownames(geneexp)
genes_to_plot <- c("CD79B", "CD3E", "BCL2", "CR2", "IFI44", "ISG15", "HSPB1", "JCHAIN", "IKZF2")

expr_long <- do.call(rbind, lapply(genes_to_plot, function(gene) {
  data.frame(
    Gene = gene,
    TLSID = colnames(geneexp),
    Expression = as.vector(geneexp[gene, ])
  )
}))

expr_long <- expr_long %>%
  left_join(tlsdis_min[, c("TLSID", "FinalClass")], by = "TLSID") %>%
  filter(!is.na(FinalClass))

setids <- data.frame(ID = genes_to_plot,
                     GeneName = c("CD79B", "CD3E", "BCL6", "CR2", "IFI44", "ISG15", "HSPB1", "CTLA4", "TIGIT"))
expr_long <- merge(expr_long, setids, by.x = "Gene", by.y = "ID")

expr_long$Gene <- factor(expr_long$Gene, genes_to_plot)
expr_long$FinalClass <- factor(expr_long$FinalClass, c("DT-TLS", "PT-TLS(R)", "PT-TLS(NR)"))
my_comp <- list(c("DT-TLS", "PT-TLS(R)"), c("PT-TLS(NR)", "PT-TLS(R)"), c("DT-TLS", "PT-TLS(NR)"))
p1 <- ggplot(expr_long, aes(x = FinalClass, y = Expression, fill = FinalClass)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "TLS score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = c("#91C299", "#ACD6EC", "#F5A889"))+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")+
  facet_wrap( ~ Gene, scales = "free_y", ncol = 9)

ggsave(file.path(subres, "V2503.figure5.spatial.class.tlsmarkers.boxplot.selected.pdf"), p1, width = 15, height = 6)


#######################plot single tls#######################
##in python scripts



































