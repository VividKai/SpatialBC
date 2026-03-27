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

#############################################################figure 5 part#############################################################
dir.create("../9.spaAnaly/NewCombined.Analy/07.public")
subres <- "../9.spaAnaly/NewCombined.Analy/07.public"

selgenes <- c("NF2", "TP53", "RB1", "EGFR", "FGFR3")
library(ComplexHeatmap)
library(circlize)
###############Inhouse (all)###############
newinfo  <- fread("newinfo.txt", header = T, data.table = F, stringsAsFactors = F)

allsams_id <- fread("rna_id.sample.info.txt", header = T, data.table = F, stringsAsFactors = F)
allsams_id <- merge(allsams_id, newinfo, by.x = "ID", by.y = "ID")

allsams <- data.frame(ID = c(allsams_id$ID, allsams_id$PostID),
                      TimePoint = rep("Pre", nrow(allsams_id)),
                      Responder = allsams_id$DownStaging,
                      Stage = allsams_id$Stage)
allsams <- na.omit(allsams)

tcga_mtr <- fread("../BulkRawData/02.MergeQuant/TPM_symbol.txt",
                  header = T, data.table = F, stringsAsFactors = F)
rownames(tcga_mtr) <- tcga_mtr$V1
tcga_mtr <- tcga_mtr[, -1]
tcga_mtr <- tcga_mtr[, allsams$ID]

exp_data <- log2(tcga_mtr + 1)
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))
basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id
base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]
# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]
# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)
# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- allsams
refdat$Basal <- basal_score
refdat$Luminal <- luminal_score
refdat$B2L <- basal_score - luminal_score
refdat$Classification <- classification
refdat$EGFRExp <- as.numeric(tcga_mtr["EGFR", ])
refdat$EGFRExp <- ifelse(refdat$EGFRExp > 50, 50, refdat$EGFRExp)
refdat$FGFR3Exp = as.numeric(tcga_mtr["FGFR3", ])
refdat$FGFR3Exp <- ifelse(refdat$FGFR3Exp > 50, 50, refdat$FGFR3Exp)
refdat$N4Exp = as.numeric(tcga_mtr["NECTIN4", ])
refdat$N4Class = ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

cnvmtr <- fread("../BulkRawData/04.GISTIC/all_thresholded.by_genes.txt", 
                header = T, stringsAsFactors = F, data.table = F)
names(cnvmtr) <- gsub(".sorted_markdup", "", names(cnvmtr))
mutmtr <- fread("../BulkRawData/02.Mutect2/Mutect2/merged_funcotated.final.maf", 
                header = T, stringsAsFactors = F, data.table = F)


# Format CNV matrix
cnv_mat <- cnvmtr
rownames(cnv_mat) <- cnv_mat$`Gene Symbol`
cnv_mat <- cnv_mat[ , !(colnames(cnv_mat) %in% c("Gene Symbol", "Locus ID", "Cytoband"))]
# Mutation dictionary
mut_dict <- split(mutmtr$Tumor_Sample_Barcode, mutmtr$Hugo_Symbol)
# Limit to samples present in refdat
valid_samples <- refdat$ID


for (gene in selgenes) {
  # CNV-altered samples (±2)
  cnv_ga <- if (gene %in% rownames(cnv_mat)) {
    sample_names <- colnames(cnv_mat)
    samples <- sample_names[cnv_mat[gene, ] %in% c(2, -2)]
    intersect(samples, valid_samples)
  } else character(0)
  
  # Mutation-altered samples
  mut_ga <- if (gene %in% names(mut_dict)) {
    intersect(mut_dict[[gene]], valid_samples)
  } else character(0)
  
  # Combined GA samples restricted to valid ones
  ga_samples <- unique(c(cnv_ga, mut_ga))
  
  # Label refdat
  refdat[, gene] <- ifelse(refdat$ID %in% ga_samples, "GA", "WT")
  refdat[, gene] <- ifelse(refdat$ID %in% colnames(cnv_mat), refdat[, gene], "NA")
}

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
res_color <- c("N" = "#F5A889", "Y" = "#ACD6EC", "NA" = "white")
sta_color <- c("T2" = "#91C299", "T3" = "#ACD6EC", "T4" = "#F5A889")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  Stage = refdat$Stage,
  Responder = refdat$Responder,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    Stage = sta_color,
    Responder = res_color,
    Stage = sta_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#FF9900", "#FF9900")),
    FGFR3Exp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#7D1EE5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.Inhouse.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()


tmprefdat <- refdat
tmprefdat$TimePoint <- factor(tmprefdat$TimePoint, c("Pre", "Post"))
my_comp <- list(c("Post", "Pre"))
p1 <- ggplot(tmprefdat, aes(x = TimePoint, y = B2L, fill = TimePoint)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tp_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.Inhouse.pre.post.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat[refdat$Responder != "NA", ]
tmprefdat$Responder <- factor(tmprefdat$Responder, c("Y", "N"))
my_comp <- list(c("N", "Y"))
p1 <- ggplot(tmprefdat, aes(x = Responder, y = B2L, fill = Responder)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = res_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.Inhouse.response.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()

library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.Inhouse.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "Classification"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Classification)*5)

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("Luminal-like", "Basal-like"))
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = c("#5C9E74","#C973A0"))+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()



###############TCGA###############
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/cnv/TCGA.BLCA.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$`Gene Symbol`
cnvmtr <- cnvmtr[, -1]

mutmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/mc3.v0.2.8.PUBLIC.maf", header = T, stringsAsFactors = F, data.table = F)
mutmtr$Tumor_Sample_Barcode <- substr(mutmtr$Tumor_Sample_Barcode, 1, 15)
mutmtr <- mutmtr[mutmtr$Variant_Classification == "Missense_Mutation", ]

mutsaminfo <- data.frame(Mut_EGFR = ifelse(names(cnvmtr) %in% mutmtr[mutmtr$Hugo_Symbol == "EGFR", ]$Tumor_Sample_Barcode,
                                           "Mut", "Others"),
                         Mut_FGFR3 = ifelse(names(cnvmtr) %in% mutmtr[mutmtr$Hugo_Symbol == "FGFR3", ]$Tumor_Sample_Barcode,
                                            "Mut", "Others"))

saminfo <- data.frame(NewID = names(cnvmtr))

for (gene in selgenes) {
  cnv_col <- ifelse(as.numeric(cnvmtr[gene, ]) %in% c(2,-2), "Amp", "Others")
  saminfo[[paste0("CNV_", gene)]] <- cnv_col
  mut_col <- ifelse(names(cnvmtr) %in% mutmtr[mutmtr$Hugo_Symbol == gene, ]$Tumor_Sample_Barcode,
                    "Mut", "Others")
  saminfo[[paste0("Mut_", gene)]] <- mut_col
  saminfo[[gene]] <- ifelse(cnv_col == "Amp" | mut_col == "Mut", "GA", "WT")
}

cellinfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/TCGA.cell.info.txt",
                  header = T, stringsAsFactors = F, data.table = F)
cellinfo <- cellinfo[, c("Case ID", "mRNA cluster")]
names(cellinfo) <- c("CaseID", "TCGACluster")
cellinfo$CaseID <- paste(cellinfo$CaseID, "-01", sep = "")

allinfo <- merge(saminfo, cellinfo, by.x = "NewID", by.y = "CaseID")

clininfo <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA/phenotype/TCGA.BLCA.sampleMap%2FBLCA_clinicalMatrix",
                  header = T, stringsAsFactors = F, data.table = F)
##prepare file for plot
clininfo$Stage_Class <- clininfo$pathologic_stage
allinfo <- merge(allinfo, clininfo[, c("sampleID", "Stage_Class")], by.x = "NewID", by.y = "sampleID")


tcga_mtr <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/03.TCGA_GDC/gene_exp/TCGA-BLCA.htseq_fpkm.gene.uniq.tsv",
                  header = T, data.table = F, stringsAsFactors = F)
rownames(tcga_mtr) <- tcga_mtr$gene
tcga_mtr <- tcga_mtr[, -1]
tcga_mtr <- tcga_mtr[, as.numeric(substr(names(tcga_mtr), 14, 15)) < 10]

exp_data <- tcga_mtr
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(tcga_mtr["EGFR", ]),
                     FGFR3Exp = as.numeric(tcga_mtr["FGFR3", ]),
                     N4Exp = as.numeric(tcga_mtr["PVRL4", ]))
refdat$Sample <- substr(refdat$Sample, 1, 15)
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "NewID")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
tp_color <- c("Post" = "#F5A889", "Pre" = "#ACD6EC")
res_color <- c("N" = "#F5A889", "Y" = "#ACD6EC", "NA" = "white")
sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("Luminal_papillary" = "#5C9E74", "Luminal" = "#A8DAB5", "Luminal_infiltrated" = "#F1C7DA", 
                "Basal_squamous" = "#C973A0", "Neuronal" = "#DEEBF7")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  TCGACluster = refdat$TCGACluster,
  Stage = refdat$Stage_Class,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    TCGACluster = tcga_color,
    Stage = sta_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.TCGA.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat
tmprefdat$TCGACluster <- factor(tmprefdat$TCGACluster, c("Luminal_papillary", "Luminal", "Luminal_infiltrated", 
                                                         "Basal_squamous", "Neuronal"))
my_comp <- list(c("Luminal_infiltrated", "Basal_squamous"), c("Luminal_papillary", "Basal_squamous"))
p1 <- ggplot(tmprefdat, aes(x = TCGACluster, y = B2L, fill = TCGACluster)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.TCGA.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat[refdat$Stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"), ]
tmprefdat$Stage <- factor(tmprefdat$Stage, c("Stage I", "Stage II", "Stage III", "Stage IV"))
my_comp <- list(c("Stage II", "Stage I"), c("Stage II", "Stage III"), c("Stage I", "Stage III"))
p1 <- ggplot(tmprefdat, aes(x = Stage, y = B2L, fill = Stage)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = sta_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.TCGA.stage.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.TCGA.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "TCGACluster"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("Luminal_papillary", "Luminal", "Luminal_infiltrated", "Basal_squamous", "Neuronal"))
# Stacked + percent
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
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/anon_abacus.fmi.data.csv", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$alias
#cnvmtr <- cnvmtr[, -c(1)]
#names(cnvmtr) <- gsub("FMI_", "", names(cnvmtr))
for (gene in colnames(cnvmtr)[-1]) {
  cnv_raw <- cnvmtr[, gene]
  cnv_status <- ifelse(is.na(cnv_raw) == FALSE, "GA", "WT")
  cnvmtr[[gene]] <- cnv_status
}

rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/ABACUS.geneexp.mtx",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/0.data/ABACUS.cliinfo.raw.txt",
                 header = T, data.table = F, stringsAsFactors = F)
allinfo <- merge(saminfo, cnvmtr, by.x = "alias", by.y = "alias")

exp_data <- rnamtr[, names(rnamtr) %in% saminfo$alias, ]
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(exp_data["EGFR", ]),
                     FGFR3Exp = as.numeric(exp_data["FGFR3", ]),
                     N4Exp = as.numeric(exp_data["NECTIN4", ]))
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "alias")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
tp_color <- c("POST" = "#F5A889", "PRE" = "#ACD6EC")
res_color <- c("No" = "#F5A889", "Yes" = "#ACD6EC")
#sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("UroA" = "#5C9E74", "UroB" = "#A8DAB5", "GU" = "#F1C7DA", 
                "SCCL" = "#C973A0", "Inf" = "#E5E5E5")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  Lund_5CAT = refdat$Lund_5CAT,
  VISIT = refdat$VISIT,
  PCR = refdat$PCR,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    Lund_5CAT = tcga_color,
    VISIT = tp_color,
    PCR = res_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(-0.5, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.ABACUS.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat[refdat$Stage != "Inf", ]
tmprefdat$Lund_5CAT <- factor(tmprefdat$Lund_5CAT, c("UroA", "UroB", "GU", 
                                                     "SCCL", "Neuronal"))
my_comp <- list(c("UroA", "SCCL"), c("UroB", "SCCL"))
p1 <- ggplot(tmprefdat, aes(x = Lund_5CAT, y = B2L, fill = Lund_5CAT)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.ABACUS.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
tmprefdat$VISIT <- factor(tmprefdat$VISIT, c("PRE", "POST"))
my_comp <- list(c("PRE", "POST"))
p1 <- ggplot(tmprefdat, aes(x = VISIT, y = B2L, fill = VISIT)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tp_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.ABACUS.pre.post.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat[refdat$Lund_5CAT != "Inf", ]
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.ABACUS.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "Lund_5CAT"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2, c("UroA", "UroB", "GU", 
                                     "SCCL", "Neuronal"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = tcga_color)+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()


###############IM010###############
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_go29294_wo30070_wo29636_FMI_binary_calls.csv", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$truncated_anonymized_sampleid
cnvmtr <- cnvmtr[, -c(2, 3)]
names(cnvmtr) <- gsub("FMI_", "", names(cnvmtr))
for (gene in colnames(cnvmtr)[-1]) {
  cnv_raw <- cnvmtr[, gene]
  cnv_status <- ifelse(is.na(cnv_raw), "NA",
                       ifelse(cnv_raw == 1, "GA", "WT"))
  cnvmtr[[gene]] <- cnv_status
}

head(cnvmtr)

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/wo29636_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)

allinfo <- merge(saminfo, cnvmtr, by.x = "trunc_anonymized_sample_ids", by.y = "truncated_anonymized_sampleid")
allinfo$NMF <- paste("NMF", allinfo$NMF, sep = "")


rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

exp_data <- rnamtr
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(exp_data["EGFR", ]),
                     FGFR3Exp = as.numeric(exp_data["FGFR3", ]),
                     N4Exp = as.numeric(exp_data["NECTIN4", ]))
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "trunc_anonymized_sample_ids")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
#tp_color <- c("Post" = "#F5A889", "Pre" = "#ACD6EC")
res_color <- c("Desert" = "#5C9E74", "Excluded" = "#A8DAB5", "Inflamed" = "#C973A0", "Unknown" = "white")
#sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("NMF1" = "#5C9E74", "NMF2" = "#A8DAB5", "NMF3" = "#F1C7DA", 
                "NMF4" = "#C973A0")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  NMF = refdat$NMF,
  CD8 = refdat$CD8_T_CELL_INFILTRATION,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    NMF = tcga_color,
    CD8 = res_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(0, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(0, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.IM010.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat
tmprefdat$NMF <- factor(tmprefdat$NMF, c("NMF1", "NMF2", "NMF3", 
                                         "NMF4"))
my_comp <- list(c("NMF1", "NMF4"), c("NMF2", "NMF4"), c("NMF3", "NMF4"))
p1 <- ggplot(tmprefdat, aes(x = NMF, y = B2L, fill = NMF)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM010.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()

tmprefdat <- refdat[refdat$CD8_T_CELL_INFILTRATION != "Unknown", ]
tmprefdat$Stage <- factor(tmprefdat$CD8_T_CELL_INFILTRATION, c("Desert", "Excluded", "Inflamed"))
my_comp <- list(c("Desert", "Inflamed"), c("Excluded", "Inflamed"))
p1 <- ggplot(tmprefdat, aes(x = Stage, y = B2L, fill = Stage)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = res_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM010.immune.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.IM010.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2,c("NMF1", "NMF2", "NMF3", 
                                    "NMF4"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = tcga_color)+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()



###############IM130###############
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_go29294_wo30070_wo29636_FMI_binary_calls.csv", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$truncated_anonymized_sampleid
cnvmtr <- cnvmtr[, -c(2, 3)]
names(cnvmtr) <- gsub("FMI_", "", names(cnvmtr))
for (gene in colnames(cnvmtr)[-1]) {
  cnv_raw <- cnvmtr[, gene]
  cnv_status <- ifelse(is.na(cnv_raw), "NA",
                       ifelse(cnv_raw == 1, "GA", "WT"))
  cnvmtr[[gene]] <- cnv_status
}

head(cnvmtr)

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/wo30070_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)

allinfo <- merge(saminfo, cnvmtr, by.x = "trunc_anonymized_sample_ids", by.y = "truncated_anonymized_sampleid")
allinfo$NMF <- paste("NMF", allinfo$NMF, sep = "")


rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

exp_data <- rnamtr
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(exp_data["EGFR", ]),
                     FGFR3Exp = as.numeric(exp_data["FGFR3", ]),
                     N4Exp = as.numeric(exp_data["NECTIN4", ]))
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "trunc_anonymized_sample_ids")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
#tp_color <- c("Post" = "#F5A889", "Pre" = "#ACD6EC")
res_color <- c("Desert" = "#5C9E74", "Excluded" = "#A8DAB5", "Inflamed" = "#C973A0", "Unknown" = "white")
#sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("NMF1" = "#5C9E74", "NMF2" = "#A8DAB5", "NMF3" = "#F1C7DA", 
                "NMF4" = "#C973A0")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  NMF = refdat$NMF,
  CD8 = refdat$CD8_T_CELL_INFILTRATION,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    NMF = tcga_color,
    CD8 = res_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(0, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(0, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.IM130.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat
tmprefdat$NMF <- factor(tmprefdat$NMF, c("NMF1", "NMF2", "NMF3", 
                                         "NMF4"))
my_comp <- list(c("NMF1", "NMF4"), c("NMF2", "NMF4"), c("NMF3", "NMF4"))
p1 <- ggplot(tmprefdat, aes(x = NMF, y = B2L, fill = NMF)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM130.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()

tmprefdat <- refdat[refdat$CD8_T_CELL_INFILTRATION != "Unknown", ]
tmprefdat$Stage <- factor(tmprefdat$CD8_T_CELL_INFILTRATION, c("Desert", "Excluded", "Inflamed"))
my_comp <- list(c("Desert", "Inflamed"), c("Excluded", "Inflamed"))
p1 <- ggplot(tmprefdat, aes(x = Stage, y = B2L, fill = Stage)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = res_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM130.immune.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.IM130.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2,c("NMF1", "NMF2", "NMF3", 
                                    "NMF4"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = tcga_color)+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()



###############IM210###############
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_go29294_wo30070_wo29636_FMI_binary_calls.csv", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$truncated_anonymized_sampleid
cnvmtr <- cnvmtr[, -c(2, 3)]
names(cnvmtr) <- gsub("FMI_", "", names(cnvmtr))
for (gene in colnames(cnvmtr)[-1]) {
  cnv_raw <- cnvmtr[, gene]
  cnv_status <- ifelse(is.na(cnv_raw), "NA",
                       ifelse(cnv_raw == 1, "GA", "WT"))
  cnvmtr[[gene]] <- cnv_status
}

head(cnvmtr)

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)

allinfo <- merge(saminfo, cnvmtr, by.x = "trunc_anonymized_sample_ids", by.y = "truncated_anonymized_sampleid")
allinfo$NMF <- paste("NMF", allinfo$NMF, sep = "")


rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

exp_data <- rnamtr
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(exp_data["EGFR", ]),
                     FGFR3Exp = as.numeric(exp_data["FGFR3", ]),
                     N4Exp = as.numeric(exp_data["NECTIN4", ]))
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "trunc_anonymized_sample_ids")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
#tp_color <- c("Post" = "#F5A889", "Pre" = "#ACD6EC")
res_color <- c("Desert" = "#5C9E74", "Excluded" = "#A8DAB5", "Inflamed" = "#C973A0", "Unknown" = "white")
#sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("NMF1" = "#5C9E74", "NMF2" = "#A8DAB5", "NMF3" = "#F1C7DA", 
                "NMF4" = "#C973A0")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  NMF = refdat$NMF,
  CD8 = refdat$CD8_T_CELL_INFILTRATION,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    NMF = tcga_color,
    CD8 = res_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(0, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(0, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.IM210.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat
tmprefdat$NMF <- factor(tmprefdat$NMF, c("NMF1", "NMF2", "NMF3", 
                                         "NMF4"))
my_comp <- list(c("NMF1", "NMF4"), c("NMF2", "NMF4"), c("NMF3", "NMF4"))
p1 <- ggplot(tmprefdat, aes(x = NMF, y = B2L, fill = NMF)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM210.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()

tmprefdat <- refdat[refdat$CD8_T_CELL_INFILTRATION != "Unknown", ]
tmprefdat$Stage <- factor(tmprefdat$CD8_T_CELL_INFILTRATION, c("Desert", "Excluded", "Inflamed"))
my_comp <- list(c("Desert", "Inflamed"), c("Excluded", "Inflamed"))
p1 <- ggplot(tmprefdat, aes(x = Stage, y = B2L, fill = Stage)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = res_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM210.immune.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.IM210.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2,c("NMF1", "NMF2", "NMF3", 
                                    "NMF4"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = tcga_color)+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()



###############IM211###############
cnvmtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29293_go29294_wo30070_wo29636_FMI_binary_calls.csv", header = T, stringsAsFactors = F, data.table = F)
rownames(cnvmtr) <- cnvmtr$truncated_anonymized_sampleid
cnvmtr <- cnvmtr[, -c(2, 3)]
names(cnvmtr) <- gsub("FMI_", "", names(cnvmtr))
for (gene in colnames(cnvmtr)[-1]) {
  cnv_raw <- cnvmtr[, gene]
  cnv_status <- ifelse(is.na(cnv_raw), "NA",
                       ifelse(cnv_raw == 1, "GA", "WT"))
  cnvmtr[[gene]] <- cnv_status
}

head(cnvmtr)

saminfo <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/go29294_clinical_biomarker.csv",
                 header = T, data.table = F, stringsAsFactors = F)

allinfo <- merge(saminfo, cnvmtr, by.x = "trunc_anonymized_sample_ids", by.y = "truncated_anonymized_sampleid")
allinfo$NMF <- paste("NMF", allinfo$NMF, sep = "")


rnamtr <- fread("/rsrch5/scratch/genomic_med/kyu3/projects/13.BC_matrix/BC2800/tpm_go29293_go29294_wo30070_wo29636.genesymbol.txt",
                header = T, stringsAsFactors = F, data.table = F)
rownames(rnamtr) <- rnamtr$V1
rnamtr <- rnamtr[, -1]
rnamtr <- rnamtr[, names(rnamtr) %in% saminfo$trunc_anonymized_sample_ids, ]

exp_data <- rnamtr
exp_data_centered <- t(apply(exp_data, 1, function(x) x - median(x, na.rm = TRUE)))

basal47 <- fread("/rsrch5/scratch/genomic_med/kyu3/datasets/01.Markers/basal47.txt", header = T, stringsAsFactors = F, data.table = F)
basal47 <- basal47[order(basal47$Basal), ]
rownames(basal47) <- basal47$id

base47_data <- exp_data_centered[rownames(exp_data_centered) %in% basal47$id, ]

# Match coefficients to expression data
coefficients <- basal47 %>% filter(id %in% rownames(base47_data))
base47_data <- base47_data[coefficients$id, ]

# Calculate scores
luminal_score <- colSums(base47_data * coefficients$Luminal)
basal_score <- colSums(base47_data * coefficients$Basal)

# Classify samples
classification <- ifelse(luminal_score > basal_score, "Luminal-like", "Basal-like")

refdat <- data.frame(Sample = names(exp_data),
                     Class = classification,
                     Basal = basal_score,
                     Luminal = luminal_score,
                     B2L = basal_score - luminal_score,
                     EGFRExp = as.numeric(exp_data["EGFR", ]),
                     FGFR3Exp = as.numeric(exp_data["FGFR3", ]),
                     N4Exp = as.numeric(exp_data["NECTIN4", ]))
refdat <- merge(refdat, allinfo, by.x = "Sample", by.y = "trunc_anonymized_sample_ids")
refdat$N4Class <- ifelse(refdat$N4Exp > median(refdat$N4Exp), "High", "Low")

refdat <- refdat[order(refdat$Basal), ]
refdat$EGFRExp <- scale(refdat$EGFRExp)
refdat$FGFR3Exp <- scale(refdat$FGFR3Exp)
refdat$EGFRExp[refdat$EGFRExp > 2] <- 2
refdat$EGFRExp[refdat$EGFRExp < -2] <- -2
refdat$FGFR3Exp[refdat$FGFR3Exp > 2] <- 2
refdat$FGFR3Exp[refdat$FGFR3Exp < -2] <- -2


library(ComplexHeatmap)
library(circlize)

class_colors <- c("Luminal-like" = "#148F28", "Basal-like" = "#EA71AE")

#threegrp <- c("#91C299", "#F5A889", "#ACD6EC")
#grpcolors <- c("#C05050", "#80A0C5", "#C09050", "#C080A0", "#A070A5")
#grpcolor2 <- c("#984EA3", "#F4B75B", "#1A65B4", "#148F28")

GAcolor <- c("GA" = "black", "WT" = "#E5E5E5", "NA" = "white")
#tp_color <- c("Post" = "#F5A889", "Pre" = "#ACD6EC")
res_color <- c("Desert" = "#5C9E74", "Excluded" = "#A8DAB5", "Inflamed" = "#C973A0", "Unknown" = "white")
#sta_color <- c("Stage I" = "#E5E5E5", "Stage II" = "#91C299", "Stage III" = "#ACD6EC", "Stage IV" = "#F5A889")
tcga_color <- c("NMF1" = "#5C9E74", "NMF2" = "#A8DAB5", "NMF3" = "#F1C7DA", 
                "NMF4" = "#C973A0")

top_anno <- HeatmapAnnotation(
  Class = refdat$Class,
  B2L = refdat$B2L,
  NMF = refdat$NMF,
  CD8 = refdat$CD8_T_CELL_INFILTRATION,
  NF2 = refdat$NF2,
  TP53 = refdat$TP53,
  RB1 = refdat$RB1,
  EGFR = refdat$EGFR,
  FGFR3 = refdat$FGFR3,
  EGFRExp = refdat$EGFRExp,
  FGFR3Exp = refdat$FGFR3Exp,
  col = list(
    Class = class_colors,
    NMF = tcga_color,
    CD8 = res_color,
    NF2 = GAcolor,
    TP53 = GAcolor,
    RB1 = GAcolor,
    EGFR = GAcolor,
    FGFR3 = GAcolor,
    B2L = colorRamp2(c(min(refdat$Basal), median(refdat$Basal), max(refdat$Basal)), c("#148F28", "#E5E5E5", "#EA71AE")),
    EGFRExp = colorRamp2(c(0, 2), c("#E5E5E5", "#FF9900")),
    FGFR3Exp = colorRamp2(c(0, 2), c("#E5E5E5", "#7D1EE5"))
  ),
  annotation_name_side = "left",
  annotation_height = unit(c(4, 4, 4, 4, 4), "mm")
)


empty_matrix <- matrix(0, nrow = 1, ncol = nrow(refdat))
colnames(empty_matrix) <- refdat$Sample


pdf(file.path(subres, "V2503.figure5.IM211.heatmap.pdf"), width = 7, height = 6)
print(Heatmap(empty_matrix,
              name = NULL,
              top_annotation = top_anno,
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_columns = FALSE,
              cluster_rows = FALSE))
dev.off()

tmprefdat <- refdat
tmprefdat$NMF <- factor(tmprefdat$NMF, c("NMF1", "NMF2", "NMF3", 
                                         "NMF4"))
my_comp <- list(c("NMF1", "NMF4"), c("NMF2", "NMF4"), c("NMF3", "NMF4"))
p1 <- ggplot(tmprefdat, aes(x = NMF, y = B2L, fill = NMF)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = tcga_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM211.subtype.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()

tmprefdat <- refdat[refdat$CD8_T_CELL_INFILTRATION != "Unknown", ]
tmprefdat$Stage <- factor(tmprefdat$CD8_T_CELL_INFILTRATION, c("Desert", "Excluded", "Inflamed"))
my_comp <- list(c("Desert", "Inflamed"), c("Excluded", "Inflamed"))
p1 <- ggplot(tmprefdat, aes(x = Stage, y = B2L, fill = Stage)) +
  #geom_violin(scale = "width", adjust = 1) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "B2L score") +
  stat_compare_means(comparisons = my_comp)+
  scale_fill_manual(values = res_color)+
  scale_y_continuous(expand = expansion(mult = c(0.1)), labels = function(x) sprintf("%.2f", x))+
  theme(legend.position = "none")

pdf(file.path(subres, "V2503.figure5.IM211.immune.pdf"), width = 3.5, height = 4)
print(p1)
dev.off()


tmprefdat <- refdat
library(ggalluvial)
pdf(file.path(subres, "V2503.figure5.IM211.N4.pdf"), width = 7, height = 3)
tmprefdat$Stage_Class <- tmprefdat[, "NMF"]
# Perform the chi-square test
chisq_test <- chisq.test(table(tmprefdat$N4Class, tmprefdat$Stage_Class))

#finalres$Stage_Class <- factor(finalres$Stage_Class, c("Stage IV", "Stage II-III"))
tmpres <- as.data.frame(prop.table(table(tmprefdat$N4Class, tmprefdat$Stage_Class), margin = 1))
# Stacked + percent
tmpres$Var2 <- factor(tmpres$Var2,c("NMF1", "NMF2", "NMF3", 
                                    "NMF4"))
# Stacked + percent
p <- ggplot(tmpres, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_flow(aes(alluvium = Var2, stratum = Var2), alpha= .8, 
            width = .6, knot.pos = 0.5) +
  geom_col(aes(fill=Var2), width = .6, color = "black")+
  scale_fill_manual(values = tcga_color)+
  xlab(paste("X-squared = ", round(chisq_test$statistic, 3), "\np-value = ", round(chisq_test$p.value, 3), sep = ""))+
  theme_classic()+ coord_flip()

print(p)
dev.off()





