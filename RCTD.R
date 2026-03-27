
library(spacexr)
library(Seurat)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
sample = args[1]

#sample = "BC01"

dir.create(paste("../9.spaAnaly/NewCombined.Analy/01.RCTD/", sample, sep = ""))

outdir <- paste("../9.spaAnaly/NewCombined.Analy/01.RCTD/", sample, sep = "")

###########################################RCTD###########################################
resdir <- "../9.spaAnaly/NewCombined.Analy/SeuratObj"
spadat <- readRDS(file.path(resdir, paste(sample, ".seurat.umap.spatial.rds", sep = "")))

reference <- readRDS("../../11.JJ_Landscape_BC/4.Results/combined.tme.epi.rds")  # Provide your reference dataset path
ref <- UpdateSeuratObject(reference)
Idents(ref) <- "ModCluster"

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$ModCluster)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)


counts <- spadat[["Spatial"]]$counts
coords <- GetTissueCoordinates(spadat)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")

spadat <- AddMetaData(spadat, metadata = RCTD@results$weights)
# Extract results

# Visualize cell type proportions
# Plot example for specific cell types; you can expand this to show spatial distribution
pdf(file.path(outdir, "cell.fraction.visu.pdf"), width = 12, height = 12)
SpatialFeaturePlot(spadat, features = colnames(RCTD@results$weights), ncol = 4, image.alpha = 0, pt.size.factor = 6)
dev.off()

write.table(RCTD@results$weights, file.path(outdir, "cell.fraction.full.txt"),
            row.names = T, col.names = NA, sep = "\t", quote = F)


