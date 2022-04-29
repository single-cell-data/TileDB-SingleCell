library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/TileDB-SingleCell/pbmc_gen/_data/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

saveRDS(pbmc, file = "/TileDB-SingleCell/pbmc_gen/_data/pbmc_test.rds")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2)

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.ids)

saveRDS(pbmc, file = "/TileDB-SingleCell/pbmc_gen/_data/pbmc.rds")

SaveH5Seurat(pbmc, filename = "/TileDB-SingleCell/pbmc_gen/_data/pbmc.h5Seurat")

Convert("/TileDB-SingleCell/pbmc_gen/_data/pbmc.h5Seurat", dest = "h5ad")





