library(tiledbsc)
library(tiledb)
library(SeuratObject)
library(SeuratDisk)

## Load data
pbmc3k <- LoadH5Seurat("/TileDB-SingleCell/pbmc_gen/_data/pbmc.h5Seurat")

## Convert a `Seurat` object to a TileDB-backed `sc_dataset`
scdataset <- SCDataset$new(uri = file.path("/TileDB-SingleCell/pbmc_gen/_data/", "sc_dataset"))

scdataset$from_seurat(object = pbmc3k)

scdataset$to_seurat(project = "SCDataset-Test")

## Convert a Seurat `Assay` to TileDB-backed `sc_group`
scgroup <- SCGroup$new(uri = file.path("/TileDB-SingleCell/pbmc_gen/_data/", "sc_group"))

scgroup$from_seurat_assay(
  object = pbmc_small[["RNA"]],
  obs = pbmc_small[[]]
)

scgroup_obs <- scdataset$scgroups$RNA$obs

obs_array <- scgroup_obs$tiledb_array(
  return_as = "tibble",
  attrs = c("nCount_RNA", "nFeature_RNA"),
  query_condition = parse_query_condition(nFeature_RNA < 2500)
)


