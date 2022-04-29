
# Overview

* Sample data:
  * The raw barcodes, genes, and counts matrix was obtained via:
    * wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
    * tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
* Inspecting HDF5 input files
  * `./desc-ann.py ../../pbmc_gen/_data/pbmc.h5ad`
* Ingesting
  * `./ingestor.py ./anndata/pbmc3k_processed.h5ad`
  * Output is in `_data/tiledb-data/pbmc_processed`
* Inspecting TileDB output groups
  * `./desc-tiledb.py ../../_data/tiledb-data/pbmc_processed`

# Build and Activate Conda Env
```
* Note: The first line of the yml file sets the new environment's name *

conda env create -f TileDB-SingleCell/pbmc_gen/singleCellDataEnv.yml

conda activate scapi-test
```

# Generate .h5ad and .h5Seurat Datasets
```
Rscript TileDB-SingleCell/pbmc_gen/_scripts/_r_dev/generate_pbmcData.R
```

# Ingestion Process - R
```
Rscript pbmc_gen/_scripts/_r_dev/single_cell_api_ingestion.R

* Output location: _data/

* Inspecting TileDB output groups
  ./desc-tiledb.py ../../pbmc_gen/_data/sc_dataset
```

# Ingestion Process - Python
```
cd apis/python

* Inspecting HDF5 input files
  ./desc-ann.py ../../pbmc_gen/_data/pbmc.h5ad
  
* Ingesting
  ./ingestor.py ../../pbmc_gen/_data/pbmc.h5ad ../../pbmc_gen/_data/tiledb-data/pbmc_processed
  * Output is in `../../pbmc_gen/_data/tiledb-data/pbmc3k_processed`
  
* Inspecting TileDB output groups
  ./desc-tiledb.py ../../pbmc_gen/_data/tiledb-data/pbmc_processed
```