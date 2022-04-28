#!/usr/bin/env python

# ================================================================
# It looks as if we need to use numpy 1.21.5 or lower with this*
#
#
# ================================================================

import numpy as np
import pandas as pd
import pylab as pl
import scanpy as sc
import scanpy as sc
import pandas as pd
from matplotlib.pyplot import rc_context


# Preprocessing and clustering 3k PBMCs
# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
from scanpy.external import tl

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = '../data/pbmc3k.h5ad'  # the file that will store the analysis results

adata = sc.read_10x_mtx(
    '../data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

print(adata)

sc.pl.highest_expr_genes(adata, n_top=20, )

# Some basic filtering of the genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Remove cells that have too many mitochondrial genes expressed or too many total counts:
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# Actually do the filtering by slicing the AnnData object.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell, so that counts become
# comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data:
sc.pp.log1p(adata)

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

'''
Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later 
use in differential testing and visualizations of gene expression. This simply freezes the state of the 
AnnData object.
*You can get back an AnnData of the object in .raw by calling .raw.to_adata().
'''
adata.raw = adata

# Do the filtering
adata = adata[:, adata.var.highly_variable]

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the
# data to unit variance.
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

# Scatter plot of the PCA coordinates, for comparission purposes of the R and Python
sc.pl.pca(adata, color='CST3')

sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)
print(adata)


# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix.
# You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take
# the following values.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Embedding the neighborhood graph
# We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below.
# It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better
# preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity
# violations. They can usually be remedied by running:


# sc.tl.tsne(adata)
sc.tl.louvain(adata)

sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')

sc.tl.umap(adata)

sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])

# As we set the .raw attribute of adata, the previous plots showed the “raw”
# (normalized, logarithmized, but uncorrected) gene expression. You can also plot the scaled
# and corrected gene expression by explicitly stating that you don’t want to use .raw.
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)

'''
Clustering the neighborhood graph
As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method 
(community detection based on optimizing modularity) by Traag *et al.* (2018). Note that Leiden clustering 
directly clusters the neighborhood graph of cells, which we already computed in the previous section.
'''
sc.tl.leiden(adata)

# Plot the clusters, which agree quite well with the result of Seurat.
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])

# Save the result.
adata.write(results_file)

# Finding marker genes
# Let us compute a ranking for the highly differential genes in each cluster. For this, by default,
# the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method
# to do so is the t-test.

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

sc.settings.verbosity = 2  # reduce the verbosity

# The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar. We recommend using the latter
# in publications, see e.g., Sonison & Robinson (2018). You might also consider much more powerful differential
# testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Save the result
adata.write(results_file)

# As an alternative, let us rank genes using logistic regression. For instance, this has been suggested by
# Natranos et al. (2018). The essential difference is that here, we use a multi-variate appraoch whereas
# conventional differential tests are uni-variate. Clark et al. (2014) has more details.
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Define a list of marker genes for later reference.
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']

# Reload the object that has been save with the Wilcoxon Rank-Sum test result.
adata = sc.read(results_file)
adata.uns['log1p']["base"] = None

# Show the 10 top ranked genes per cluster 0, 1, …, 7 in a dataframe.
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

# Get a table with the scores and groups.
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'pvals']}).head(5)

# Compare to a single cluster:
sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)

# If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin.
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

# Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups):
adata = sc.read(results_file)
adata.uns['log1p']["base"] = None  # Temp fix to get around the error: if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
                                   # KeyError: 'base'
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

# If you want to compare a certain gene across groups, use the following.
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')

# Actually mark the cell types.
new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')

# Now that we annotated the cell types, let us visualize the marker genes.
sc.pl.dotplot(adata, marker_genes, groupby='leiden');

# There is also a very compact violin plot.
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90);


# During the course of this analysis, the AnnData accumlated the following annotations.
print(adata)

adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing
                                               # and subsequent reading

#adata.raw.to_adata().write('./data/pbmc3k_withoutX.h5ad')


