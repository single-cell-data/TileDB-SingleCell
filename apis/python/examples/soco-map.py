#!/usr/bin/env python -i

import tiledb
import tiledbsc
import tiledbsc
import tiledbsc.io as io
import sys, os, time

import anndata
import anndata as ad  # so we can type it either way
import pandas
import pandas as pd  # so we can type it either way
import numpy
import numpy as np  # so we can type it either way
import scipy

from typing import List, Dict

if len(sys.argv) == 1:
    soco_path = "soma-collection"
elif len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsc.SOMACollection(
    soco_path, ctx=tiledb.Ctx({"py.init_buffer_bytes": 4 * 1024**3})
)

# ----------------------------------------------------------------
# obs_ids_per_soma = soco.map(
#     lambda soma: soma.obs.ids_from_attribute_filter(
#         query_string='cell_type == "pericyte cell"', attrs=["cell_type"]
#     )
# )

# lens = {name: len(obs_ids) for name, obs_ids in obs_ids_per_soma.items()}
# for name, length in lens.items():
#   print(length, name)
# 0 esophagus-epithelium
# 0 acute-covid19-cohort
# 0 healthy-human-liver-integrated
# 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
# 0 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
# 0 longitudinal-profiling-49
# 0 autoimmunity-pbmcs
# 0 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
# 121 ileum
# 0 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
# 0 human-kidney-tumors-wilms
# 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
# 0 retinal-ganglion-cells-in-human-retina
# 14 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
# 0 49-years-old-male-airway-wash-3-days-post-intubation
# 0 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
# 0 4056cbab-2a32-4c9e-a55f-c930bc793fb6
# 0 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562

# X_dfs = soco.map(lambda soma, obs_ids: soma.X.data.df(obs_ids=obs_ids), obs_ids_per_soma)

# ----------------------------------------------------------------
# print("TWO-SIDED QUERY")
# soco_attribute_filter_and_store(
#     soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
#     output_h5ad_path="mini-atlas-two-sided.h5ad",
#     output_soma_path="mini-atlas-two-sided",
#     obs_attr_names=["cell_type"],
#     obs_query_string='cell_type == "B cell"',
#     var_attr_names=["feature_name"],
#     var_query_string='feature_name == "MT-CO3"',
# )

# obs_ids_per_soma = soco.map(
#     lambda soma: soma.obs.ids_from_attribute_filter(
#         query_string='cell_type == "B cell"', attrs=["cell_type"]
#     )
# )
# print()
# for name, obs_ids in obs_ids_per_soma.items():
#     print(len(obs_ids), name)
# var_ids_per_soma = soco.map(
#     lambda soma: soma.var.ids_from_attribute_filter(
#         query_string='feature_name == "MT-CO3"', attrs=["feature_name"]
#     )
# )
# print()
# for name, var_ids in var_ids_per_soma.items():
#     print(len(var_ids), name)

# ovids_per_soma = {}
# for name in obs_ids_per_soma.keys():
#     ovids_per_soma[name] = {
#         "obs_ids": obs_ids_per_soma[name],
#         "var_ids": var_ids_per_soma[name],
#     }
# X_dfs = soco.map(
#     lambda soma, ovids: soma.X.data.df(
#         obs_ids=ovids["obs_ids"], var_ids=ovids["var_ids"]
#     ),
#     ovids_per_soma,
# )
# print()
# for name, X_df in X_dfs.items():
#     dense = tiledbsc.util.triples_to_dense_df(X_df)
#     print(X_df.shape, dense.shape, name)

# ----------------------------------------------------------------
# bruce:

# obs_id_map = soco.filter(dim='obs', 'cell_type == "foobar"')   # returns { soma_name: [obs_id, ...] }
# obs_df_slices = []
# for soma_name, obs_ids in obs_id_map.items():
#     df = soco[soma_name].obs.df[obs_ids]   # returns slice of obs as pandas dataframe
#     obs_df_slices.append(df)
# full_result = pandas.concat(obs_df_slices, join="outer")

obs_ids_per_soma = soco.map(
    lambda soma: soma.obs.ids_from_attribute_filter(
        query_string='cell_type == "B cell"', attrs=["cell_type"]
    )
)
print()
for name, obs_ids in obs_ids_per_soma.items():
    print(len(obs_ids), name)

obs_df_slices = []
print()
for soma_name, obs_ids in obs_ids_per_soma.items():
    print("...", len(obs_ids), soma_name)
    df = soco[soma_name].obs.df(ids=obs_ids)  # returns slice of obs as pandas dataframe
    obs_df_slices.append(df)
full_result = pandas.concat(obs_df_slices, join="outer")
print(full_result)
