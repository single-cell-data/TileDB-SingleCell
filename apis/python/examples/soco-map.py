#!/usr/bin/env python -i

import sys

import pandas
import tiledb

import tiledbsc

# ================================================================
# SOCO for most of these examples: 18 SOMAs, total 9.2GB
# ================================================================

# ----------------------------------------------------------------
# Use this for full-SOMA, all-at-once queries:
ctx = tiledb.Ctx({"py.init_buffer_bytes": 4 * 1024**3})

# ================================================================
# SOCO for most of these examples: 18 SOMAs, total 9.2GB
# ================================================================

# ----------------------------------------------------------------
# Use this for full-SOMA, all-at-once queries:
ctx = tiledb.Ctx({"py.init_buffer_bytes": 4 * 1024**3})

if len(sys.argv) == 1:
    soco_path = "soma-collection"
elif len(sys.argv) == 2:
    soco_path = sys.argv[1]
elif len(sys.argv) == 3:  # anything after argv[2], like '.'
    soco_path = sys.argv[1]
    # Use this to test return_incomplete=True:
    ctx = tiledb.Ctx({})
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsc.SOMACollection(soco_path, ctx=ctx)


# ----------------------------------------------------------------
def show_cell_type_counts(soco):
    obs_ids_per_soma = soco.map(
        lambda soma: soma.obs.ids_from_query(
            query_string='cell_type == "pericyte cell"', attrs=["cell_type"]
        )
    )

    lens = {name: len(obs_ids) for name, obs_ids in obs_ids_per_soma.items()}
    for name, length in lens.items():
        print(length, name)

    # >>> t1=time.time(); show_cell_type_counts(soco); t2=time.time(); print(); print(t2-t1, 'seconds')

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

    # 0.306948184967041 seconds


# ----------------------------------------------------------------
def show_X_dfs(soco):
    obs_ids_per_soma = soco.map(
        lambda soma: soma.obs.ids_from_query(
            query_string='cell_type == "pericyte cell"', attrs=["cell_type"]
        )
    )

    X_dfs = soco.map(
        lambda soma, obs_ids: soma.X.data.df(obs_ids=obs_ids), obs_ids_per_soma
    )

    # print(X_dfs)
    for name, X_df in X_dfs.items():
        print(name, X_df.shape)

    # >>> t1=time.time(); show_X_dfs(soco); t2=time.time(); print(); print(t2-t1, 'seconds')
    # (smaller SOCO here)

    # subset-soma-01 (248619, 1)
    # subset-soma-04 (225615, 1)
    # subset-soma-03 (248130, 1)
    # subset-soma-02 (242145, 1)

    # 0.42031097412109375 seconds


# ----------------------------------------------------------------
def bruce_example(soco):
    # bruce:

    # obs_id_map = soco.filter(dim='obs', 'cell_type == "foobar"')   # returns { soma_name: [obs_id, ...] }
    # obs_df_slices = []
    # for soma_name, obs_ids in obs_id_map.items():
    #     df = soco[soma_name].obs.df[obs_ids]   # returns slice of obs as pandas dataframe
    #     obs_df_slices.append(df)
    # full_result = pandas.concat(obs_df_slices, join="outer")

    obs_ids_per_soma = soco.map(
        lambda soma: soma.obs.ids_from_query(
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
        df = soco[soma_name].obs.df(
            ids=obs_ids
        )  # returns slice of obs as pandas dataframe
        obs_df_slices.append(df)
    full_result = pandas.concat(obs_df_slices, join="outer")
    print(full_result)

    # >>> t1=time.time(); bruce_example(soco); t2=time.time(); print(); print(t2-t1, 'seconds')

    # 773 esophagus-epithelium
    # 6131 acute-covid19-cohort
    # 1780 healthy-human-liver-integrated
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # 0 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # 1453 longitudinal-profiling-49
    # 510 autoimmunity-pbmcs
    # 0 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # 3183 ileum
    # 124 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # 0 human-kidney-tumors-wilms
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # 0 retinal-ganglion-cells-in-human-retina
    # 0 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # 14 49-years-old-male-airway-wash-3-days-post-intubation
    # 0 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # 306 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # 0 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562
    #
    # ... 773 esophagus-epithelium
    # ... 6131 acute-covid19-cohort
    # ... 1780 healthy-human-liver-integrated
    # ... 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # ... 0 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # ... 1453 longitudinal-profiling-49
    # ... 510 autoimmunity-pbmcs
    # ... 0 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # ... 3183 ileum
    # ... 124 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # ... 0 human-kidney-tumors-wilms
    # ... 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # ... 0 retinal-ganglion-cells-in-human-retina
    # ... 0 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # ... 14 49-years-old-male-airway-wash-3-days-post-intubation
    # ... 0 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # ... 306 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # ... 0 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562
    #                                           assay_ontology_term_id cell_type_ontology_term_id  ...     sex                        tissue
    # obs_id                                                                                       ...
    # AAACCTGAGTCGTTTG-1-HCATisStab7646031                 EFO:0009899                 CL:0000236  ...    male       epithelium of esophagus
    # AAACCTGTCTGTGCAA-1-HCATisStabAug177184862            EFO:0009899                 CL:0000236  ...  female       epithelium of esophagus
    # AAACGGGAGACGACGT-1-HCATisStabAug177376564            EFO:0009899                 CL:0000236  ...  female       epithelium of esophagus
    # AAACGGGGTCTCTTAT-1-HCATisStab7646031                 EFO:0009899                 CL:0000236  ...    male       epithelium of esophagus
    # AAACGGGGTGTCAATC-1-HCATisStabAug177184862            EFO:0009899                 CL:0000236  ...  female       epithelium of esophagus
    # ...                                                          ...                        ...  ...     ...                           ...
    # H12_F10_RT_BC_99_Lig_BC_172                          EFO:0010550                 CL:0000010  ...  female  cultured cell (cell culture)
    # H12_F10_RT_BC_99_Lig_BC_226                          EFO:0010550                 CL:0000010  ...  female  cultured cell (cell culture)
    # H12_F10_RT_BC_99_Lig_BC_337                          EFO:0010550                 CL:0000010  ...  female  cultured cell (cell culture)
    # H12_F10_RT_BC_9_Lig_BC_21                            EFO:0010550                 CL:0000010  ...  female  cultured cell (cell culture)
    # H12_F10_RT_BC_9_Lig_BC_274                           EFO:0010550                 CL:0000010  ...  female  cultured cell (cell culture)
    #
    # [435398 rows x 16 columns]

    # 5.637082099914551 seconds


# ----------------------------------------------------------------
def two_sided_query(soco):
    print("TWO-SIDED QUERY")

    obs_ids_per_soma = soco.map(
        lambda soma: soma.obs.ids_from_query(
            query_string='cell_type == "B cell"', attrs=["cell_type"]
        )
    )
    print()
    for name, obs_ids in obs_ids_per_soma.items():
        print(len(obs_ids), name)
    var_ids_per_soma = soco.map(
        lambda soma: soma.var.ids_from_query(
            query_string='feature_name == "MT-CO3"', attrs=["feature_name"]
        )
    )
    print()
    for name, var_ids in var_ids_per_soma.items():
        print(len(var_ids), name)

    ovids_per_soma = {}
    for name in obs_ids_per_soma.keys():
        ovids_per_soma[name] = {
            "obs_ids": obs_ids_per_soma[name],
            "var_ids": var_ids_per_soma[name],
        }
    X_dfs = soco.map(
        lambda soma, ovids: soma.X.data.df(
            obs_ids=ovids["obs_ids"], var_ids=ovids["var_ids"]
        ),
        ovids_per_soma,
    )
    print()
    for name, X_df in X_dfs.items():
        dense = tiledbsc.util.triples_to_dense_df(X_df)
        print(X_df.shape, dense.shape, name)

    # >>> t1=time.time(); two_sided_query(soco); t2=time.time(); print(); print(t2-t1, 'seconds')

    # TWO-SIDED QUERY
    #
    # 773 esophagus-epithelium
    # 6131 acute-covid19-cohort
    # 1780 healthy-human-liver-integrated
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # 0 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # 1453 longitudinal-profiling-49
    # 510 autoimmunity-pbmcs
    # 0 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # 3183 ileum
    # 124 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # 0 human-kidney-tumors-wilms
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # 0 retinal-ganglion-cells-in-human-retina
    # 0 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # 14 49-years-old-male-airway-wash-3-days-post-intubation
    # 0 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # 306 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # 0 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562
    #
    # 1 esophagus-epithelium
    # 1 acute-covid19-cohort
    # 1 healthy-human-liver-integrated
    # 1 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # 1 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # 1 longitudinal-profiling-49
    # 1 autoimmunity-pbmcs
    # 1 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # 1 ileum
    # 1 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # 1 human-kidney-tumors-wilms
    # 1 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # 1 retinal-ganglion-cells-in-human-retina
    # 1 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # 1 49-years-old-male-airway-wash-3-days-post-intubation
    # 1 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # 1 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # 1 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562
    #
    # (754, 1) (754, 1) esophagus-epithelium
    # (5795, 1) (5795, 1) acute-covid19-cohort
    # (1742, 1) (1742, 1) healthy-human-liver-integrated
    # (40105, 1) (40105, 1) molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # (10564, 1) (10564, 1) scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # (1453, 1) (1453, 1) longitudinal-profiling-49
    # (510, 1) (510, 1) autoimmunity-pbmcs
    # (48954, 1) (48954, 1) individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # (3166, 1) (3166, 1) ileum
    # (123, 1) (123, 1) 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # (4344, 1) (4344, 1) human-kidney-tumors-wilms
    # (23669, 1) (23669, 1) molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # (8778, 1) (8778, 1) retinal-ganglion-cells-in-human-retina
    # (643, 1) (643, 1) spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # (14, 1) (14, 1) 49-years-old-male-airway-wash-3-days-post-intubation
    # (88448, 1) (88448, 1) individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # (305, 1) (305, 1) 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # (112907, 1) (112907, 1) massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562

    # 23.38720679283142 seconds


# ----------------------------------------------------------------
def paginated_query(soco):
    print("PAGINATED QUERY")

    obs_ids_per_soma = soco.map(
        lambda soma: soma.obs.ids_from_query(
            query_string='cell_type == "B cell"', attrs=["cell_type"]
        )
    )
    print()
    for name, obs_ids in obs_ids_per_soma.items():
        print(len(obs_ids), name)
    print()

    for soma in soco:
        print(soma.name + ":")
        with soma.X.data._open() as A:
            obs_ids = obs_ids_per_soma[soma.name]
            iterable = A.query(return_incomplete=True).df[obs_ids, :]

            for result in iterable:
                print("  ", result.shape)

            # TODO:
            # o Work this into the Python SOMA API
            # o Come up with aggregators that can use ("reduce") the df pieces that come back

    # >>> t1=time.time(); paginated_query(soco); t2=time.time(); print(); print(t2-t1, 'seconds')

    # PAGINATED QUERY
    #
    # 773 esophagus-epithelium
    # 6131 acute-covid19-cohort
    # 1780 healthy-human-liver-integrated
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus
    # 0 scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-
    # 1453 longitudinal-profiling-49
    # 510 autoimmunity-pbmcs
    # 0 individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al
    # 3183 ileum
    # 124 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e
    # 0 human-kidney-tumors-wilms
    # 0 molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex
    # 0 retinal-ganglion-cells-in-human-retina
    # 0 spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9
    # 14 49-years-old-male-airway-wash-3-days-post-intubation
    # 0 individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al
    # 306 4056cbab-2a32-4c9e-a55f-c930bc793fb6
    # 0 massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562
    # esophagus-epithelium:
    #    (516082, 3)
    # acute-covid19-cohort:
    #    (7629480, 3)
    # healthy-human-liver-integrated:
    #    (1189534, 3)
    # molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimers-disease-superior-frontal-gyrus:
    #    (91548313, 3)
    # scrna-seq-data-analysis-of-endothelium-enriched-mesenteric-arterial-tissues-from-human-donors-:
    #    (23846766, 3)
    # longitudinal-profiling-49:
    #    (2198966, 3)
    # autoimmunity-pbmcs:
    #    (816874, 3)
    # individual-single-cell-rna-seq-pbmc-data-from-arunachalam-et-al:
    #    (103908557, 3)
    # ileum:
    #    (2159794, 3)
    # 0cfab2d4-1b79-444e-8cbe-2ca9671ca85e:
    #    (301876, 3)
    # human-kidney-tumors-wilms:
    #    (14120269, 3)
    # molecular-characterization-of-selectively-vulnerable-neurons-in-alzheimer’s-disease-entorhinal-cortex:
    #    (45384499, 3)
    # retinal-ganglion-cells-in-human-retina:
    #    (52490106, 3)
    # spatiotemporal-analysis-of-human-intestinal-development-at-single-cell-resolution-fetal-a9:
    #    (2433947, 3)
    # 49-years-old-male-airway-wash-3-days-post-intubation:
    #    (29392, 3)
    # individual-single-cell-rna-seq-pbmc-data-from-schulte-schrepping-et-al:
    #    (147904532, 3)
    # 4056cbab-2a32-4c9e-a55f-c930bc793fb6:
    #    (615667, 3)
    # massively-multiplex-chemical-transcriptomics-at-single-cell-resolution--k562:
    #    (124254874, 3)

    # 117.86341905593872 seconds
