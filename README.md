# Cell Painting morphological profiling for mammary carcinogenicity risk assessment

Code and results underlying our study on using Cell Painting (JUMP-CP) morphological
profiles to assess mammary carcinogenicity risk, through descriptive analysis (UMAP,
hierarchical clustering) and a Guilt-By-Association (GBA) similarity-based scoring
framework applied prospectively to EDLists endocrine-disruptor candidates.

If you use this code or these results, please cite:

> Achebouche R., Taboureau O. *[title]*. [journal], [year]. doi:[XXXX]

Archived at Zenodo: doi: 10.5281/zenodo.21512685

---

## Repository structure

```
.
├── data/
│   ├── chemicals/
│   │   ├── KAY_RUDEL_2024_CombinedSuppTables_rev2.xlsx   # source: Kay & Rudel (2024) MC/Non-MC compendium
│   │   ├── prepare_chemical_lists.py                     # -> BreastCancerChemList_{MC,NOMC}_v2.csv
│   │   ├── BreastCancerChemList_MC_v2.csv                 # 279 mammary carcinogens
│   │   ├── BreastCancerChemList_NOMC_v2.csv               # 850 non-mammary-carcinogens
│   │   └── ED_Lists/
│   │       ├── ED_list1.xlsx, ED_list2.xlsx, ED_list3.xlsx  # source: EDLists.org (List I/II/III)
│   │       ├── prepare_data.py                            # -> ED_lists_combined.csv (PubChemPy lookup)
│   │       └── ED_lists_combined.csv
│   ├── metadata/
│   │   ├── comp_metadata_jumpcp.csv                       # JUMP-CP compound metadata (JCP2022 <-> INCHIKEY)
│   │   └── jumpcp_int_feature_names.csv                   # harmony_N -> real Cell Painting feature name
│   ├── jump-preprocessed/                                 # raw JUMP-CP parquet (not tracked - see Data availability)
│   ├── preprocess_JUMPCP.py                                # -> chemlist_mc_nomc.csv + data/processed/int_ns_h/*.csv
│   ├── chemlist_mc_nomc.csv                                 # combined MC + Non-MC reference list (1,105 compounds)
│   ├── processed/int_ns_h/                                  # per-compound, per-source median profiles (748 features)
│   │   ├── mc_jump_profiles.csv, nomc_jump_profiles.csv, nomc_genotox_jump_profiles.csv
│   │   ├── dmso_jump_profiles.csv
│   │   └── edlist1_jump_profiles.csv, edlist2_jump_profiles.csv
│   └── results_tables/                                      # analysis outputs (see below)
├── analyses/
│   ├── UMAP/umap_exploratory_analysis.ipynb                 # Methods 2.3.1, Figure 1
│   ├── ClusterMap/clustermap.ipynb                           # Methods 2.3.2, Figure 2
│   └── Guilt-By-Association/guilt-by-association_2_clean.ipynb  # Methods 2.4, Results 3.2, Table 1
├── figures/                                                   # all figures, referenced by the notebooks above
├── requirements.txt
├── LICENSE
└── README.md
```

## Pipeline / reproduction order

1. **`data/chemicals/prepare_chemical_lists.py`** - builds the primary MC / Non-MC chemical
   list from the Kay & Rudel (2024) compendium, resolving INCHIKEY/SMILES via PubChemPy.
   (Already run; outputs are checked in. Re-running requires network access to PubChem and
   may return slightly different matches if PubChem's data has changed.)
2. **`data/chemicals/ED_Lists/prepare_data.py`** - concatenates the three EDLists.org
   candidate lists and resolves structural identifiers via PubChemPy. Same caveat as above.
3. **`data/preprocess_JUMPCP.py`** - extracts per-compound, per-source median JUMP-CP
   profiles for each chemical list from the raw Harmony-corrected parquet file. Requires
   `data/jump-preprocessed/profiles_wellpos_cc_var_mad_int_featselect_harmony.parquet`
   (see *Data availability* below) and `data/metadata/comp_metadata_jumpcp.csv` (included).
   Run with `python3 data/preprocess_JUMPCP.py`.
4. **`analyses/UMAP/umap_exploratory_analysis.ipynb`** - UMAP projection of MC / non-genotoxic
   Non-MC / DMSO profiles (Figure 1).
5. **`analyses/ClusterMap/clustermap.ipynb`** - hierarchical clustering heatmap of the same
   compound set, colored by category and by hormonal activity (Figure 2), plus a
   genotoxicity-colored supplementary variant.
6. **`analyses/Guilt-By-Association/guilt-by-association_2_clean.ipynb`** - the full GBA
   pipeline: PCA / Random-Forest feature selection, dynamic-tree-cut clustering, five
   similarity metrics x three risk-score aggregation strategies x three feature
   representations, leave-one-out cross-validation, and prospective risk scoring of the
   EDList I / II candidates (Table 1, Results 3.2). Includes a feature-selection stability
   diagnostic (Section 6.1) explaining why the RF-selected-feature scenario is evaluated
   with a nested (non-leaky) cross-validation.

Steps 1-3 only need to be rerun if you want to regenerate the derived CSVs from scratch;
their outputs are already included in this repository (except the large raw parquet).
Steps 4-6 can be run directly against the included `data/processed/int_ns_h/*.csv` files.

## Data availability

This repository includes every file the analysis notebooks read directly, **except** the
raw JUMP-CP parquet file consumed by `data/preprocess_JUMPCP.py`
(`data/jump-preprocessed/profiles_wellpos_cc_var_mad_int_featselect_harmony.parquet`,
~2.5 GB), which is excluded from version control for size reasons. Cell Painting
morphological profiles (cpg0016) are publicly available from the JUMP Cell Painting
Consortium: https://github.com/jump-cellpainting/datasets, processed with the
`jump-profiling-recipe` workflow (https://github.com/broadinstitute/jump-profiling-recipe)
using the well-position-correction / cell-count-variance-normalization / MAD /
inverse-normal-transform / variance-threshold-feature-selection / Harmony-batch-correction
configuration described in Methods 2.2. The per-compound profiles derived from it
(`data/processed/int_ns_h/*.csv`) **are** included, so the analysis notebooks (steps 4-6
above) can be run without downloading the raw parquet.

Mammary carcinogenicity annotations are from Kay & Rudel (2024),
https://doi.org/10.1289/ehp13233. Endocrine-disruptor candidate lists are from EDLists,
https://edlists.org.

## Environment

Python 3.13. Install dependencies with:

```bash
pip install -r requirements.txt
```

## License

MIT License - see [LICENSE](LICENSE).

## Authors

Rayane Achebouche, Olivier Taboureau - Université Paris Cité, CNRS, Unité de Biologie
Fonctionnelle et Adaptative, INSERM U1133, Paris, France.

This study was supported by the Fondation pour la Recherche Médicale (FRM), grant
ECO202306017371.
