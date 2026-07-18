"""
Build per-compound JUMP-CP morphological profiles for the MC / Non-MC and
EDLists chemical sets.

Starting from the batch-corrected JUMP-CP parquet file (well-position
correction, cell-count variance normalization, MAD scaling, inverse-normal
transform, variance-threshold feature selection, Harmony batch correction
- see Methods 2.2), this script extracts, for each chemical list, the
per-(compound, imaging source) median profile and saves it as a CSV under
data/processed/int_ns_h/.

Run with: conda run -n lincs python3 data/preprocess_JUMPCP.py
"""
import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).parent
PARQUET_PATH = DATA_DIR / "jump-preprocessed" / "profiles_wellpos_cc_var_mad_int_featselect_harmony.parquet"
METADATA_PATH = DATA_DIR / "metadata" / "comp_metadata_jumpcp.csv"
OUT_DIR = DATA_DIR / "processed" / "int_ns_h"

DMSO_INCHIKEY = "IAZDPXIOMUYVGZ-UHFFFAOYSA-N"


def load_chemical_data():
    """Load the MC / Non-MC reference set and the EDLists List I / II candidates."""
    mc_chem_df = pd.read_csv(DATA_DIR / "chemicals/BreastCancerChemList_MC_v2.csv").dropna(subset=["INCHIKEY"])
    nomc_chem_df = pd.read_csv(DATA_DIR / "chemicals/BreastCancerChemList_NOMC_v2.csv").dropna(subset=["INCHIKEY"])
    mc_chem_df = mc_chem_df.rename(columns={"MammaryTumorEvidence": "Category", "preferred_name": "Name"})
    nomc_chem_df = nomc_chem_df.rename(columns={"MammaryTumorEvidence": "Category", "preferred_name": "Name"})
    nomc_chem_df = nomc_chem_df.copy()
    nomc_chem_df["Category"] = "NOMC"
    nomc_chem_genotox_df = nomc_chem_df[nomc_chem_df["Genotoxicity"] == "negative"]

    # EDLists List I and List II candidates (List III is excluded from the
    # analysis - see Methods 2.1 - so only these two are processed here).
    edlists_dataset = pd.read_csv(DATA_DIR / "chemicals/ED_Lists/ED_lists_combined.csv").sort_values(
        by=["ID (for name)", "Year"]).drop_duplicates(
        subset=["InChIKey"], keep="last").dropna(subset=["InChIKey"])
    edlist1_chem_df = edlists_dataset[edlists_dataset["List category"] == "List I"].copy()
    edlist2_chem_df = edlists_dataset[edlists_dataset["List category"] == "List II"].copy()

    return mc_chem_df, nomc_chem_df, nomc_chem_genotox_df, edlist1_chem_df, edlist2_chem_df


def get_jump_profiles(jumpcp_df, metadata_df, inchikeys, feature_cols):
    """Median profile per (compound, imaging source) for the given INCHIKEYs."""
    jump_ids = metadata_df[metadata_df["INCHIKEY"].isin(inchikeys)]["JCP2022"].tolist()
    subset = jumpcp_df[jumpcp_df["Metadata_JCP2022"].isin(jump_ids)][
        ["Metadata_JCP2022", "Metadata_Source"] + feature_cols
    ]
    return subset.groupby(["Metadata_JCP2022", "Metadata_Source"], observed=True).median()


def main():
    mc_chem_df, nomc_chem_df, nomc_chem_genotox_df, edlist1_chem_df, edlist2_chem_df = load_chemical_data()
    metadata_df = pd.read_csv(METADATA_PATH, low_memory=False)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Combined MC + Non-MC chemical reference list (1,105 compounds: 276 MC,
    # 829 Non-MC - Methods 2.1), used by the ClusterMap and GBA notebooks.
    chemlist_df = pd.concat([mc_chem_df, nomc_chem_df])
    chemlist_df.to_csv(DATA_DIR / "chemlist_mc_nomc.csv", index=False)
    print(f"Chemical list:     {chemlist_df.shape}")

    jumpcp_df = pd.read_parquet(PARQUET_PATH)
    feature_cols = [c for c in jumpcp_df.columns if c.startswith("harmony_")]

    dmso_df = get_jump_profiles(jumpcp_df, metadata_df, [DMSO_INCHIKEY], feature_cols)
    dmso_df.to_csv(OUT_DIR / "dmso_jump_profiles.csv")
    print(f"DMSO:              {dmso_df.shape}")

    mc_df = get_jump_profiles(jumpcp_df, metadata_df, mc_chem_df["INCHIKEY"], feature_cols)
    mc_df.to_csv(OUT_DIR / "mc_jump_profiles.csv")
    print(f"MC:                {mc_df.shape}")

    nomc_df = get_jump_profiles(jumpcp_df, metadata_df, nomc_chem_df["INCHIKEY"], feature_cols)
    nomc_df.to_csv(OUT_DIR / "nomc_jump_profiles.csv")
    print(f"NOMC:              {nomc_df.shape}")

    nomc_genotox_df = get_jump_profiles(jumpcp_df, metadata_df, nomc_chem_genotox_df["INCHIKEY"], feature_cols)
    nomc_genotox_df.to_csv(OUT_DIR / "nomc_genotox_jump_profiles.csv")
    print(f"NOMC non-genotox:  {nomc_genotox_df.shape}")

    edlist1_df = get_jump_profiles(jumpcp_df, metadata_df, edlist1_chem_df["InChIKey"], feature_cols)
    edlist1_df.to_csv(OUT_DIR / "edlist1_jump_profiles.csv")
    print(f"EDList I:          {edlist1_df.shape}")

    edlist2_df = get_jump_profiles(jumpcp_df, metadata_df, edlist2_chem_df["InChIKey"], feature_cols)
    edlist2_df.to_csv(OUT_DIR / "edlist2_jump_profiles.csv")
    print(f"EDList II:         {edlist2_df.shape}")


if __name__ == "__main__":
    main()
