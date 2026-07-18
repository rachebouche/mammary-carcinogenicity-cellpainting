#!/usr/bin/env python3
"""
Build the primary MC / Non-MC chemical reference lists.

Source: Kay & Rudel (2024) mammary carcinogenicity compendium
(KAY_RUDEL_2024_CombinedSuppTables_rev2.xlsx). Table S1 lists mammary
carcinogens (MC), Table S5 lists tested compounds with no evidence of
mammary tumour induction (Non-MC). Each compound is looked up in PubChem
(via PubChemPy, trying DTXSID, then CASRN, then preferred name) to attach
an InChIKey and isomeric SMILES.

Outputs:
  BreastCancerChemList_MC_v2.csv
  BreastCancerChemList_NOMC_v2.csv

Usage:
    python3 prepare_chemical_lists.py
"""
import pandas as pd
import pubchempy as pcp

SOURCE_XLSX = "KAY_RUDEL_2024_CombinedSuppTables_rev2.xlsx"
KEEP_COLUMNS = [
    "CASRN", "DTXSID", "preferred_name", "MammaryTumorEvidence",
    "MammaryTumorRefs", "Genotoxicity", "HormoneSummary",
]


def add_inchikey_and_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Attach INCHIKEY/SMILES from PubChem, trying DTXSID, CASRN, then name."""
    inchikeys, smiles = [], []
    for _, row in df.iterrows():
        match = None
        for identifier in (row["DTXSID"], row["CASRN"], row["preferred_name"]):
            try:
                match = pcp.get_compounds(identifier, "name")[0]
                break
            except Exception:
                continue
        inchikeys.append(match.inchikey if match else None)
        smiles.append(match.isomeric_smiles if match else None)
    df = df.copy()
    df["INCHIKEY"] = inchikeys
    df["SMILES"] = smiles
    return df


def main():
    mc_raw = pd.read_excel(SOURCE_XLSX, sheet_name="Excel Table S1", header=1)
    nomc_raw = pd.read_excel(SOURCE_XLSX, sheet_name="Excel Table S5", header=1)

    mc_df = mc_raw[mc_raw["MammaryTumorEvidence"] == "MC"][KEEP_COLUMNS]
    nomc_df = nomc_raw[nomc_raw["MammaryTumorEvidence"] == "Bioassay_noMC"][KEEP_COLUMNS]

    mc_df = add_inchikey_and_smiles(mc_df)
    nomc_df = add_inchikey_and_smiles(nomc_df)

    mc_df.to_csv("BreastCancerChemList_MC_v2.csv", index=False)
    nomc_df.to_csv("BreastCancerChemList_NOMC_v2.csv", index=False)
    print(f"MC:     {len(mc_df)} compounds -> BreastCancerChemList_MC_v2.csv")
    print(f"Non-MC: {len(nomc_df)} compounds -> BreastCancerChemList_NOMC_v2.csv")


if __name__ == "__main__":
    main()
