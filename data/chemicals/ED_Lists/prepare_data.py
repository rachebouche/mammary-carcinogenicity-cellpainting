#!/usr/bin/env python3
"""
Prepare ED list data:
  1. Concatenate the three Excel list files.
  2. Add a "List category" column: the lowest-numbered list a compound belongs
     to across all three files (List I > List II > List III in priority).
  3. Query PubChem (via pubchempy) for IsomericSMILES and InChIKey for each
     compound, trying CAS number first, then each ";" separated name as fallback.
  4. Save the enriched table to a CSV file.

Usage (dataproc conda env):
    conda run -n dataproc python prepare_data.py
"""

import time
import pandas as pd
import pubchempy as pcp

# ---------------------------------------------------------------------------
# 1. Load and concatenate the three source files
# ---------------------------------------------------------------------------

FILES = {
    "List I":   "ED_list1.xlsx",
    "List II":  "ED_list2.xlsx",
    "List III": "ED_list3.xlsx",
}

frames = []
for list_name, path in FILES.items():
    df = pd.read_excel(path)
    df["Source list"] = list_name
    # Normalise the "appears on lists" column name across files
    for col in df.columns:
        if "appears on lists" in col.lower():
            df = df.rename(columns={col: "Appears on lists"})
    # Normalise the year column name
    for col in df.columns:
        if col.lower() in ("year", "status year"):
            df = df.rename(columns={col: "Year"})
    frames.append(df)

combined = pd.concat(frames, ignore_index=True)

print(f"Combined table: {combined.shape[0]} rows, {combined.shape[1]} columns")
print("Columns:", combined.columns.tolist())

# ---------------------------------------------------------------------------
# 2. Add "List category": lowest-numbered list each compound belongs to
# ---------------------------------------------------------------------------

LIST_ORDER = {"List I": 1, "List II": 2, "List III": 3}
LIST_LABEL = {1: "List I", 2: "List II", 3: "List III"}

# Compound identity key: CAS number when available, otherwise name
def compound_key(row):
    cas = str(row["CAS no."]).strip()
    if cas and cas.lower() not in ("nan", ""):
        return cas
    return str(row["Name and abbreviation"]).strip()

combined["_compound_key"] = combined.apply(compound_key, axis=1)

# For each compound key, find the minimum list rank across all its rows
min_list_rank = (
    combined["Source list"]
    .map(LIST_ORDER)
    .groupby(combined["_compound_key"])
    .transform("min")
)
combined["List category"] = min_list_rank.map(LIST_LABEL)
combined.drop(columns=["_compound_key"], inplace=True)

print("\nList category distribution:")
print(combined["List category"].value_counts().sort_index())

# ---------------------------------------------------------------------------
# 3. Query PubChem for each unique compound
# ---------------------------------------------------------------------------

def query_pubchem(cas: str, name: str) -> dict:
    """Return dict with IsomericSMILES and InChIKey, or empty strings."""
    result = {"IsomericSMILES": "", "InChIKey": ""}

    # Helper to extract the two properties from a list of Compound objects
    def extract(compounds):
        if compounds:
            c = compounds[0]
            result["IsomericSMILES"] = c.isomeric_smiles or ""
            result["InChIKey"] = c.inchikey or ""
            return True
        return False

    # Try CAS number first (most unambiguous)
    cas_clean = str(cas).strip() if pd.notna(cas) else ""
    if cas_clean and cas_clean.lower() not in ("nan", ""):
        try:
            compounds = pcp.get_compounds(cas_clean, namespace="name")
            if extract(compounds):
                return result
        except Exception:
            pass
        time.sleep(0.2)  # be polite to the API

    # Fallback: try each name listed (separated by ";") until one matches
    name_clean = str(name).strip() if pd.notna(name) else ""
    if name_clean and name_clean.lower() not in ("nan", ""):
        for candidate in name_clean.split(";"):
            candidate = candidate.strip()
            if not candidate:
                continue
            try:
                compounds = pcp.get_compounds(candidate, namespace="name")
                if extract(compounds):
                    return result
            except Exception:
                pass
            time.sleep(0.2)

    return result


# Build a lookup table keyed on (CAS, Name) to avoid duplicate API calls
unique_ids = combined[["CAS no.", "Name and abbreviation"]].drop_duplicates()
total = len(unique_ids)
print(f"\nQuerying PubChem for {total} unique compound entries …")

lookup: dict[tuple, dict] = {}
for idx, (_, row) in enumerate(unique_ids.iterrows(), start=1):
    key = (row["CAS no."], row["Name and abbreviation"])
    lookup[key] = query_pubchem(row["CAS no."], row["Name and abbreviation"])
    if idx % 20 == 0 or idx == total:
        found = sum(1 for v in lookup.values() if v["IsomericSMILES"])
        print(f"  {idx}/{total} processed — {found} matched so far")

# ---------------------------------------------------------------------------
# 4. Attach the new columns and save
# ---------------------------------------------------------------------------

combined["IsomericSMILES"] = combined.apply(
    lambda r: lookup.get((r["CAS no."], r["Name and abbreviation"]), {}).get("IsomericSMILES", ""),
    axis=1,
)
combined["InChIKey"] = combined.apply(
    lambda r: lookup.get((r["CAS no."], r["Name and abbreviation"]), {}).get("InChIKey", ""),
    axis=1,
)

output_path = "ED_lists_combined.csv"
combined.to_csv(output_path, index=False)

matched = (combined["IsomericSMILES"] != "").sum()
print(f"\nDone. {matched}/{len(combined)} rows have a SMILES.")
print(f"Output saved to: {output_path}")
