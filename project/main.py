import pandas as pd
from pathlib import Path
import numpy as np

# =========================================================
# 📌 PATH SETUP
# Purpose: Ensure correct access to the data directory
# regardless of where the script is executed
# =========================================================

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent

input_file = BASE_DIR / "data" / "pubchem_brd4_bioactivity_protein.csv"
output_file = BASE_DIR / "data" / "brd4_full_dataset.csv"

# =========================================================
# 📌 CACHE CHECK
# If processed data already exists, skip reprocessing
# and load the existing dataset to save time
# =========================================================

if output_file.exists():
    print(f"--- {output_file.name} already exists. Loading cached data ---")
    final_df = pd.read_csv(output_file)

else:
    print(f"--- Processing raw dataset: {input_file.name} ---")

    # =========================================================
    # 📌 LOAD RAW DATA
    # PubChem bioactivity dataset
    # =========================================================

    df = pd.read_csv(
        input_file,
        sep=",",
        engine="python",
        on_bad_lines="skip",
        encoding="utf-8"
    )

    # =========================================================
    # 📌 COLUMN SELECTION
    # Only keep relevant bioactivity fields
    # =========================================================

    selected_columns = [
        "Activity",
        "Activity_Type",
        "Activity_Value",
        "Compound_CID"
    ]

    df_filtered = df[selected_columns].copy()

    # =========================================================
    # 📌 FILTER VALID ASSAY TYPES
    # Keep only standardized activity measurements
    # =========================================================

    df_filtered = df_filtered[
        df_filtered["Activity_Type"].isin(["IC50", "Ki", "Kd", "Unspecified"])
    ]

    # Remove missing activity values
    df_filtered = df_filtered[df_filtered["Activity_Value"].notna()]

    # =========================================================
    # 📌 UNIT ASSUMPTION
    # All values are assumed to be in µM
    # =========================================================

    df_filtered["uM"] = df_filtered["Activity_Value"]

    # =========================================================
    # 📌 ACTIVITY NORMALIZATION (p-value)
    # Log-scale transformation for model stability
    # =========================================================

    df_filtered["p_value"] = 6 - np.log10(df_filtered["Activity_Value"])

    # =========================================================
    # 📌 BINARY LABELING
    # Threshold: 1 µM
    # Active: ≤ 1 µM
    # Inactive: > 1 µM
    # =========================================================

    df_filtered["Label"] = (df_filtered["uM"] <= 1.0).astype(int)

    # =========================================================
    # 📌 UNSPECIFIED HANDLING
    # Assign labels based on threshold
    # =========================================================

    mask = df_filtered["Activity"] == "Unspecified"

    df_filtered.loc[mask, "Final_Activity"] = np.where(
        df_filtered.loc[mask, "uM"] <= 1.0,
        "Active",
        "Inactive"
    )

    # =========================================================
    # 📌 CLEANUP + DEDUPLICATION
    # Keep most reliable measurement per compound
    # =========================================================

    final_df = (
        df_filtered
        .sort_values(by="p_value", ascending=False)
        .drop_duplicates(subset=["Compound_CID"])
    )

    # =========================================================
    # 📌 MODEL DATA PREPARATION (PLACEHOLDER)
    # SELFIES will be generated in next pipeline step
    # =========================================================

    full_df = final_df.copy()
    full_df["SELFIES"] = None  # to be generated from SMILES

    model_df = full_df[["SELFIES", "Label"]].copy()

    model_df.to_csv(BASE_DIR / "data" / "brd4_model_dataset.csv", index=False)

    # Save full processed dataset
    output_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_file, index=False)

    print("--- Processing completed successfully ---")

print("Done.")