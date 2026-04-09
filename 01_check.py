import pandas as pd
from openff.toolkit import Molecule

INPUT_CSV = "s01.csv"

df = pd.read_csv(INPUT_CSV)

required_cols = {"id", "Capped_SMILES"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

ok_rows = []
bad_rows = []

for _, row in df.iterrows():
    mol_id = str(row["id"])
    smiles = str(row["Capped_SMILES"]).strip()

    try:
        mol = Molecule.from_smiles(smiles, name=mol_id)

        ok_rows.append(
            {
                "id": mol_id,
                "input_smiles": smiles,
                "canonical_smiles": mol.to_smiles(),
                "n_atoms": mol.n_atoms,
                "name": mol.name,
            }
        )

    except Exception as e:
        bad_rows.append(
            {
                "id": mol_id,
                "input_smiles": smiles,
                "error": str(e),
            }
        )

ok_df = pd.DataFrame(ok_rows)
bad_df = pd.DataFrame(bad_rows)

ok_df.to_csv("molecules_clean.csv", index=False)
bad_df.to_csv("molecules_failed.csv", index=False)

print(f"Total rows: {len(df)}")
print(f"Passed: {len(ok_df)}")
print(f"Failed: {len(bad_df)}")

if len(ok_df) > 0:
    print("\nFirst successful entry:")
    print(ok_df.iloc[0].to_dict())

if len(bad_df) > 0:
    print("\nFirst failed entry:")
    print(bad_df.iloc[0].to_dict())