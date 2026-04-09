from pathlib import Path
import pandas as pd
from openff.toolkit import Molecule
import os

INPUT_CSV = "molecules_clean.csv"
OUTDIR = Path("step2_output")
OUTDIR.mkdir(exist_ok=True)

print("Current working directory:", os.getcwd())
print("Reading:", Path(INPUT_CSV).resolve())

df = pd.read_csv(INPUT_CSV)

mapped_rows = []
atom_rows = []

for _, row in df.iterrows():
    mol_id = str(row["id"])
    smiles = str(row["input_smiles"]).strip()

    mol = Molecule.from_smiles(smiles, name=mol_id)

    mapped_smiles = mol.to_smiles(
        isomeric=True,
        explicit_hydrogens=True,
        mapped=True
    )

    mapped_mol = Molecule.from_mapped_smiles(mapped_smiles)
    mapped_mol.name = mol_id

    mapped_rows.append({
        "id": mol_id,
        "mapped_smiles": mapped_smiles
    })

    for atom in mapped_mol.atoms:
        neighbors = ",".join(str(n.molecule_atom_index) for n in atom.bonded_atoms)

        atom_rows.append({
            "id": mol_id,
            "atom_index": atom.molecule_atom_index,
            "element": atom.symbol,
            "formal_charge": atom.formal_charge.m,
            "is_aromatic": atom.is_aromatic,
            "neighbors": neighbors
        })

mapped_path = OUTDIR / "mapped_smiles.csv"
atom_path = OUTDIR / "atom_table.csv"

pd.DataFrame(mapped_rows).to_csv(mapped_path, index=False)
pd.DataFrame(atom_rows).to_csv(atom_path, index=False)

print("Wrote:", mapped_path.resolve())
print("Wrote:", atom_path.resolve())