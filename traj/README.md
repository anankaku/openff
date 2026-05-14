# Peptoid Topology Generation, System Assembly, and Folding@home Preparation

This repository documents a workflow for preparing peptoid molecular dynamics (MD) systems starting from SMILES strings and ending with simulation-ready files that can be tested locally in OpenMM or GROMACS before being prepared for Folding@home production MD.

The overall goal is not to calculate a topology on Folding@home. Instead, the topology, coordinates, force-field parameters, and simulation settings are prepared locally first. Folding@home then runs production MD trajectories using validated input files.

---

## 1. Project Goal

The purpose of this workflow is to prepare peptoid systems for force-field benchmarking and large-scale MD sampling.

For each peptoid structure, the pipeline should produce:

```text
SMILES / molecular structure
→ 3D conformer
→ force-field parameters
→ topology and coordinates
→ assembled MD system
→ local minimization / short MD validation
→ Folding@home-ready simulation package
→ trajectory analysis
```

The expected analysis may include:

* Peptoid conformational sampling
* Backbone torsional distributions
* cis/trans population analysis using the omega dihedral
* Comparison between force fields
* Comparison against DFT torsional profiles or literature MD results

---

## 2. Recommended Repository Structure

```text
peptoid-fah-pipeline/
├── README.md
├── environment.yml
├── molecules.csv
│
├── 00_input/
│   ├── peptoid_001.smi
│   ├── peptoid_001.sdf
│   └── peptoid_001.pdb
│
├── 01_build_structure/
│   └── smiles_to_sdf.py
│
├── 02_openff_parameterization/
│   ├── openff_vacuum_minimize.py
│   ├── openff_short_md.py
│   └── outputs/
│
├── 03_amber_gromacs_parameterization/
│   ├── run_antechamber.sh
│   ├── tleap.in
│   ├── convert_amber_to_gromacs.py
│   └── outputs/
│
├── 04_system_assembly/
│   ├── assemble_openmm_system.py
│   ├── minimization.py
│   └── short_md.py
│
├── 05_validation/
│   ├── md.log
│   ├── minimized.pdb
│   └── traj.dcd
│
├── 06_fah_package/
│   ├── system.xml
│   ├── integrator.xml
│   ├── state.xml
│   ├── system.top
│   ├── system.gro
│   ├── mdp/
│   └── README_fah_package.md
│
└── 07_analysis/
    ├── analyze_omega.py
    ├── analyze_dihedrals.py
    └── plot_cis_trans.py
```

---

## 3. Software Requirements

This workflow may use either the OpenFF/OpenMM route or the AmberTools/GROMACS route.

### Core tools

* Python
* RDKit
* OpenFF Toolkit
* OpenMM
* AmberTools
* ParmEd
* GROMACS

### Optional tools

* MDTraj
* MDAnalysis
* NumPy
* Pandas
* Matplotlib

---

## 4. Conda Environment

Create a conda environment for OpenFF and OpenMM:

```bash
conda create -n openff -c conda-forge \
    python=3.11 \
    openff-toolkit \
    openff-interchange \
    openmm \
    rdkit \
    parmed \
    mdtraj \
    mdanalysis \
    numpy pandas matplotlib -y

conda activate openff
```

For AmberTools:

```bash
conda create -n ambertools -c conda-forge ambertools parmed -y
conda activate ambertools
```

Check installation:

```bash
python -c "import openmm; print(openmm.__version__)"
python -c "from openff.toolkit import Molecule; print('OpenFF OK')"
antechamber -h
parmchk2 -h
```

---

## 5. Input Format

Prepare a CSV file containing peptoid IDs and capped SMILES strings.

Example: `molecules.csv`

```csv
id,Capped_SMILES
1,CN(C)C(=O)CN(C)Cc1ccccc1F
```

Each molecule should have a unique ID so that output files can be named consistently.

---

## 6. Step 1 — Generate 3D Structure from SMILES

This step converts a SMILES string into a 3D structure file such as SDF and PDB.

Create `01_build_structure/smiles_to_sdf.py`:

```python
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

INPUT_CSV = "../molecules.csv"
OUTPUT_DIR = Path("../00_input")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

molecules = pd.read_csv(INPUT_CSV)

for _, row in molecules.iterrows():
    mol_id = str(row["id"])
    smiles = row["Capped_SMILES"]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES for molecule {mol_id}: {smiles}")

    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol, randomSeed=2026)
    if status != 0:
        raise RuntimeError(f"3D embedding failed for molecule {mol_id}")

    AllChem.MMFFOptimizeMolecule(mol)

    sdf_path = OUTPUT_DIR / f"peptoid_{mol_id}.sdf"
    pdb_path = OUTPUT_DIR / f"peptoid_{mol_id}.pdb"
    smi_path = OUTPUT_DIR / f"peptoid_{mol_id}.smi"

    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    Chem.MolToPDBFile(mol, str(pdb_path))

    with open(smi_path, "w") as f:
        f.write(smiles + "\n")

    print(f"Generated files for peptoid {mol_id}")
    print(f"  SDF: {sdf_path}")
    print(f"  PDB: {pdb_path}")
```

Run:

```bash
cd 01_build_structure
python smiles_to_sdf.py
```

Expected output:

```text
00_input/peptoid_1.sdf
00_input/peptoid_1.pdb
00_input/peptoid_1.smi
```

---

## 7. Step 2A — OpenFF Parameterization for OpenMM

This route uses OpenFF to assign force-field parameters and create an OpenMM `System`.

Create `02_openff_parameterization/openff_vacuum_minimize.py`:

```python
from pathlib import Path
from openff.toolkit import Molecule, ForceField
from openmm import unit
from openmm import LangevinIntegrator
from openmm.app import Simulation, PDBReporter, StateDataReporter

INPUT_SDF = Path("../00_input/peptoid_1.sdf")
OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

mol = Molecule.from_file(str(INPUT_SDF))

forcefield = ForceField("openff-2.2.0.offxml")
topology = mol.to_topology()
system = forcefield.create_openmm_system(topology)

integrator = LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.001 * unit.picoseconds,
)

simulation = Simulation(
    topology.to_openmm(),
    system,
    integrator,
)

simulation.context.setPositions(mol.conformers[0])

state = simulation.context.getState(getEnergy=True)
print("Initial potential energy:", state.getPotentialEnergy())

simulation.minimizeEnergy()

state = simulation.context.getState(getEnergy=True, getPositions=True)
print("Minimized potential energy:", state.getPotentialEnergy())

with open(OUTPUT_DIR / "minimized.pdb", "w") as f:
    from openmm.app import PDBFile
    PDBFile.writeFile(
        topology.to_openmm(),
        state.getPositions(),
        f,
    )
```

Run:

```bash
cd 02_openff_parameterization
python openff_vacuum_minimize.py
```

Expected success criteria:

* OpenFF can read the molecule.
* OpenMM can create a system.
* Initial potential energy can be evaluated.
* Energy minimization finishes without crashing.

This is the first important validation step.

---

## 8. Step 2B — AmberTools / GAFF2 Parameterization

This route generates AMBER topology and coordinate files, which can later be converted to GROMACS.

Create `03_amber_gromacs_parameterization/run_antechamber.sh`:

```bash
#!/bin/bash
set -e

INPUT_SDF="../00_input/peptoid_1.sdf"
OUTPUT_DIR="outputs"
mkdir -p ${OUTPUT_DIR}

antechamber \
    -i ${INPUT_SDF} \
    -fi sdf \
    -o ${OUTPUT_DIR}/peptoid_1.mol2 \
    -fo mol2 \
    -c bcc \
    -s 2

parmchk2 \
    -i ${OUTPUT_DIR}/peptoid_1.mol2 \
    -f mol2 \
    -o ${OUTPUT_DIR}/peptoid_1.frcmod
```

Run:

```bash
cd 03_amber_gromacs_parameterization
bash run_antechamber.sh
```

Expected output:

```text
outputs/peptoid_1.mol2
outputs/peptoid_1.frcmod
```

---

## 9. Step 3 — Generate AMBER Topology with tleap

Create `03_amber_gromacs_parameterization/tleap.in`:

```text
source leaprc.gaff2

MOL = loadmol2 outputs/peptoid_1.mol2
loadamberparams outputs/peptoid_1.frcmod

saveamberparm MOL outputs/peptoid_1.prmtop outputs/peptoid_1.inpcrd
savepdb MOL outputs/peptoid_1_amber.pdb

quit
```

Run:

```bash
tleap -f tleap.in
```

Expected output:

```text
outputs/peptoid_1.prmtop
outputs/peptoid_1.inpcrd
outputs/peptoid_1_amber.pdb
```

These files define the molecule using AMBER-style topology and coordinates.

---

## 10. Step 4 — Convert AMBER Files to GROMACS Files

Create `03_amber_gromacs_parameterization/convert_amber_to_gromacs.py`:

```python
import parmed as pmd
from pathlib import Path

OUTPUT_DIR = Path("outputs")

amber = pmd.load_file(
    str(OUTPUT_DIR / "peptoid_1.prmtop"),
    str(OUTPUT_DIR / "peptoid_1.inpcrd"),
)

amber.save(str(OUTPUT_DIR / "peptoid_1.top"), overwrite=True)
amber.save(str(OUTPUT_DIR / "peptoid_1.gro"), overwrite=True)

print("Generated GROMACS files:")
print(OUTPUT_DIR / "peptoid_1.top")
print(OUTPUT_DIR / "peptoid_1.gro")
```

Run:

```bash
python convert_amber_to_gromacs.py
```

Expected output:

```text
outputs/peptoid_1.top
outputs/peptoid_1.gro
```

---

## 11. Step 5 — Assemble the Simulation System

System assembly means building a complete MD system, not just a single molecule.

A complete system may include:

```text
peptoid
+ water box
+ ions
+ topology
+ force-field parameters
+ simulation settings
```

At this stage, the goal is to create files that can be minimized and simulated locally before being prepared for Folding@home.

---

## 12. Step 5A — OpenMM System Assembly

For OpenMM, the final target files may include:

```text
system.xml
integrator.xml
state.xml
```

A simple vacuum system can be serialized first. Solvated systems should be tested after the vacuum system works.

Create `04_system_assembly/serialize_openmm_system.py`:

```python
from pathlib import Path
from openff.toolkit import Molecule, ForceField
from openmm import unit, LangevinIntegrator, XmlSerializer
from openmm.app import Simulation

INPUT_SDF = Path("../00_input/peptoid_1.sdf")
OUTPUT_DIR = Path("../06_fah_package/openmm")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

mol = Molecule.from_file(str(INPUT_SDF))
forcefield = ForceField("openff-2.2.0.offxml")

topology = mol.to_topology()
system = forcefield.create_openmm_system(topology)

integrator = LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

simulation = Simulation(topology.to_openmm(), system, integrator)
simulation.context.setPositions(mol.conformers[0])
simulation.minimizeEnergy()

state = simulation.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getParameters=True,
)

with open(OUTPUT_DIR / "system.xml", "w") as f:
    f.write(XmlSerializer.serialize(system))

with open(OUTPUT_DIR / "integrator.xml", "w") as f:
    f.write(XmlSerializer.serialize(integrator))

with open(OUTPUT_DIR / "state.xml", "w") as f:
    f.write(XmlSerializer.serialize(state))

print("Serialized OpenMM files written to:", OUTPUT_DIR)
```

Run:

```bash
cd 04_system_assembly
python serialize_openmm_system.py
```

Expected output:

```text
06_fah_package/openmm/system.xml
06_fah_package/openmm/integrator.xml
06_fah_package/openmm/state.xml
```

---

## 13. Step 5B — GROMACS System Assembly

For GROMACS, the final target files include:

```text
system.top
system.gro
minim.mdp
nvt.mdp
prod.mdp
```

Example minimization file: `06_fah_package/gromacs/mdp/minim.mdp`

```text
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000

nstlist         = 10
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
```

Test minimization locally:

```bash
gmx grompp \
    -f mdp/minim.mdp \
    -c peptoid_1.gro \
    -p peptoid_1.top \
    -o em.tpr \
    -maxwarn 2

gmx mdrun -deffnm em
```

Expected output:

```text
em.gro
em.edr
em.log
em.trr
```

If this step fails, fix the topology, coordinates, or force-field assignment before moving forward.

---

## 14. Step 6 — Run a Short Local MD Test

Before preparing Folding@home files, the system must pass a local short MD test.

### OpenMM short MD test

Create `04_system_assembly/short_md.py`:

```python
from pathlib import Path
from openff.toolkit import Molecule, ForceField
from openmm import unit, LangevinIntegrator
from openmm.app import Simulation, DCDReporter, StateDataReporter

INPUT_SDF = Path("../00_input/peptoid_1.sdf")
OUTPUT_DIR = Path("../05_validation")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

mol = Molecule.from_file(str(INPUT_SDF))
forcefield = ForceField("openff-2.2.0.offxml")

topology = mol.to_topology()
system = forcefield.create_openmm_system(topology)

integrator = LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

simulation = Simulation(topology.to_openmm(), system, integrator)
simulation.context.setPositions(mol.conformers[0])

simulation.minimizeEnergy()

simulation.reporters.append(
    StateDataReporter(
        str(OUTPUT_DIR / "md.log"),
        100,
        step=True,
        potentialEnergy=True,
        temperature=True,
        speed=True,
    )
)

simulation.reporters.append(
    DCDReporter(str(OUTPUT_DIR / "traj.dcd"), 100)
)

simulation.step(5000)  # 10 ps with a 2 fs timestep

print("Short MD completed successfully.")
```

Run:

```bash
python short_md.py
```

Expected output:

```text
05_validation/md.log
05_validation/traj.dcd
```

Success criteria:

* No NaN energies
* No exploding coordinates
* Temperature remains physically reasonable
* Energy does not diverge immediately
* Trajectory can be opened and visualized

---

## 15. Step 7 — Prepare Folding@home-Ready Files

Folding@home should be treated as the production MD platform. The topology and assembled system are inputs.

### What should be prepared locally

For OpenMM-style Folding@home projects:

```text
system.xml
integrator.xml
state.xml
```

For GROMACS-style Folding@home projects:

```text
system.top
system.gro
production.mdp
frame0.tpr
```

### Important note

A Folding@home project normally requires project-side setup, such as:

* Project number
* Work server configuration
* Project XML file
* Work unit generation
* Internal testing
* Beta testing
* Public release

This is usually handled by the PI, project manager, or Folding@home project-side infrastructure. The goal of this repository is to prepare and validate the molecular system before it is packaged into work units.

---

## 16. Example Folding@home Package Layout

```text
06_fah_package/
├── README_fah_package.md
│
├── openmm/
│   ├── system.xml
│   ├── integrator.xml
│   └── state.xml
│
└── gromacs/
    ├── peptoid_1.top
    ├── peptoid_1.gro
    ├── mdp/
    │   ├── minim.mdp
    │   ├── nvt.mdp
    │   └── prod.mdp
    └── frame0.tpr
```

---

## 17. Step 8 — Analyze Folding@home or Local MD Trajectories

After trajectories are generated, analyze backbone torsions and cis/trans populations.

For peptoids, the omega dihedral is especially important:

```text
omega near 0°      → cis-like
omega near ±180°   → trans-like
```

Create `07_analysis/analyze_omega.py`:

```python
import mdtraj as md
import numpy as np
import pandas as pd
from pathlib import Path

TRAJ = "../05_validation/traj.dcd"
TOP = "../00_input/peptoid_1.pdb"
OUTPUT = Path("omega_analysis.csv")

traj = md.load(TRAJ, top=TOP)

# Replace these atom indices with the actual omega dihedral atoms.
# MDTraj uses zero-based atom indexing.
omega_atoms = [0, 1, 2, 3]

angles_rad = md.compute_dihedrals(traj, [omega_atoms])[:, 0]
angles_deg = np.degrees(angles_rad)

# Normalize to [-180, 180]
angles_deg = ((angles_deg + 180) % 360) - 180

cis_mask = np.abs(angles_deg) < 45
trans_mask = np.abs(np.abs(angles_deg) - 180) < 45

n_cis = np.sum(cis_mask)
n_trans = np.sum(trans_mask)

kct = n_cis / n_trans if n_trans > 0 else np.nan

summary = pd.DataFrame({
    "n_frames": [len(angles_deg)],
    "n_cis": [n_cis],
    "n_trans": [n_trans],
    "Kct_cis_over_trans": [kct],
})

summary.to_csv(OUTPUT, index=False)
print(summary)
```

The atom indices must be updated for each molecule.

---

## 18. Validation Checklist

Before sending files for Folding@home preparation, check the following:

### Structure generation

* [ ] SMILES can be parsed.
* [ ] 3D conformer is generated.
* [ ] Hydrogens are present.
* [ ] Structure is chemically reasonable.

### Parameterization

* [ ] OpenFF or GAFF2 parameters are assigned.
* [ ] No missing atom types.
* [ ] No missing bonds, angles, or dihedrals.
* [ ] Charges are assigned.

### Topology generation

* [ ] OpenMM system can be created.
* [ ] AMBER `prmtop/inpcrd` files are generated if using AmberTools.
* [ ] GROMACS `top/gro` files are generated if using GROMACS.

### System assembly

* [ ] System can be minimized.
* [ ] Short MD runs without crashing.
* [ ] Energies are finite.
* [ ] Temperature is stable.
* [ ] Trajectory can be analyzed.

### Folding@home preparation

* [ ] OpenMM XML files or GROMACS TPR files are available.
* [ ] Initial state is minimized or equilibrated.
* [ ] Simulation protocol is documented.
* [ ] Expected returned trajectory files are specified.
* [ ] Project-side Folding@home setup is coordinated with the PI or project manager.

---

## 19. Summary of the Workflow

```text
1. Prepare SMILES input
2. Generate 3D structure using RDKit
3. Parameterize the peptoid using OpenFF or AmberTools/GAFF2
4. Generate topology and coordinate files
5. Assemble a complete MD system
6. Test minimization locally
7. Run short local MD
8. Prepare OpenMM or GROMACS files for Folding@home
9. Run production MD on Folding@home
10. Analyze trajectories for torsions, cis/trans populations, and conformational sampling
```

---

## 20. Current Minimum Deliverable

The first deliverable should be one successfully validated peptoid system.

Minimum files:

```text
peptoid_1.sdf
peptoid_1.pdb
minimized.pdb
md.log
traj.dcd
system.xml
```

or, for GROMACS:

```text
peptoid_1.mol2
peptoid_1.frcmod
peptoid_1.prmtop
peptoid_1.inpcrd
peptoid_1.top
peptoid_1.gro
em.gro
em.log
```

Once this works for one peptoid, the same workflow can be applied to the full set of peptoids for force-field benchmarking.

---

## 21. Notes

* Topology generation is a preparation step, not the Folding@home calculation itself.
* Folding@home runs production MD trajectories.
* The topology, coordinates, and force-field parameters are inputs to production MD.
* Local validation is required before large-scale production runs.
* For benchmarking, the same molecule should be tested across different force fields under consistent simulation conditions.

---

## 22. References

* Harris, B. S.; Bejagam, K. K.; Baer, M. D. *Development of a Systematic and Extensible Force Field for Peptoids (STEPs).* J. Phys. Chem. B 2023.
* Muli, C. S.; Xie, D.; Post, C. B.; Trader, D. J. *NMR-Guided Studies to Establish the Binding Interaction between a Peptoid and Protein.* J. Am. Chem. Soc. 2025.
* OpenFF Toolkit Documentation.
* OpenMM Documentation.
* AmberTools Documentation.
* GROMACS Documentation.
* Folding@home Work Server Documentation.
