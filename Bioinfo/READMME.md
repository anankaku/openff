## Can OpenFF accurately capture the conformational preferences of peptoids, especially their backbone torsions and cis/trans amide behavior, compared with QM reference data?
- the project aims to understand if openff can descibe peptoid which is

### create openff env
- jupyter notebook: `1_read_sdf.py`
```bash
conda create -n openff_env python=3.10 -y
conda activate openff_env

conda install -c conda-forge openff-toolkit rdkit openmm -y
pip install nglview matplotlib
```
### read sdf and find the most conformation with phi and psi

```bash
import numpy as np
import matplotlib.pyplot as plt
from openff.toolkit import Molecule

def dihedral_angle(p0, p1, p2, p3):
    """Return dihedral angle in degrees."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so projection is stable
    b1 /= np.linalg.norm(b1)

    # vectors perpendicular to b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    angle = np.degrees(np.arctan2(y, x))
    return angle

mol = Molecule.from_file("s06.sdf")
mol.generate_conformers(n_conformers=200)

# phi atom index
dihedral_atoms = [1, 3, 17, 18]

angles = []

for conf in mol.conformers:
    coords = conf.m_as("angstrom")   # shape: (n_atoms, 3)
    i, j, k, l = dihedral_atoms
    angle = dihedral_angle(coords[i], coords[j], coords[k], coords[l])
    angles.append(angle)

print(angles)

plt.hist(angles, bins=18)
plt.xlabel("Dihedral Angle (degrees)")
plt.ylabel("Count")
plt.title("Dihedral Distribution")
plt.show()
```
### count the most conformation
```bash
import numpy as np

angles_array = np.array(angles)

counts, bins = np.histogram(angles_array, bins=18)

# find 
top2 = np.argsort(counts)[-2:]

peak_angles = [(bins[i] + bins[i+1]) / 2 for i in top2]

print("Top dihedral angles:")
for p in peak_angles:
    print(f"≈ {p:.2f}°")
```
### use openbabel create input file for Gaussian.
- install openbabel
```bash
conda install openbabel
```
- create input file
```bash
obabel my_peptoid.sdf -O my_peptoid.com --gen3d
```
### QM calculation

- edit input file for gaussian
- check multiplicity before submitting the job on que
```bash
python gaussian.py
```
- submit the input file
```bash
sbatch run_single.sh /file/name.com
```
### scan dihedral with Gaussian

- assign dihedral
- example:
```bash
%nprocshared=8
%mem=16GB
# opt=modredundant b3lyp/6-31g(d,p)

QM dihedral opt

0 1
 C                  0.70644900   -2.74280700   -0.00609300
 C                  1.81039800   -1.76208300    0.35536100
 O                  2.89271900   -2.17414000    0.77667500
......
......
......
 H                  6.99802400    0.92863900   -0.71248500

D 2 4 18 19 S 36 10.000000 # scan dihedral 36 steps with each 10 degree
```
- submit multiple job
```bash
bash run_all.sh
```