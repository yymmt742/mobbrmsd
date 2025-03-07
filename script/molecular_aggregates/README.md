# Molecular aggregates

### Overview

Here we perform superposing of molecular aggregates.
The model molecules are water (H2O), ammonia (NH3), and methane (CH4),
and the system contains M homologous molecules.
Molecular aggregation systems are constructed by applying rigid body transformations to molecular templates containing hydrogen atoms.




The procedure is as follows:

- A single molecular structure is loaded with a template.
- The rigid-body transformation is performed for each molecule
  using random rotations from SO(3),
  and a translation vector based on the normal distribution N(0, M 1/3).
- The wall-clock times taken by mobbRMSD to find the best atomic mapping and the RMSD values are reported.
- For comparison, the following two previous methods are also run simultaneously:

* GetBestRMS (implementation in [RDKit](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html#rdkit.Chem.rdMolAlign.GetBestRMS))
* LP-RMSD (implementation in [MDtraj](https://mdtraj.org/1.9.4/api/generated/mdtraj.lprmsd.html))

### Prerequisites

The environment python should have the following dependent libraries available:

* RDKit
* mdtraj
* numpy
* scipy
* pandas
* matplotlib

### Test molecules

* water (H2O)
* ammonia (NH3)
* methane (CH4)

You can run the calculations via:

```bash
source runtest.sh
```
