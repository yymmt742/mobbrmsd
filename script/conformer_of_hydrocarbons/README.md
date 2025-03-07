# Conformer of hydrocarbons

### Overview

Here we perform conformer structure fitting using saturated hydrocarbons.
Equivalent atoms of saturated hydrocarbons often occur due to C-C bond rotations and molecular symmetry.
For example, The hydrogen atoms of a methyl group can be mapped to three by C-C bonding.

In order to find the optimal superposition between conformers,
optimal atomic mappings between equivalent atoms must be found,
a problem that often leads to combinatorial explosion.
The number of atomic mappings increases as an exponential power with a base of 3 relative to the number of methyl groups.

### Prerequisites

The environment python should have the following dependent libraries available:

* RDKit
* mdtraj
* numpy
* pandas

The following command line tools should be available from the bash environment:

* openbabel (obabel)

### Test molecules

* ethane (C2H6)
* propane (C3H8)
* isobutane (C4H10)
* neopentane (C5H12)
* tetramethylbutane (C6H14)
* tetramethylhexane (C8H16)
* tetramethylpentane (C10H22)
* hexamethylpentane  (C11H24)
* octamethylhexane   (C14H30)

You can run the calculations via:

```bash
source runtest.sh
```
