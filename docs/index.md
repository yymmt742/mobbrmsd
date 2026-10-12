# mobbrmsd

**Molecular-Oriented RMSD with Branch-and-Bound**

[![CI](https://github.com/yymmt742/mobbrmsd/actions/workflows/ci.yml/badge.svg)](https://github.com/yymmt742/mobbrmsd/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-yymmt742.github.io%2Fmobbrmsd-blue)](https://yymmt742.github.io/mobbrmsd/)

`mobbrmsd` is an implementation of molecular-oriented RMSD (moRMSD)
using a branch-and-bound algorithm for efficient and exact
superposition of molecular assemblies.

It is designed for molecular systems containing multiple molecules,
including molecular liquids, solvated systems, and molecular assemblies.

## Features

- Exact molecular-oriented RMSD — permutation-invariant molecule comparison
- Branch-and-bound optimization — no brute-force enumeration
- Multi-species systems — each species can carry intramolecular symmetry
- Python API and command-line interface
- Fortran core with OpenMP / OpenBLAS support

## Getting started

- [Installation](installation.md)
- [Quick Start](quickstart.md)
- [Usage](usage/index.md)
- [API Reference](api/index.md)

## Algorithm

See [Theory](theory/index.md) for the mathematical definition
and algorithmic details.

## Reference

Yamamoto, Y.
"Algorithm for Efficient Superposition and Clustering of Molecular Assemblies Using the Branch-and-Bound Method"
*J. Chem. Inf. Model.* **2025**.
<https://doi.org/10.1021/acs.jcim.4c02217>
