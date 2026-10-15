# Fortran Core API

The Fortran core is documented with [FORD](https://forddocs.readthedocs.io/).
It covers the internal modules that implement the branch-and-bound search,
the RMSD kernel, and the Python–Fortran interface layer.

[Open Fortran API →](https://yymmt742.github.io/mobbrmsd/fortran-api/){ .md-button .md-button--primary }

## Source coverage

| Directory | Description |
|-----------|-------------|
| `src/mobbrmsd/` | Branch-and-bound solver, RMSD kernels |
| `src/python_interface/` | Fortran ↔ Python wrapper (`iso_c_binding`) |

## Dependency graphs

FORD generates call-graph and module-dependency diagrams (via Graphviz) for
every public procedure.  They are embedded in each entity's documentation page.

!!! note "Versioning"
    The Fortran API reference is rebuilt and deployed together with the Python
    docs on every push to `main` and on every release tag.
