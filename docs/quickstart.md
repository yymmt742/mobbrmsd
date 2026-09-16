# Quick Start

## Run the demonstration

After installing `mobbrmsd`, run:

```bash
python -m mobbrmsd demo
```

## Run an RMSD calculation

The command-line interface accepts a JSON input file:

```bash
python -m mobbrmsd run -i input.json
```

For details of the input format, see
[Input format](usage/input.md).

## Python interface

The Python API provides data classes for describing molecular systems
and functions for calculating molecular-oriented RMSD.

See the [API Reference](api/index.md).

