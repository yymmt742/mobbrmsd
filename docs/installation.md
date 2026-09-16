# Installation

## Requirements

For the Python interface:

- Python >= 3.8
- pip

The native implementation additionally requires a Fortran compiler.

## Install from GitHub

```bash
pip install git+https://github.com/yymmt742/mobbrmsd.git
```

## Install from source

Clone the repository:

```bash
git clone https://github.com/yymmt742/mobbrmsd.git
cd mobbrmsd
```

Then install:

```bash
pip install .
```

## Development installation

For development, install the package in editable mode:

```bash
pip install -e .
```

## Documentation development

Install the documentation dependencies:

```bash
python -m pip install mkdocs mkdocs-material "mkdocstrings[python]"
```

Then start the local documentation server:

```bash
mkdocs serve
```

The documentation will be available at:

```text
http://127.0.0.1:8000/
```
