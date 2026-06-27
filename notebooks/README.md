# OSET — Oblate Spheroid Excitation Theory

Symbolic derivations and numerical verification for:

> **Oblate Spheroid Excitation Theory: A Unified, Lattice-Free Foundation for
> Plastic Deformation from Which Dislocations Emerge as Collective Excitations**
> Albert Linda and K. A. Padmanabhan.

This repository accompanies the manuscript's *Data Availability* statement. It
contains a single Jupyter notebook that reproduces, symbolically and
numerically, every closed-form result reported in the paper.

## Contents

- [OSET_derivations.ipynb](OSET_derivations.ipynb) — Full symbolic (SymPy) +
  numerical (NumPy/SciPy) verification of the OSET derivations. Sections 1–17
  check the final boxed formulas (closed forms, limiting cases, numeric tables)
  against the manuscript; Sections 18–35 re-derive each formula step by step
  from its stated starting equations.

## Requirements

- Python 3.12 (tested; 3.10+ should work)
- Packages listed in [requirements.txt](requirements.txt): NumPy, SciPy, SymPy,
  and Jupyter.

Install them with:

```bash
pip install -r requirements.txt
```

## Running

Open the notebook interactively:

```bash
jupyter notebook OSET_derivations.ipynb
```

or execute it headlessly:

```bash
jupyter nbconvert --to notebook --execute OSET_derivations.ipynb
```