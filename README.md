# Electron Density Decomposition onto Atom-Centered Basis

[![DOI](https://img.shields.io/badge/DOI-10.1039%2FC9SC02696G-blue)](https://doi.org/10.1039/C9SC02696G)

This code supports the paper
> A. Fabrizio, A. Grisafi, B. Meyer, M. Ceriotti, and C. Corminboeuf,<br>
> “Electron density learning of non-covalent systems”<br>
> [Chem. Sci. **10**, 9492 (2019)](https://doi.org/10.1039/C9SC02696G)

It is written to compute the density matrix of a given molecule withihn KS-DFT
and project its electron density onto a pre-selected Atom-Centered Basis.

## Requirements
* `python >= 3.6`
* `numpy >= 1.16`
* [`pyscf >= 1.6`](https://github.com/pyscf/pyscf)

## Usage

### 1. Decompose the electron density into atomic contributions
```
python deco.py --mol molecule --basis=basis --auxbasis=auxbasis --func=func
```
Computes the density matrix of a molecule at KS-DFT level 
and writes the decomposition coefficients (coeff.npy),
the density matrix (dm.npy) and the Overlap/Coulomb integrals (S.npy, J.npy).

#### Command-line arguments
* molecule: `.xyz` file with molecular geometry
* basis: AO basis
* auxbasis: density-fitting basis.
* func: KS-DFT functional to perform the computation.

#### Examples
```
python deco.py --mol ./Example/water.xyz --basis='cc-pvqz' --auxbasis='cc-pvqz-jkfit' --func='pbe'
python deco.py --mol ./Example/water.xyz --basis='6-31g' --auxbasis='cc-pvdz-jkfit' --func='pbe0'

```

