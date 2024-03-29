# eduHF - educational Hartree-Fock program

This repository is an educational implementation of a Hartree-Fock (HF) program which can be used for simple quantum chemical calculations. The program is mainly written in Python an uses `numba` to compile the integral routines just-in-time (`jit`). At the moment, only a restricted HF calculation is possible.

*Note: The code is **not** optimized! It serves merely as an easy to follow implementation*

## 0. Theoretical Background

An explanation on the Hartree-Fock method, geometry optimization and DIIS is given in `docs/theory.pdf`. Literature recommendations are given there too.

The program uses the McMurchie-Davidson Scheme for the evaluation of the necessary integrals and their derivatives with respect to atomic positions. A short intrduction is given in `docs/MMD.pdf`.

## 1. Installation
Download the code to your local machine. To use the package, you might need to add the location of the package to your `PYTHONPATH` by adding
```
export PYTHONPATH=$PYTHONPATH:/path/to/package/eduHF
```
(and changing the path accordingly) to your `.bashrc`.

## 2. Usage

For basic usage, please see the `example_usage.ipynb`. Else, you can always have a look into the code, it is not too complex.