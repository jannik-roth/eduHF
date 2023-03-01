# eduHF - educational Hartree-Fock program

This repository is an educational implementation of a Hartree-Fock (HF) program which can be used for simple quantum chemical calculations.

*Note: The code is **not** optimized! It serves merely as an easy to follow implementation*

The program uses the McMurchie-Davidson Scheme for the evaluation of the necessary integrals. A short intrduction is given in `docs/MMD.pdf`. At the moment, only a restricted HF calculation is possible.

## 1. Installation
Download the code to your local machine. To use the package, you might need to add the location of the package to your `PYTHONPATH` by adding
```
export PYTHONPATH=$PYTHONPATH:/path/to/package/eduHF
```
(and changing of course the path accordingly) to your `.bashrc`.

The package uses `f2py`. It is therefore necessary to first compile the parts written in Fortran. This can be achieved by executing the following commands.
```
cd path/to/package/eduHF/eduHF
python3 -m numpy.f2py -c slater_expansion.f90 -m slater_expansion
python3 -m numpy.f2py -c mcmurchie_davidson.f90 -m mcmurchie_davidson
```
Alternatively, you can execute the `build_fortran.sh` script which executes the commands for you.

## 2. Usage

For basic usage, please see the `example_usage.ipynb`. Else, you can always have a look into the code, it is not too complex and should be easily readible.