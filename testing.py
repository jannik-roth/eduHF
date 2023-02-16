import geometry
import BasisFunction
import numpy as np
import scf

at = geometry.Atom('H', np.array([1.0, 2.0, 3.0]))
at2 = geometry.Atom('H', np.array([1.0, 3.0, 4.0]))
mol = geometry.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.5), ('2s', 0.74)]}
basisset = BasisFunction.BasisSet(base)
basis = BasisFunction.Basis(mol, basisset, 1)

SCF = scf.SCF(mol, basis)
S, T, V, ERIs = SCF.prepare_integrals()
print(S)
#print(T)
#print(V)
#print(ERIs)
