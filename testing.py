import geometry
import BasisFunction
import numpy as np
import scf

at = geometry.Atom('H', np.array([0.0, 0.0, 0.0]))
at2 = geometry.Atom('H', np.array([1.4, 0.0, 0.0]))
mol = geometry.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.24), ('2s', 1.80)]}
basisset = BasisFunction.BasisSet(base)
basis = BasisFunction.Basis(mol, basisset, 3)

SCF = scf.SCF(mol, basis)
SCF.prepare_integrals()
SCF.run_scf()

