import numpy as np
import eduHF

at = eduHF.Atom('H', np.array([0.0, 0.0, 0.0]))
at2 = eduHF.Atom('H', np.array([1.4, 0.0, 0.0]))
mol = eduHF.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.24), ('2s', 1.80)]}
basisset = eduHF.BasisSet(base)
basis = eduHF.Basis(mol, basisset, 3)

SCF = eduHF.SCF(mol, basis)
SCF.prepare_integrals()
SCF.run_scf()

