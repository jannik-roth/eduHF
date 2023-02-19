import numpy as np
import eduHF

at = eduHF.Atom('H', np.array([0.0, 0.0, 0.0]))
at2 = eduHF.Atom('F', np.array([1.733, 0.0, 0.0]))
mol = eduHF.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.24)],
        'Li': [('1s', 2.69), ('2s', 0.75)],
        'Be': [('1s', 3.68), ('2s', 1.10)],
        'B' : [('1s', 4.68), ('2s', 1.45), ('2p', 1.45)],
        'C' : [('1s', 5.67), ('2s', 1.72), ('2p', 1.72)],
        'N' : [('1s', 6.67), ('2s', 1.95), ('2p', 1.95)],
        'O' : [('1s', 7.66), ('2s', 2.25), ('2p', 2.25)],
        'F' : [('1s', 8.65), ('2s', 2.55), ('2p', 2.55)],
        }
basisset = eduHF.BasisSet(base)
basis = eduHF.Basis(mol, basisset, 3)

SCF = eduHF.SCF(mol, basis, max_iter=20)
print(basis.nbf)
SCF.prepare_integrals()
SCF.run_scf()

