import geometry
import BasisFunction
import slater_expansion
import numpy as np

#h1 = [['H', [1, 2, 3]],
#    ]

#print(BasisFunction._Lvecs_from_nl("4d"))

#testing = [BasisFunction.SlaterFunction(2.0,"2p",0), BasisFunction.SlaterFunction(1.0,"1s",0), BasisFunction.SlaterFunction(3.0,"3s",0)]
#print(testing)

#for Slater in testing:
#    print(BasisFunction.Slater2CGaussians(Slater, 3))

#print(BasisFunction.Slaters2CGaussians(testing, 1))

at = geometry.Atom('H', np.array([1.0, 2.0, 3.0]))
at2 = geometry.Atom('H', np.array([1.0, 3.0, 4.0]))
mol = geometry.Molecule([at, at2], -1)

base = {'H' : [('1s', 1.5), ('2p', 3.0)]}
basisset = BasisFunction.BasisSet(base)

basis = BasisFunction.Basis(mol, basisset, 3)
print(basis.nof_slater_funcs())
print(basis.nof_gauss_funcs())