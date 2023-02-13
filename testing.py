import geometry
import BasisFunction
import slater_expansion
import numpy as np
import mcmurchiedavidson

#h1 = [['H', [1, 2, 3]],
#    ]

#print(BasisFunction._Lvecs_from_nl("4d"))

#testing = [BasisFunction.SlaterFunction(2.0,"2p",0), BasisFunction.SlaterFunction(1.0,"1s",0), BasisFunction.SlaterFunction(3.0,"3s",0)]
#print(testing)

#for Slater in testing:
#    print(BasisFunction.Slater2CGaussians(Slater, 3))

#print(BasisFunction.Slaters2CGaussians(testing, 1))

def kinetic_item(bf1 : BasisFunction.ContractedGaussianFunction, bf2 : BasisFunction.ContractedGaussianFunction):
    return mcmurchiedavidson.kinetic(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)

def overlap_item(bf1 : BasisFunction.ContractedGaussianFunction, bf2 : BasisFunction.ContractedGaussianFunction):
    return mcmurchiedavidson.overlap(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)

def overlap_matrix(basis : BasisFunction.Basis):
    nbf = basis.nof_gauss_funcs()
    S = np.eye(nbf)
    for i in range(nbf):
        for j in range(i+1, nbf):
            S[i, j] = overlap_item(basis.basis_gauss[i], basis.basis_gauss[j])
            S[j, i] = S[i, j]
    return S

def kinetic_matrix(basis : BasisFunction.Basis):
    nbf = basis.nof_gauss_funcs()
    T = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(i, nbf):
            T[i, j] = kinetic_item(basis.basis_gauss[i], basis.basis_gauss[j])
            T[j, i] = T[i, j]
    return T



at = geometry.Atom('H', np.array([1.0, 2.0, 3.0]))
at2 = geometry.Atom('H', np.array([1.0, 3.0, 4.0]))
mol = geometry.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.5), ('2s', 0.74)]}
basisset = BasisFunction.BasisSet(base)
basis = BasisFunction.Basis(mol, basisset, 1)

print(overlap_matrix(basis))
print(kinetic_matrix(basis))
