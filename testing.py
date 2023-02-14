import geometry
import BasisFunction
import slater_expansion
import numpy as np
import mcmurchiedavidson

def kinetic_item(bf1 : BasisFunction.ContractedGaussianFunction,
                 bf2 : BasisFunction.ContractedGaussianFunction):
    return mcmurchiedavidson.kinetic(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                     bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)

def overlap_item(bf1 : BasisFunction.ContractedGaussianFunction,
                 bf2 : BasisFunction.ContractedGaussianFunction):
    return mcmurchiedavidson.overlap(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                     bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)

def potential_1e_item(bf1 : BasisFunction.ContractedGaussianFunction,
                      bf2 : BasisFunction.ContractedGaussianFunction,
                      at : geometry.Atom):
    return mcmurchiedavidson.potential_1e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, 
                                          bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                          at.xyz)

def potential_2e_item(bf1 : BasisFunction.ContractedGaussianFunction,
                      bf2 : BasisFunction.ContractedGaussianFunction,
                      bf3 : BasisFunction.ContractedGaussianFunction,
                      bf4 : BasisFunction.ContractedGaussianFunction):
    return mcmurchiedavidson.potential_2e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                          bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                          bf3.coeffs, bf3.alphas, bf3.l_vec, bf3.xyz,
                                          bf4.coeffs, bf4.alphas, bf4.l_vec, bf4.xyz,)

def overlap_matrix(basis : BasisFunction.Basis):
    nbf = basis.nof_gauss_funcs()
    S = np.eye(nbf)
    for mu in range(nbf):
        for nu in range(mu+1, nbf):
            S[mu, nu] = overlap_item(basis.basis_gauss[mu], basis.basis_gauss[nu])
            S[nu, mu] = S[mu, nu]
    return S

def kinetic_matrix(basis : BasisFunction.Basis):
    nbf = basis.nof_gauss_funcs()
    T = np.zeros((nbf,nbf))
    for mu in range(nbf):
        for nu in range(mu, nbf):
            T[mu, nu] = kinetic_item(basis.basis_gauss[mu], basis.basis_gauss[nu])
            T[nu, mu] = T[mu, nu]
    return T

def potential_1e_matrix(basis : BasisFunction.Basis, mol : geometry.Molecule):
    nbf = basis.nof_gauss_funcs()
    V = np.zeros((nbf, nbf))
    for mu in range(nbf):
        for nu in range(mu, nbf):
            for at in mol.geometry:
                V[mu, nu] += - at.charge * potential_1e_item(basis.basis_gauss[mu], basis.basis_gauss[nu], at)
                V[nu, mu] = V[mu, nu]
    return V

def potential_2e_tensor(basis : BasisFunction.Basis):
    nbf = basis.nof_gauss_funcs()
    eris = np.zeros((nbf, nbf, nbf, nbf))
    for mu in range(nbf):
        for nu in range(nbf):
            for lamda in range(nbf):
                for sigma in range(nbf):
                    eris[mu, nu, lamda, sigma] = potential_2e_item(basis.basis_gauss[mu],
                                                                   basis.basis_gauss[nu],
                                                                   basis.basis_gauss[lamda],
                                                                   basis.basis_gauss[sigma],)
    return eris
            



at = geometry.Atom('H', np.array([1.0, 2.0, 3.0]))
at2 = geometry.Atom('H', np.array([1.0, 3.0, 4.0]))
mol = geometry.Molecule([at, at2], 0)

base = {'H' : [('1s', 1.5), ('2s', 0.74)]}
basisset = BasisFunction.BasisSet(base)
basis = BasisFunction.Basis(mol, basisset, 1)

print(overlap_matrix(basis))
print(kinetic_matrix(basis))
print(potential_1e_matrix(basis, mol))
print(potential_2e_tensor(basis))
