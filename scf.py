import numpy as np
import geometry
import BasisFunction
import mcmurchiedavidson

class SCF:
    def __init__(self,
                 mol : geometry.Molecule,
                 basis : BasisFunction.Basis):
        self.mol = mol
        self.basis = basis

    def prepare_integrals(self):
        S = self.overlap_matrix()
        T = self.kinetic_matrix()
        V = self.potential_1e_matrix()
        ERIs = self.potential_2e_tensor()
        return S, T, V, ERIs

    def overlap_matrix(self):
        nbf = self.basis.nof_gauss_funcs()
        S = np.eye(nbf)
        for mu in range(nbf):
            for nu in range(mu+1, nbf):
                S[mu, nu] = SCF.overlap_item(self.basis.basis_gauss[mu], self.basis.basis_gauss[nu])
                S[nu, mu] = S[mu, nu]
        return S

    def kinetic_matrix(self):
        nbf = self.basis.nof_gauss_funcs()
        T = np.zeros((nbf,nbf))
        for mu in range(nbf):
            for nu in range(mu, nbf):
                T[mu, nu] = SCF.kinetic_item(self.basis.basis_gauss[mu], self.basis.basis_gauss[nu])
                T[nu, mu] = T[mu, nu]
        return T
    
    def potential_1e_matrix(self):
        nbf = self.basis.nof_gauss_funcs()
        V = np.zeros((nbf, nbf))
        for mu in range(nbf):
            for nu in range(mu, nbf):
                for at in self.mol.geometry:
                    V[mu, nu] += - at.charge * SCF.potential_1e_item(self.basis.basis_gauss[mu], self.basis.basis_gauss[nu], at)
                    V[nu, mu] = V[mu, nu]
        return V

    def potential_2e_tensor(self):
        nbf = self.basis.nof_gauss_funcs()
        eris = np.zeros((nbf, nbf, nbf, nbf))
        for mu in range(nbf):
            for nu in range(nbf):
                for lamda in range(nbf):
                    for sigma in range(nbf):
                        eris[mu, nu, lamda, sigma] = SCF.potential_2e_item(self.basis.basis_gauss[mu],
                                                                           self.basis.basis_gauss[nu],
                                                                           self.basis.basis_gauss[lamda],
                                                                           self.basis.basis_gauss[sigma])
        return eris
    
    @staticmethod
    def kinetic_item(bf1 : BasisFunction.ContractedGaussianFunction,
                     bf2 : BasisFunction.ContractedGaussianFunction):
        return mcmurchiedavidson.kinetic(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                         bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)
    @staticmethod
    def overlap_item(bf1 : BasisFunction.ContractedGaussianFunction,
                     bf2 : BasisFunction.ContractedGaussianFunction):
        return mcmurchiedavidson.overlap(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                         bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)
    @staticmethod
    def potential_1e_item(bf1 : BasisFunction.ContractedGaussianFunction,
                          bf2 : BasisFunction.ContractedGaussianFunction,
                          at : geometry.Atom):
        return mcmurchiedavidson.potential_1e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, 
                                              bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                              at.xyz)
    @staticmethod
    def potential_2e_item(bf1 : BasisFunction.ContractedGaussianFunction,
                          bf2 : BasisFunction.ContractedGaussianFunction,
                          bf3 : BasisFunction.ContractedGaussianFunction,
                          bf4 : BasisFunction.ContractedGaussianFunction):
        return mcmurchiedavidson.potential_2e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                              bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                              bf3.coeffs, bf3.alphas, bf3.l_vec, bf3.xyz,
                                              bf4.coeffs, bf4.alphas, bf4.l_vec, bf4.xyz,)