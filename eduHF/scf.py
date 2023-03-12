import numpy as np
from .mcmurchie_davidson import *
from .geometry import *
from .basisfunction import *


class SCF:
    def __init__(self,
                 mol : Molecule,
                 basis : Basis):
        self.mol = mol
        self.basis = basis

        self.S, self.T, self.V = None, None, None
        self.X, self.H, self.F_prime = None, None, None
        self.C, self.P = None, None
        self.ERIs = None
        self.epsilons = None
        self.scf_energy = None
        self.req_iterations = None
        self.error = None

        self.setup_int = False
        self.converged = False

        self.print_width = 60

    def run_scf(self, 
                max_iter : int = 50,
                convergence_crit : float = 1e-6,
                convergence_type : str = 'com'):
        
        self.max_iter = max_iter
        self.convergence_type = convergence_type
        self.convergence_crit = convergence_crit
        
        if not self.setup_int:
            raise ValueError(f"Integrals are not set up! Please run {self.__class__.__name__}.setup_int() before running the SCF procedure")
            return
        
        self.print_scf_setup()
        
        # get X = S^(-1/2)
        self.obtain_X()
        self.inital_Fock_guess(type='zero density')

        print("{:<10} {:<15} {:<15} {:<10}".format('iteration','energy','error', 'converged'))
        
        for i in range(self.max_iter):
            # F' = X^T F X
            self.F_prime = self.X.T @ self.F @ self.X
            # C' F' = C' epsilon
            self.epsilons, C_prime = np.linalg.eigh(self.F_prime)
            # C = X C'
            self.C = self.X @ C_prime
            # P_mn = \sum_i C_mi C_ni
            self.build_P()
            # Fock build
            self.build_Fock()
            # calc energy
            self.scf_energy = self.calc_HF_energy()

            self.converged, self.error = self.check_convergence()
            print("{:<10} {:<15} {:<15} {:<10}".format(i, f"{self.scf_energy:,.7e}", f"{self.error:,.7e}", str(self.converged)))

            if self.converged:
                self.req_iterations = i
                self.print_conv_info()
                break
        
        if self.converged == False:
            self.print_not_conv_info()

    def print_scf_setup(self):
        print("#"*self.print_width)
        print("SCF setup".center(self.print_width))
        print("#"*self.print_width)
        print(f"  max_iter          -> {self.max_iter}")
        print(f"  convergence_crit  -> {self.convergence_crit:e}")
        print(f"  convergence_type  -> {self.convergence_type}")
        print("Starting SCF now".center(self.print_width))
        print("#"*self.print_width)

    def print_not_conv_info(self):
        print("#"*self.print_width)
        print(f'    SCF did NOT converge after {self.max_iter} iterations'.center(self.print_width))
        print("#"*self.print_width)

    def print_conv_info(self):
        print("#"*self.print_width)
        print("SCF converged".center(self.print_width))
        print(f'req. {self.req_iterations} iterations'.center(self.print_width))
        print("#"*self.print_width)
        print(f'SCF energy    = {self.scf_energy:,.10f}')
        print(f'Nuc Repulsion = {self.mol.core_pot:,.10f}')
        print(f'Final Energy  = {(self.scf_energy + self.mol.core_pot):,.10f}')

    def check_convergence(self):
        if self.convergence_type == 'com':
            val = np.linalg.norm(self.X.T @ (self.F @ self.P @ self.S - self.S @ self.P @ self.F) @ self.X)
            return val < self.convergence_crit, val
        else:
            raise ValueError(f"I do not know the convergence_type {self.convergence_type}")

    def calc_HF_energy(self):
        return 0.5 * np.sum(np.multiply(self.P, self.T + self.V + self.F))

    def build_Fock(self):
        J = np.einsum('ls,mnls->mn', self.P, self.ERIs, optimize=True)
        K = -0.5 * np.einsum('ls,mlsn->mn', self.P, self.ERIs, optimize=True)
        self.F = self.T + self.V + J + K

    def build_P(self):
        for mu in range(self.basis.nbf):
            for nu in range(self.basis.nbf):
                self.P[mu,nu] = 0
                for i in range(int(self.mol.noe / 2)):
                    self.P[mu, nu] += self.C[nu, i] * self.C[mu, i]
        self.P = 2*self.P

    def inital_Fock_guess(self, type='zero density'):
        if (type == 'zero density'):
            self.P = np.zeros_like(self.S)
            self.F = self.T + self.V
        else:
            raise ValueError(f"I do not know the initial Fock guess type: {type}")

    def obtain_X(self):
        eps, vecs = np.linalg.eigh(self.S)
        self.X = vecs @ np.diag(eps**(-0.5)) @ vecs.T

    def prepare_integrals(self):
        self.S = self.overlap_matrix()
        self.T = self.kinetic_matrix()
        self.V = self.potential_1e_matrix()
        self.ERIs = self.potential_2e_tensor()
        self.setup_int = True

    def overlap_matrix(self):
        nbf = self.basis.nbf
        S = np.eye(nbf)
        for mu in range(nbf):
            for nu in range(mu+1, nbf):
                S[mu, nu] = SCF.overlap_item(self.basis(mu), self.basis(nu))
                S[nu, mu] = S[mu, nu]
        return S
    
    def overlap_der_matrix(self, center, dim):
        nbf = self.basis.nbf
        S_der = np.zeros((nbf,nbf))
        for mu in range(nbf):
            for nu in range(mu+1, nbf):
                if ((self.basis(mu).center == center) and (self.basis(nu).center == center)):
                    S_der[mu, nu] = 0.0
                elif self.basis(mu).center == center:
                    S_der[mu, nu] = SCF.overlap_der_item(self.basis(mu), self.basis(nu), 0, dim)
                elif self.basis(nu).center == center:
                    S_der[mu, nu] = SCF.overlap_der_item(self.basis(mu), self.basis(nu), 1, dim)
                else:
                    S_der[mu, nu] = 0.0
                S_der[nu, mu] = S_der[mu, nu]
        return S_der

    def kinetic_matrix(self):
        nbf = self.basis.nbf
        T = np.zeros((nbf,nbf))
        for mu in range(nbf):
            for nu in range(mu, nbf):
                T[mu, nu] = SCF.kinetic_item(self.basis(mu), self.basis(nu))
                T[nu, mu] = T[mu, nu]
        return T
    
    def kinetic_der_matrix(self, center, dim):
        nbf = self.basis.nbf
        T_der = np.zeros((nbf,nbf))
        for mu in range(nbf):
            for nu in range(mu+1, nbf):
                if ((self.basis(mu).center == center) and (self.basis(nu).center == center)):
                    T_der[mu, nu] = 0.0
                elif self.basis(mu).center == center:
                    T_der[mu, nu] = SCF.kinetic_der_item(self.basis(mu), self.basis(nu), 0, dim)
                elif self.basis(nu).center == center:
                    T_der[mu, nu] = SCF.kinetic_der_item(self.basis(mu), self.basis(nu), 1, dim)
                else:
                    T_der[mu, nu] = 0.0
                T_der[nu, mu] = T_der[mu, nu]
        return T_der
    
    def potential_1e_matrix(self):
        nbf = self.basis.nbf
        V = np.zeros((nbf, nbf))
        for mu in range(nbf):
            for nu in range(mu, nbf):
                for at in self.mol.geometry:
                    V[mu, nu] += - at.charge * SCF.potential_1e_item(self.basis(mu), self.basis(nu), at)
                    V[nu, mu] = V[mu, nu]
        return V
    
    def potential_1e_der_matrix(self, center, dim):
        nbf = self.basis.nbf
        V_der = np.zeros((nbf, nbf))
        for mu in range(nbf):
            for nu in range(mu, nbf):
                for idx, at in enumerate(self.mol.geometry):
                    if ((self.basis(mu).center == center) and (self.basis(nu).center == center) and (idx == center)):
                        # all at the same center
                        V_der[mu, nu] += 0.0
                    elif ((self.basis(mu).center == center) and (idx == center)):
                        V_der[mu, nu] -= -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 1, dim)
                    elif ((self.basis(nu).center == center) and (idx == center)):
                        V_der[mu, nu] -= -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 0, dim)
                    elif ((self.basis(mu).center == center) and (self.basis(nu).center == center)):
                        V_der[mu, nu] += -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 0, dim)
                        V_der[mu, nu] += -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 1, dim)
                    elif (idx == center):
                        # use chain rule, factor 2.0 because we differentiate wrt to mu AND nu
                        V_der[mu, nu] -= -2.0 * at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 0, dim)
                    elif (self.basis(mu).center == center):
                        V_der[mu, nu] += -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 0, dim)
                    elif (self.basis(nu).center == center):
                        V_der[mu, nu] += -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 1, dim)
                    else:
                        V_der[mu, nu] += 0.0
                V_der[nu, mu] = V_der[mu, nu]
        return V_der

    def potential_2e_tensor(self):
        nbf = self.basis.nbf
        eris = np.zeros((nbf, nbf, nbf, nbf))
        for mu in range(nbf):
            for nu in range(nbf):
                for lamda in range(nbf):
                    for sigma in range(nbf):
                        eris[mu, nu, lamda, sigma] = SCF.potential_2e_item(self.basis(mu),
                                                                           self.basis(nu),
                                                                           self.basis(lamda),
                                                                           self.basis(sigma))
        return eris
    
    @staticmethod
    def kinetic_item(bf1 : ContractedGaussianFunction,
                     bf2 : ContractedGaussianFunction):
        return kinetic(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                       bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)
    
    @staticmethod
    def kinetic_der_item(bf1 : ContractedGaussianFunction,
                         bf2 : ContractedGaussianFunction,
                         center : int,
                         dim : int):
        return kinetic_der(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                           bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                           center, dim)
    @staticmethod
    def overlap_item(bf1 : ContractedGaussianFunction,
                     bf2 : ContractedGaussianFunction):
        return overlap(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                       bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz)
    @staticmethod
    def overlap_der_item(bf1 : ContractedGaussianFunction,
                         bf2 : ContractedGaussianFunction,
                         center : int,
                         dim : int):
        return overlap_der(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                           bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                           center, dim)
    @staticmethod
    def potential_1e_item(bf1 : ContractedGaussianFunction,
                          bf2 : ContractedGaussianFunction,
                          at : Atom):
        return potential_1e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, 
                            bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                            at.xyz)
    
    @staticmethod
    def potential_1e_der_item(bf1 : ContractedGaussianFunction,
                              bf2 : ContractedGaussianFunction,
                              at : Atom,
                              center : int,
                              dim : int):
        return potential_1e_der(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz, 
                                bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                at.xyz, center, dim)
    @staticmethod
    def potential_2e_item(bf1 : ContractedGaussianFunction,
                          bf2 : ContractedGaussianFunction,
                          bf3 : ContractedGaussianFunction,
                          bf4 : ContractedGaussianFunction):
        return potential_2e(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                            bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                            bf3.coeffs, bf3.alphas, bf3.l_vec, bf3.xyz,
                            bf4.coeffs, bf4.alphas, bf4.l_vec, bf4.xyz)