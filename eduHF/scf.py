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

    def geom_opt(self,
                 max_iter_geom : int = 50,
                 opt_method : str = 'gd',
                 convergence_crit_geom : float = 1e-5,
                 convergence_type_geom : str = 'rms',
                 info : int = 10,
                 scf_params : dict = {'info' : 1},
                 opt_method_params : dict = {}):
        
        self.max_iter_geom = max_iter_geom
        self.opt_method = opt_method
        self.convergence_crit_geom = convergence_crit_geom
        self.convergence_type_geom = convergence_type_geom
        self.req_iterations_geom = None
        self.converged_geom = False

        if info > 5:
            self.print_geom_opt_setup()

        for iter in range(self.max_iter_geom):
            if info > 7:
                self.print_geom_opt_step(iter)
            self.prepare_integrals()
            self.run_scf(**scf_params)
            self.grad = self.get_gradient()

            self.converged_geom, error = self.check_convergence_geom()

            if info > 7:
                self.print_grad_error(error)

            if self.converged_geom:
                if info > 1:
                    self.req_iterations_geom = iter
                    self.print_conv_geom_info()
                break
            
            step = self.build_step(iter, **opt_method_params)
            print(step)
            self.make_step(step)

    def print_conv_geom_info(self):
        print("#"*self.print_width)
        print(f"GEOM OPT CONVERGED".center(self.print_width))
        print(f'req. {self.req_iterations_geom} iterations'.center(self.print_width))
        print("#"*self.print_width)

    def print_grad_error(self, error):
        print(f"Current Gradient".center(self.print_width))
        print("{:<10} {:<10} {:<10} {:<10}".format('symbol','x','y', 'z'))
        for idx, at in enumerate(self.mol.geometry):
            print("{:<10} {:<10} {:<10} {:<10}".format(at.symbol, f"{self.grad[idx*3]:.7f}", f"{self.grad[idx*3+1]:.7f}", f"{self.grad[idx*3+2]:.7f}"))
        print(f"Error in geom opt: {error:e}")

    def print_geom_opt_step(self, iter):
        print("#"*self.print_width)
        print(f"GEOM OPT ITER {iter}".center(self.print_width))
        print(f"Current Geometry".center(self.print_width))
        self.print_mol(self.mol)

    @staticmethod
    def print_mol(mol):
        print("{:<10} {:<10} {:<10} {:<10}".format('symbol','x','y', 'z'))
        for at in mol.geometry:
            print("{:<10} {:<10} {:<10} {:<10}".format(at.symbol, f"{at.xyz[0]:.7f}", f"{at.xyz[1]:.7f}", f"{at.xyz[2]:.7f}"))


    def print_geom_opt_setup(self):
        print("#"*self.print_width)
        print("Geom opt setup".center(self.print_width))
        print("#"*self.print_width)
        print(f"  max_iter_geom          -> {self.max_iter_geom}")
        print(f"  opt_method             -> {self.opt_method}")
        print(f"  convergence_crit_geom  -> {self.convergence_crit_geom:e}")
        print(f"  convergence_type_geom  -> {self.convergence_type_geom}")
        print("Starting geom opt now".center(self.print_width))
        print("#"*self.print_width)
            
    def make_step(self, step):
        for at in range(self.mol.nofatoms):
            self.mol.geometry[at].xyz += step[at*3:at*3+3]

    def build_step(self, iter, **kwargs):
        if self.opt_method == 'gd':
            alpha = kwargs.get('alpha', 1.0)
            return - alpha * self.grad
        if self.opt_method == 'bfgs':
            if iter == 0:
                self.inv_hessian = np.eye(len(self.grad))
            else:
                diff_grad = self.grad - self.grad_old 
                tmp = np.eye((len(self.grad))) - np.outer(self.step_old, diff_grad) / np.dot(self.step_old, diff_grad)
                self.inv_hessian = tmp @ self.inv_hessian @ tmp.T  + np.outer(self.step_old, self.step_old) / np.dot(self.step_old, diff_grad)       

            alpha = kwargs.get('alpha', 1.0)
            self.step = - alpha * self.inv_hessian @ self.grad
            self.grad_old = self.grad
            self.step_old = self.step

            return self.step

    def check_convergence_geom(self):
        if self.convergence_type_geom == 'rms':
            rms = np.linalg.norm(self.grad)
            return rms < self.convergence_crit_geom, rms

    def run_scf(self, 
                max_iter : int = 50,
                convergence_crit : float = 1e-6,
                convergence_type : str = 'com',
                info : int = 10,
                diis : bool = True,
                diis_size : int = 6):
        
        self.max_iter = max_iter
        self.convergence_type = convergence_type
        self.convergence_crit = convergence_crit
        self.info = info
        self.mol.core_pot = self.mol.core_potential()
        self.diis = diis
        self.diis_size = diis_size

        
        if not self.setup_int:
            raise ValueError(f"Integrals are not set up! Please run {self.__class__.__name__}.setup_int() before running the SCF procedure")
        
        if info > 5:
            self.print_scf_setup()
        
        # get X = S^(-1/2)
        self.obtain_X()
        self.inital_Fock_guess(type='zero density')
        self.F_prime = self.X.T @ self.F @ self.X

        if info > 7:
            print("{:<10} {:<15} {:<15} {:<10}".format('iteration','energy','error', 'converged'))
        
        for i in range(self.max_iter):
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
            if info > 7:
                print("{:<10} {:<15} {:<15} {:<10}".format(i, f"{self.scf_energy:,.7e}", f"{self.error:,.7e}", str(self.converged)))

            if self.converged:
                self.req_iterations = i
                if info > 1:
                    self.print_conv_info()
                break

            # F' = X^T F X
            self.F_prime = self.X.T @ self.F @ self.X

            if self.diis == True:
                self.save_fock_prime(i)
                self.diis_build_Fock_prime()
                
        if self.converged == False:
            if info > 0:
                self.print_not_conv_info()

    def diis_build_Fock_prime(self):
        current_diis_size = self.diis_error_com.size
        if current_diis_size > 2:
            B = -1.0 * np.ones((current_diis_size + 1, current_diis_size + 1))
            B[-1, -1] = 0.0
            for i in range(current_diis_size):
                for j in range(i, current_diis_size):
                    B[i, j] = np.trace(self.diis_error_mat_save[i,:,:].T @ self.diis_error_mat_save[j,:,:])
                    B[j, i] = B[i, j]
                
            rhs = np.zeros(current_diis_size + 1)
            rhs[-1] = -1.0
            sol = np.linalg.solve(B, rhs)

            self.F_prime = np.einsum('i,ijk', sol[:-1], self.diis_fock_save)


    def save_fock_prime(self, i):
        if i == 0:
            self.diis_fock_save = np.zeros((0, self.basis.nbf, self.basis.nbf))
            self.diis_error_mat_save = np.zeros((0, self.basis.nbf, self.basis.nbf))
            self.diis_error_com = np.zeros(0)
        else:
            err_mat = self.X.T @ (self.F @ self.P @ self.S - self.S @ self.P @ self.F) @ self.X
            err_com = np.linalg.norm(err_mat)

            if self.diis_error_com.shape[0] < self.diis_size:
                self.diis_fock_save = np.vstack((self.diis_fock_save, self.F_prime[None]))
                self.diis_error_mat_save = np.vstack((self.diis_error_mat_save, err_mat[None]))
                self.diis_error_com = np.append(self.diis_error_com, err_com)
            else:
                idx_to_replace = np.argmax(self.diis_error_com)
                self.diis_fock_save[idx_to_replace] = self.F_prime
                self.diis_error_mat_save[idx_to_replace] = err_mat
                self.diis_error_com[idx_to_replace] = err_com            

    def get_gradient(self):
        # check converged
        Q = self.build_Q()
        grad = np.zeros((self.mol.nofatoms * 3))

        for at in range(self.mol.nofatoms):
            for dim in range(3):
                grad[at * 3 + dim] = self.get_derivative(at, dim, Q)
        return grad

    def get_derivative(self, atom, dim, Q):
        H_core_d = self.kinetic_der_matrix(atom, dim) + self.potential_1e_der_matrix(atom, dim)
        ERIS_d = self.potential_2e_der_tensor(atom, dim)
        S_d = self.overlap_der_matrix(atom, dim)
        core_pot_d = self.mol.core_potential_der(atom, dim)

        der = (np.einsum('nm,mn', self.P, H_core_d) 
               + 0.5 * (np.einsum('nm,ls,mnsl', self.P, self.P, ERIS_d) 
                        - 0.5 * np.einsum('nm,ls,mlsn', self.P, self.P, ERIS_d)) 
               - np.einsum('nm,mn', Q, S_d)
               + core_pot_d)
        return der

    def build_Q(self):
        nbf = self.basis.nbf
        Q = np.zeros((nbf, nbf))
        for mu in range(nbf):
            for nu in range(nbf):
                for i in range(int(self.mol.noe / 2)):
                    Q[mu, nu] += self.epsilons[i] * self.C[mu, i] * self.C[nu, i]
        Q *= 2
        return Q

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
                        V_der[mu, nu] -= -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 0, dim)
                        V_der[mu, nu] -= -at.charge * SCF.potential_1e_der_item(self.basis(mu), self.basis(nu), at, 1, dim)
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
    
    def potential_2e_der_tensor(self, center, dim):
        nbf = self.basis.nbf
        eris_der = np.zeros((nbf, nbf, nbf, nbf))
        for mu in range(nbf):
            for nu in range(nbf):
                for lamda in range(nbf):
                    for sigma in range(nbf):
                        if ((center == self.basis(mu).center) + (center == self.basis(nu).center)
                             + (center == self.basis(lamda).center) + (center == self.basis(sigma).center) == 0):
                            eris_der[mu, nu, lamda, sigma] = 0.0
                        elif ((center == self.basis(mu).center) + (center == self.basis(nu).center)
                             + (center == self.basis(lamda).center) + (center == self.basis(sigma).center) == 4):
                            eris_der[mu, nu, lamda, sigma] = 0.0
                        elif ((center == self.basis(mu).center) + (center == self.basis(nu).center)
                             + (center == self.basis(lamda).center) + (center == self.basis(sigma).center) == 3):
                            if not (center == self.basis(mu).center):
                                eris_der[mu, nu, lamda, sigma] -= SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            0, dim)
                            elif not (center == self.basis(nu).center):
                                eris_der[mu, nu, lamda, sigma] -= SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            1, dim)
                            elif not (center == self.basis(lamda).center):
                                eris_der[mu, nu, lamda, sigma] -= SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            2, dim)
                            elif not (center == self.basis(sigma).center):
                                eris_der[mu, nu, lamda, sigma] -= SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            3, dim)
                        elif ((center == self.basis(mu).center) + (center == self.basis(nu).center)
                             + (center == self.basis(lamda).center) + (center == self.basis(sigma).center) <= 2):
                            if (center == self.basis(mu).center):
                                eris_der[mu, nu, lamda, sigma] += SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            0, dim)
                            if (center == self.basis(nu).center):
                                eris_der[mu, nu, lamda, sigma] += SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            1, dim)
                            if (center == self.basis(lamda).center):
                                eris_der[mu, nu, lamda, sigma] += SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            2, dim)
                            if (center == self.basis(sigma).center):
                                eris_der[mu, nu, lamda, sigma] += SCF.potential_2e_der_item(self.basis(mu),
                                                                                            self.basis(nu),
                                                                                            self.basis(lamda),
                                                                                            self.basis(sigma),
                                                                                            3, dim)
        return eris_der
    
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
    
    @staticmethod
    def potential_2e_der_item(bf1 : ContractedGaussianFunction,
                              bf2 : ContractedGaussianFunction,
                              bf3 : ContractedGaussianFunction,
                              bf4 : ContractedGaussianFunction,
                              center : int,
                              dim : int):
        return potential_2e_der(bf1.coeffs, bf1.alphas, bf1.l_vec, bf1.xyz,
                                bf2.coeffs, bf2.alphas, bf2.l_vec, bf2.xyz,
                                bf3.coeffs, bf3.alphas, bf3.l_vec, bf3.xyz,
                                bf4.coeffs, bf4.alphas, bf4.l_vec, bf4.xyz,
                                center, dim)