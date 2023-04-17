import numpy as np
from dataclasses import dataclass
from scipy.special import factorial2
from .geometry import Molecule
from . import slater_expansion

@staticmethod
def _Lvecs_from_nl(nl_quant: str):
    conv = {"s" : [np.array([0, 0, 0], dtype=int)],
            "p" : [np.array([1, 0, 0], dtype=int), np.array([0, 1, 0], dtype=int), np.array([0, 0, 1], dtype=int)],
            "d" : [np.array([2, 0, 0], dtype=int), np.array([1, 1, 0], dtype=int), np.array([1, 0, 1], dtype=int),
                   np.array([0, 2, 0], dtype=int), np.array([0, 1, 1], dtype=int), np.array([0, 0, 2], dtype=int)],
            "f" : [np.array([3, 0, 0], dtype=int), np.array([0, 3, 0], dtype=int), np.array([0, 0, 3], dtype=int),
                   np.array([2, 1, 0], dtype=int), np.array([2, 0, 1], dtype=int),
                   np.array([1, 2, 0], dtype=int), np.array([0, 2, 1], dtype=int),
                   np.array([1, 0, 2], dtype=int), np.array([0, 1, 2], dtype=int),
                   np.array([1, 1, 1], dtype=int)],
            "g" : [np.array([4, 0, 0], dtype=int), np.array([0, 4, 0], dtype=int), np.array([0, 0, 4], dtype=int),
                   np.array([3, 1, 0], dtype=int), np.array([3, 0, 1], dtype=int),
                   np.array([1, 3, 0], dtype=int), np.array([0, 3, 1], dtype=int),
                   np.array([1, 0, 3], dtype=int), np.array([0, 1, 3], dtype=int),
                   np.array([2, 1, 1], dtype=int), np.array([1, 2, 1], dtype=int), np.array([1, 1, 2], dtype=int),
                   np.array([2, 2, 0], dtype=int), np.array([2, 0, 2], dtype=int), np.array([0, 2, 2], dtype=int)]
    }
    return conv[nl_quant[1]]

@dataclass
class SlaterFunction:
    zeta: float
    nl_quant: str
    center: int
    xyz: np.array

@dataclass
class ContractedGaussianFunction:
    alphas: np.array
    coeffs: np.array
    l_vec: np.array
    center: int
    xyz: np.array
    normalized: bool

    def __post_init__(self):
        if not self.normalized:
            self._normalize_CGF()
            self.normalized = True

    def _normalize_CGF(self):
        L = np.sum(self.l_vec)

        for i in range(len(self.alphas)):
            self.coeffs[i] *= (2.0 / np.pi)**0.75 * 2.0**L * (self.alphas[i] ** ((2.0*L+3.0)/4.0)) / np.sqrt(factorial2(2*self.l_vec[0]-1)*factorial2(2*self.l_vec[1]-1)*factorial2(2*self.l_vec[2]-1))
        
        tmp = 0.0
        for ai, ci in zip(self.alphas, self.coeffs):
            for aj, cj in zip(self.alphas, self.coeffs):
                tmp += (ci * cj) / (ai+aj)**(L + 1.5)
        
        fac1 = 1.0 / np.sqrt(np.pi**1.5 * factorial2(2*self.l_vec[0]-1) * factorial2(2*self.l_vec[1]-1) * factorial2(2*self.l_vec[2]-1) * tmp / (2.0**L))
        self.coeffs *= fac1

class Basis:
    def __init__(self, mol: Molecule, basis_set: dict, normalized: bool):
        self.mol = mol
        self.basis = []

        if not self.mol.list_atom_types().issubset(set(basis_set.keys())):
            raise ValueError("Basis set does not cover all necessary atom types!")

        for idx, atom in enumerate(self.mol.geometry):
            for bf in basis_set[atom.symbol]:
                l_vecs = _Lvecs_from_nl(bf[0])
                for l_vec in l_vecs:
                    self.basis.append(ContractedGaussianFunction(np.array(bf[1]), np.array(bf[2]), l_vec, idx, atom.xyz, normalized))

        self.nbf = len(self.basis)

    @classmethod
    def from_slater(cls, mol : Molecule, dict_slater : dict, ng):
        return cls(mol, cls._slater_to_gauss_dict(dict_slater, ng), True) 
    
    @classmethod
    def from_file(cls, mol : Molecule, filename : str):
        return cls(mol, cls._gauss_dict_from_file(filename), False)

    def __call__(self, idx : int):
        return self.basis[idx]

    @staticmethod
    def _slater_to_gauss_dict(dict_s, ng):
        dict_g = {}
        for key, val in dict_s.items():
            bfs_for_atom = []
            for tup in val:
                coeffs, alphas = slater_expansion.slater_expansion(tup[1], tup[0], ng)
                bfs_for_atom.append((tup[0], alphas, coeffs))
            dict_g[key] = bfs_for_atom
        return dict_g
    
    @staticmethod
    def _gauss_dict_from_file(filename : str) -> dict:
        dict_g = {}
        with open(filename) as fl:
            lines = fl.read().splitlines()
            for l_num in range(1, len(lines)-1):
                if (lines[l_num-1] == '*' and lines[l_num+1] == '*'):
                    symb = lines[l_num].split()[0].title()
                    counter_bfs = np.arange(1, 6)
                    counter = 2
                    funcs = []
                    while (not (lines[l_num + counter] == '*')):
                        length, func_type =  int(lines[l_num + counter].split()[0]), lines[l_num + counter].split()[1].lower()
                        l_quant = slater_expansion._letter_to_l(func_type)
                        alphas = np.zeros(length)
                        coeffs = np.zeros(length)
                        for i in range(length):
                            alphas[i] = float(lines[l_num + counter + i + 1].split()[0].replace('D', 'E'))
                            coeffs[i] = float(lines[l_num + counter + i + 1].split()[1].replace('D', 'E')) #* Basis.norm2(l_quant, alphas[i])

                        nl_quant = str(counter_bfs[l_quant]) + str(func_type)
                        funcs.append((nl_quant, alphas, coeffs))            

                        counter_bfs[l_quant] += 1
                        counter += length + 1
                    dict_g[symb] = funcs
        return dict_g