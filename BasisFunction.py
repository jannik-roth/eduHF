import numpy as np
import slater_expansion
from dataclasses import dataclass
import geometry

class BasisFunction:
    def __init__(self, zeta: float, xyz, nl_quant: str, ng: int):
        self.zeta = zeta
        self.xyz = np.asarray(xyz)
        self.nl_quant = nl_quant
        self.ng = ng

        self.coeffs = np.zeros(ng)
        self.alphas = np.zeros(ng)
        slater_expansion.slater_exp(self.alphas, self.coeffs, self.ng, self.nl_quant)

@staticmethod
def _Lvecs_from_nl(nl_quant: str):
    conv = {"s" : [np.array([0, 0, 0], dtype=int)],
            "p" : [np.array([1, 0, 0], dtype=int), np.array([0, 1, 0], dtype=int), np.array([0, 0, 1], dtype=int)],
            "d" : [np.array([2, 0, 0], dtype=int), np.array([0, 2, 0], dtype=int), np.array([0, 0, 2], dtype=int),
                   np.array([1, 1, 0], dtype=int), np.array([0, 1, 1], dtype=int), np.array([1, 0, 1], dtype=int)],
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

@staticmethod 
def _ls_from_nl(nl_quant : str):
    conv = {"s" : 1,
            "p" : 3,
            "d" : 6,
            "f" : 10,
            "g" : 15}
    return conv[nl_quant[1]]

@dataclass
class SlaterFunction:
    zeta: float
    nl_quant: str
    center: int
    xyz: np.array

@dataclass
class ContractedGaussianFunction:
    ng: int
    alphas: np.array
    coeffs: np.array
    l_vec: np.array
    center: int
    xyz: np.array

def Slater2CGaussians(slater: SlaterFunction, ng: int):
    cgaussians = []

    l_vecs = _Lvecs_from_nl(slater.nl_quant)
    coeffs = np.zeros(ng)
    alphas = np.zeros(ng)
    slater_expansion.slater_exp(alphas, coeffs, ng, slater.nl_quant)

    for l_vec in l_vecs:
        cgaussians.append(ContractedGaussianFunction(ng, alphas, coeffs, l_vec, slater.center, slater.xyz))
    
    return cgaussians

def Slaters2CGaussians(slaters: list[SlaterFunction], ng: int):
    cgaussians = []
    for slater in slaters:
        cgaussians.append(Slater2CGaussians(slater, ng))
    return cgaussians


class BasisSet:
    def __init__(self, basis_dict= dict[ str: list[tuple[str, float]]]):
        self.basis = {}
        for at, funcs in basis_dict.items():
            self.basis[at] = funcs

    def __call__(self, at):
        return self.basis[at]
    
    def available_atom_types(self):
        return set(self.basis.keys())
    
class Basis:
    def __init__(self, mol: geometry.Molecule, basis_set: BasisSet, ng: int):
        self.mol = mol
        self.basis_set = basis_set
        self.ng = ng
        self._check_compatible()

        self.basis_slater = self._build_basis_slater()
        print(self.basis_slater)
        self.basis_gauss = self._build_basis_gauss()
        print(self.basis_gauss)

    def _check_compatible(self):
        print(self.mol.list_atom_types())
        print(self.basis_set.available_atom_types())
        if not self.mol.list_atom_types().issubset(self.basis_set.available_atom_types()):
            raise ValueError("BasisSet does not cover all necessary atom types!!!")
    
    def _build_basis_slater(self) -> list[SlaterFunction]:
        basis_sl = []
        for idx, atom in enumerate(self.mol.geometry):
            bfs = self.basis_set(atom.symbol)
            for (nl, zeta) in bfs:
                basis_sl.append(SlaterFunction(zeta, nl, idx, atom.xyz))
        return basis_sl
    
    def nof_slater_funcs(self):
        return len(self.basis_slater)
    
    def _build_basis_gauss(self) -> list[ContractedGaussianFunction]:
        basis = []
        for sl_func in self.basis_slater:
            basis.extend(Slater2CGaussians(sl_func, self.ng))
        return basis
    
    @staticmethod
    def Slater2CGaussians(slater: SlaterFunction, ng: int):
        cgaussians = []

        l_vecs = _Lvecs_from_nl(slater.nl_quant)
        coeffs = np.zeros(ng)
        alphas = np.zeros(ng)
        slater_expansion.slater_exp(alphas, coeffs, ng, slater.nl_quant)
        for l_vec in l_vecs:
            cgaussians.append(ContractedGaussianFunction(ng, alphas, coeffs, l_vec, slater.center, slater.xyz))
        return cgaussians
    
    def nof_gauss_funcs(self):
        return len(self.basis_gauss)


