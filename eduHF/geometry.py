import numpy as np
from dataclasses import dataclass, field

@dataclass
class Atom:
    symbol: str
    xyz : np.array
    charge: float = field(init=False)

    def __post_init__(self):
        self.charge = self._asymb_to_charge(self.symbol)

    @staticmethod
    def _asymb_to_charge(asymb: str) -> int:
        conv = {'H'  : 1,
                'He' : 2,
                'Li' : 3,
                'Be' : 4,
                'B'  : 5,
                'C'  : 6,
                'N'  : 7,
                'O'  : 8,
                'F'  : 9,
                'Ne' : 10}
        return conv[asymb]

@dataclass
class Molecule:
    geometry: list[Atom] = field(default_factory=list)
    charge: int = 0
    nofatoms: int = field(init=False)
    noe: int = field(init=False)
    core_pot: float = field(init=False)

    def __post_init__(self):
        self.nofatoms = len(self.geometry)
        self.noe = 0
        for atom in self.geometry:
            self.noe += atom.charge
        self.noe -= self.charge
        self.core_pot = self.core_potential()

    def list_atom_types(self):
        return set([atom.symbol for atom in self.geometry])
    
    def core_potential(self):
        res = 0.0
        for a in range(self.nofatoms):
            for b in range(a+1, self.nofatoms):
                res += self.geometry[a].charge * self.geometry[b].charge * 1.0/(self._distance(a, b))
        return res
    
    def core_potential_der(self, center, dim):
        res = 0.0
        for b in range(self.nofatoms):
            if not (center == b):
                if (dim == 0):
                    res += self.geometry[b].charge * (self.geometry[b].xyz[0] - self.geometry[center].xyz[0]) / self._distance(center, b)**3.0
                elif (dim == 1):
                    res += self.geometry[b].charge * (self.geometry[b].xyz[1] - self.geometry[center].xyz[1]) / self._distance(center, b)**3.0
                elif (dim == 2):
                    res += self.geometry[b].charge * (self.geometry[b].xyz[2] - self.geometry[center].xyz[2]) / self._distance(center, b)**3.0
        res *= self.geometry[center].charge
        return res
    def _distance(self, a : int, b : int):
        xyza = self.geometry[a].xyz
        xyzb = self.geometry[b].xyz
        return np.sqrt(np.sum((xyza - xyzb)**2))
