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
    nofelectrons: int = field(init=False)

    def __post_init__(self):
        self.nofatoms = len(self.geometry)
        self.nofelectrons = 0
        for atom in self.geometry:
            self.nofelectrons += atom.charge
        self.nofelectrons -= self.charge

    def list_atom_types(self):
        return set([atom.symbol for atom in self.geometry])
