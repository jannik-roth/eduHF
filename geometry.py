import numpy as np

class geometry:
    def __init__(self, inp_geom: list[str, list]):
        self.number_of_atoms = len(inp_geom)
        self.xyz = np.zeros((3, len(inp_geom)))
        self.atom_types = np.empty(len(inp_geom), dtype=object)
        self.charges = np.zeros(len(inp_geom))

        for idx, atom in enumerate(inp_geom):
            self.atom_types[idx] = atom[0]
            self.xyz[:,idx] = np.asarray(atom[1], float)
            self.charges[idx] = self._asymb_to_charge(atom[0])
    
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


        