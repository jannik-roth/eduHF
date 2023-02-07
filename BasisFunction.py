import numpy as np
import slater_expansion

class BasisFunction:
    def __init__(self, zeta: float, xyz, lm_quant: str, ng: int):
        self.zeta = zeta
        self.xyz = np.asarray(xyz)
        self.lm_quant = lm_quant
        self.ng = int(ng)

        self.coeffs = np.zeros(ng)
        self.alphas = np.zeros(ng)
        slater_expansion.slater_exp(self.alphas, self.coeffs, self.ng, self.lm_quant)
        print(self.alphas)
        print(self.coeffs)

