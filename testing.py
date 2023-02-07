import geometry
import BasisFunction
import slater_expansion
import numpy as np

h1 = [['H', [1, 2, 3]],
    ]

test = geometry.geometry(h1)

test = BasisFunction.BasisFunction(1.0, [1.3, 2.5, 1.7], "5d", 3)

s = 3
a = np.ones(s)
slater_expansion.testing(a, s)
print(a)