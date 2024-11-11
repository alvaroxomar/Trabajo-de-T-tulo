import numpy as np
from bip1 import algoritmo_rho_v, Rn

#t = algoritmo_rho_v(Rn)

def mod(x1,x2):
    return np.sqrt(x1**2+x2**2)

m = np.vectorize(mod)

g = m(np.ones((10,10))*4, np.ones((10,10)))