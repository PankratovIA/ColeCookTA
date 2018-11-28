import numpy as np
#import matplotlib.pyplot as plt
#from sympy import Symbol, diff
#from sympy.abc import x


np.set_printoptions(formatter={'float': '{: 0.3f}'.format}, suppress=True)

WIDTH = 1.0 

Minf = 1.01
delta = 0.01

K = (1 - Minf**2)*(delta **(-2.0/3))

gamma = 7.0/5

Phi_0 = 0.25
alpha = 45

alpha = np.radians(alpha)

dPhi_0 = np.tan(alpha)

Phi_a = 0.8

Nx = 11
dx = WIDTH / (Nx - 1)

def analytical_solution(x):
    xu = Phi_a - Phi_0 - K * WIDTH / (gamma + 1)
    xu /= (dPhi_0 - K/(gamma+1))
    assert(0<=xu<=WIDTH)
    #print("xu =", xu)
    
    if x <= xu:
        ans = dPhi_0 * x + Phi_0
    else:
        ans = Phi_a + K*(x - WIDTH) / (gamma+1)
    
    return ans

def parabola(x):
    """
        First approximation
    """
    A_ = (Phi_a - Phi_0 - dPhi_0*WIDTH) / (WIDTH ** 2.0)
    return A_ * x * x + dPhi_0 * x + Phi_0

if __name__ == "__main__":
    print("FDM 1D")
    
    print("K =", K)
    print("gamma =", gamma)
    
    p = WIDTH / 2
    print(analytical_solution(p))
    
    x = np.linspace(0, WIDTH, Nx)
    phi = np.array(list(map(parabola, x)))
    
    print("x =", x)
    print("phi =", phi)
    
    
