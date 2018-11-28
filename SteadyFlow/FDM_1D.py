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

Nx = 51
dx = WIDTH / (Nx - 1)

eps = 1e-9

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

def createMatrix(phi):
    A = np.zeros((Nx, Nx))
    b = np.zeros((Nx, 1))
    # phi(0) = Phi_0
    A[0][0] = 1.0
    b[0] = Phi_0
    
    # dphi(0) = dPhi_0
    A[1][1] = 1.0 / dx
    A[1][0] = -1.0 / dx
    b[1] = dPhi_0
    
    for row in range(2, Nx-1):
        dphi = (phi[row] - phi[row - 1]) / dx
        if dphi < dPhi_0:
            # K - (gamma + 1) * phi_x = 0
            #print("phi_x")
            A[row][row] = -1.0 / dx
            A[row][row+1] = 1.0 / dx
            b[row] = K / (gamma + 1)
        else:
            # phi_xx = 0
            #print("phi_xx")
            A[row][row-1] = 1.0 / (dx * dx)
            A[row][row] = -2.0 / (dx * dx)
            A[row][row+1] = 1.0 / (dx * dx)
            # b[row] = 0

    A[Nx-1][Nx-1] = 1.0
    b[Nx-1] = Phi_a
    return A, b


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
    
    print("dphi = ", [(x[1] - x[0]) / dx for x in zip(phi, phi[1:])])
    
    cnt = 0
    while cnt < 200:
        cnt += 1
        print("Iteration {0}".format(cnt))
        A, b = createMatrix(phi)
        #print("A =\n", A)
        #print("b =\n", b)
    
        cur = np.linalg.solve(A, b)
        #print("cur = \n{0}".format(cur))
        
        mx = 0
        for idx in range(Nx):
            diff = abs(phi[idx] - cur[idx])
            mx = max(mx, diff)
            phi[idx] = cur[idx]
        print("mx = {0}".format(mx, "%.3f"))
        if mx < eps:
            break
    print("iter = ", cnt)
    print("Result phi = \n{0}".format(phi))

    
    
