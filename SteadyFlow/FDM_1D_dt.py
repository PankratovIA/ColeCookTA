import numpy as np
#import matplotlib.pyplot as plt
#from sympy import Symbol, diff
#from sympy.abc import x


np.set_printoptions(formatter={'float': '{: 0.5f}'.format}, suppress=True)

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

Nx = 16
dx = WIDTH / (Nx - 1)

dt = 1e-2

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
    
    phi[1] = phi[0] + dPhi_0 * dx
    phi[-1] = Phi_a
    
    print("dphi = ", [(x[1] - x[0]) / dx for x in zip(phi, phi[1:])])
    
    cnt = 0
    while cnt < 2000:
        cnt += 1
        print("Iteration {0}".format(cnt))
        
        cur = np.zeros(Nx)
        
        mx = 0
        for idx in range(2, Nx-1):
            # phi_x = (phi[idx + 1] - phi[idx - 1]) / (2 * dx)
            # phi_xx = (phi[idx + 1] - 2.0 * phi[idx] + phi[idx - 1]) / (dx * dx)
            # cur[idx] = phi[idx] + dt * phi_xx * (K - (gamma + 1) * phi_x)
            
            # Blinkov Yu.A., p. 157
            first = 0
            tmp = phi[idx + 1] - 2 * phi[idx] + phi[idx - 1]
            tmp2 = (phi[idx] - 2 * phi[idx - 1] + phi[idx - 2])
            
            first += tmp * tmp2 * K
            first -= tmp * tmp2 * \
            0.5 * (gamma + 1) *(phi[idx + 1] - 1 * phi[idx] + phi[idx - 1] - phi[idx - 2])
            
            mul = -tmp * tmp2 * 0.5 * (gamma + 1)
            
            second = 0 
            tmp = phi[idx] - 2 * phi[idx - 1] + phi[idx - 2]
            tmp2 = (phi[idx + 1] - 2 * phi[idx] + phi[idx - 1])
            
            second += tmp * tmp2 * K
            second -= tmp * tmp2 * \
            0.5 * (gamma + 1) *(phi[idx + 1] - 1 * phi[idx] + phi[idx - 1] - phi[idx - 2])
            
            mul += -tmp * tmp2 * 0.5 * (gamma + 1)
            
            
            cur[idx] = phi[idx] + dt * (first + second) / dx
            #cur[idx] = (first + second) * mul
            
            mx = max(mx, abs(cur[idx] - phi[idx]))
        
        phi[2: Nx-1] = cur[2: Nx-1]

        print("mx = {0}".format(mx, "%.3f"))
        print("phi =", phi)
        if mx < eps:
            break
    print("iter = ", cnt)
    print("Result phi = \n{0}".format(phi))

    
    
