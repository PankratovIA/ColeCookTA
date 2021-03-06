import numpy as np
#import matplotlib.pyplot as plt
#from sympy import Symbol, diff
#from sympy.abc import x


np.set_printoptions(formatter={'float': '{: 0.7f}'.format}, suppress=True)

WIDTH = 1.0 

Minf = 1.01
delta = 0.01

K = (1 - Minf**2)*(delta **(-2.0/3))
#K = 1
K = np.abs(K)
K = 0
print(K)
#exit(0)

gamma = 7.0/5

Phi_0 = 0.25
alpha = 45

alpha = np.radians(alpha)

dPhi_0 = np.tan(alpha)

Phi_a = 0.8

Nx = 6
dx = WIDTH / (Nx - 1)

dt = 0.5e-4
dt = 0.04*dx**2

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
#K = 1.0
#gamma = 1.4
#x_s = 0.61
#fi_x0 = 0.93
#assert fi_x0 > K/(gamma+1)
#fi_x1 = 2*K/(gamma+1) - fi_x0

#def fi(x):
  #if x <= x_s:
    #return fi_x0*x
  #else:
    #return fi_x0*x_s + fi_x1*(x-x_s)

#a = np.array([[0.25, -0.5, 1.0],
           #[-1.0,  1.0, 0.0],
           #[0.25,  0.5, 1.0]])
#b = np.array([0.0, fi_x0, fi(1.0)])
#k = np.linalg.solve(a, b)

#def fi_0(x):
  #return (k[0]*(x-0.5) + k[1])*(x-0.5) + k[2]

def parabola(x):
    """
        First approximation
    """
    A_ = (Phi_a - Phi_0 - dPhi_0*WIDTH) / (WIDTH ** 2.0)
    return A_ * x * x + dPhi_0 * x + Phi_0
    #return fi_0(x)

#[ 0.0000000  0.1246857  0.2493714  0.3740571  0.4987429  0.5408763
  #0.5352380  0.5296000]


if __name__ == "__main__":
    print("FDM 1D")
    
    print("K =", K)
    print("gamma =", gamma)
    
    x = np.linspace(0, WIDTH, Nx)
    phi_prev = phi = np.array(list(map(parabola, x)))
    
    print("x =", x)
    print("phi =", phi)
    
    #phi[1] = phi[0] + dPhi_0 * dx
    #phi[-1] = Phi_a
    
#     Result phi = 
# [ 0.2500000  0.3836735  0.5084070  0.6011152  0.6849887  0.7394815
#   0.7871739  0.8000000]
    
    cnt = 0
    while cnt <150000:
        cnt += 1
        print("Iteration {0}".format(cnt))
        
        cur = np.zeros(Nx)
        
        mx = 0
        
        A = np.zeros((Nx, Nx))
        b = np.zeros((Nx, 1))
        A[0][0] = 1.0
        b[0] = 0
        
        #A[1][2] = 1.0 / (2 * dx)
        #A[1][0] = -1.0 / (2 * dx)
        
        A[1][1] = 1.0
        b[1] = 0
        
        A[-1][-1] = 1.0
        b[-1] = 0
        
        
        for j in range(2, Nx-1):
            # Blinkov Yu.A., p. 157
            phi_xx = phi[j + 1] - 2 * phi[j] + phi[j - 1]
            phi_xx2 = (phi[j] - 2 * phi[j - 1] + phi[j - 2])
            
            b[j] = -phi_xx * (phi_xx2 * (K - (0.5 / dx)*(gamma+1) *\
            (phi[j+1] - phi[j] + phi[j-1] - phi[j-2]))*dt/dx\
             - (-phi[j-1]+phi_prev[j-1])*dx/(2*dt))
            
            A[j][j] += -phi_xx
            A[j][j-2] += phi_xx
            A[j][j-1] += -phi_xx * dx / (2*dt)

            b[j] += -phi_xx2 * (phi_xx * (K - (0.5 / dx)*(gamma+1) *\
            (phi[j+1] - phi[j] + phi[j-1] - phi[j-2]))*dt/dx\
            - (-phi[j]+phi_prev[j])*dx/(2*dt))
            
            A[j][j+1] += -phi_xx2
            A[j][j-1] += phi_xx2
            A[j][j] += -phi_xx2 * dx / (2*dt)           
        
        #print("A =", A)
        #print("b =", b)
        
        cur = np.linalg.solve(A, b)
        
        cur = cur.reshape(Nx)
        
        #print("cur =", cur)
        print("phi =", phi)
        #print("cur - phi =", cur - phi)
        
        mx = max(map(abs, cur))
        
        phi_prev = phi.copy()
        phi = phi.copy() + cur.copy()

        print("mx = {0}".format(mx, "%.3f"))
        print("phi =", phi)
        if mx < eps:
            break
    print("iter = ", cnt)
    print("Result phi = \n{0}".format(phi))

    
    
