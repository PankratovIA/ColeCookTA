WIDTH = 1.0 
HEIGTH = 1.0

Vinf = 1.0

LEFT = WIDTH / 3.0
RIGHT = 2.0 * WIDTH / 3.0


EPS = 1e-6

M = 3

Ny = 100

dy = HEIGTH / Ny

from sympy.abc import x, y
from sympy import pi, sin, cos, integrate

import numpy as np
from random import seed, random

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

def psi():
    A_ = (Phi_a - Phi_0 - dPhi_0*WIDTH) / (WIDTH ** 2.0)
    return A_ * x * x + dPhi_0 * x + Phi_0

def N(m):
    return (WIDTH - x) * (x ** m)
    
def phi(a):
    #print(a, len(a))
    #assert(len(a) == M)
    ans = psi()
    for idx, cur in enumerate(a):
        ans += cur * N(idx+1)
    return ans
    
   
if __name__ == "__main__":
    print("WRM_1D Transonic steady flow")
    print(N(3))
    print(N(3).diff(x).diff(x))
    
    res = integrate(sin(x)*sin(y), (x, 0, WIDTH), (y, 0, HEIGTH))
    print(res)
    
    print(psi())
    
    seed()
    tmp = [random() for _ in range(M)]
    print(tmp)
    print(phi(tmp))
    print(phi(tmp).diff(x))
    #exit(0)
    a = np.zeros(M)
    print(phi(a))
    
    #exit(0)
    
    a = np.zeros(M)
    
    cnt = 0
    while cnt < 50:
        cnt += 1
        print("Iteration = ", cnt)
        A = np.zeros((M, M))
        b = np.zeros((M, 1))
        for s in range(1, M+1):
            #print("s = ", s)
            tmp = K * N(s).diff(x) + (gamma+1) * N(s) * phi(a).diff(x).diff(x)
            b[s-1] = -integrate(psi().diff(x) * tmp, (x, 0, WIDTH))
            for m in range(1, M+1):
                #print("m = ", m)
                A[s-1][m-1] = integrate(N(m).diff(x) * tmp, (x, 0, WIDTH))
            
        print("A = \n", A)
        print("det(K) = \n", np.linalg.det(A))
        print("b = \n", b)
    
        a_prev = a.reshape((M, 1))
        print("a_prev = ", a_prev)
        a = np.linalg.solve(A, b)
        print("a = ", a)
    
        diff = np.linalg.norm(a - a_prev)
        print("diff =", diff)        
        
        a = a.reshape(M)
        
        if diff < EPS:
            break
        
        
        
    print("a =", a)
    
    ans = phi(a)
    
    print("ans =", ans)
    
    
    #X = np.linalg.
