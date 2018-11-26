WIDTH = 1.0 
HEIGTH = 1.0

Vinf = 1.0

LEFT = WIDTH / 3.0
RIGHT = 2.0 * WIDTH / 3.0


EPS = 1e-9

M = 3

Ny = 100

dy = HEIGTH / Ny

from sympy.abc import x, y
from sympy import pi, sin, cos, integrate

import numpy as np
from random import seed, random

def psi():
    return Vinf * x

def N(m):
    return cos(pi * m * x / WIDTH) * cos(pi * m * y / HEIGTH)
    
def phi(a):
    #print(a, len(a))
    #assert(len(a) == M)
    ans = psi()
    for idx, cur in enumerate(a):
        ans += cur * N(idx+1)
    return ans
    
def airfoil(left, right):
    assert(0 < left < right < WIDTH)
    coord = [left, 0.5 * (left + right), right]
    A = np.zeros((3, 3))
    for row in range(3):
        A[row] = np.array([coord[row] ** 2, coord[row], 1])
    B = np.array([[0], [3*dy/2], [0]])
    #print(A)
    # A *[a, b, c]^T = B
    #print("x =", x, type(x))
    a, b, c = np.linalg.solve(A, B)
    #res = np.linalg.solve(A, B)
    #print(a, b, c)
    #print(res)
    ans = a[0] * x**2 + b[0] * x + c[0]
    #print(ans, type(ans))
    #print(f[0], type(ans))
    return ans
    
if __name__ == "__main__":
    print("WRM Transonic steady flow")
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
    
    print(airfoil(0.3, 0.6))
    #exit(0)
    
    airfoil_func = airfoil(LEFT, RIGHT).diff(x)
    print("airfoil_func", airfoil_func)
    
    cnt = 0
    while cnt < 20:
        cnt += 1
        print("Iteration = ", cnt)
        K = np.zeros((M, M))
        F = np.zeros((M, 1))
        for s in range(1, M+1):
            print("s = ", s)
            func = airfoil_func * N(s)
            func = func.subs(y, dy)
            #tmp = integrate(func, (x, LEFT, RIGHT))
            #print("tmp =", tmp)
            F[s-1] = -integrate(func, (x, LEFT, RIGHT))
            for m in range(1, M+1):
                print("m = ", m)
                func = phi(a) * N(m).diff(x).diff(x) - N(m).diff(y).diff(y)
                func *= N(s)
                #print(func)
                K[s-1][m-1] = integrate(func, (x, 0, WIDTH), (y, dy, HEIGTH))
            
                func = N(m).diff(y) * N(s)
                #print(func)
                func = func.subs(y, dy)
                #print(func)
                K[s-1][m-1] -= integrate(func, (x, 0, WIDTH))
                #print("K = \n", K)
                #print("F = \n", F)
    
        print("K = \n", K)
        print("det(K) = \n", np.linalg.det(K))
        #print("F = \n", F)
    
        a_prev = a.reshape((M, 1))
        print("a_prev = ", a_prev)
        a = np.linalg.solve(K, F)
        print("a = ", a)
    
        diff = np.linalg.norm(a - a_prev)
        print("diff =", diff)        
        
        a = a.reshape(M)
        
        if diff < EPS:
            break
        
        
        
    print("a =", a)
    
    ans = phi(a)
    
    print("ans =", ans)
    
    BC = ans.diff(y)

    print("d ans / dy =", BC)
    
    print("d ans / dy (y=H) =", BC.subs(y, HEIGTH))  
    
    derivY_dy = BC.subs(y, dy)
    
    print("d ans / dy (y=dy) =", derivY_dy)
    
    #X = np.linalg.
