import numpy as np

WIDTH = 1.0
HEIGTH = 1.0

Nx = 3
Ny = 3
SIZE = Nx * Ny
       

dx = WIDTH / Nx
dy = HEIGTH / Ny

eps = 1e-9

phi = np.zeros((Ny, Nx))

def cellNum(row, col):
    assert 0<=row<Ny
    assert 0<=col<Nx
    return row * Ny + col
    
def createMatrix(Nx, Ny):
    size = Nx * Ny
    A = 2* np.eye(size)
    return A

if __name__ == "__main__":
    print("FDM Transonic steady flow")
    print("dx = {0}".format(dx, "%f0.3"))
    print("dy = {0}".format(dy, "%f0.3"))
    print("phi = \n{0}".format(phi))
    print("BC >>>")
    for i in range(Ny):
        phi[i][0] = 1.0
    print("phi = \n{0}".format(phi))
    print("BC <<<")
    
    iter = 0
    while 1:
        iter += 1
        print("Iteration {0}".format(iter))
        A = createMatrix(Nx, Ny)
        print("A = \n{0}".format(A))
        print("det( A ) = {0:0.3f}".format(np.linalg.det(A)))
        
        b = np.zeros((SIZE, 1))
        for row in range(Ny):
            num = cellNum(row, 0)
            A[num] = np.zeros((1, SIZE))
            A[num][num] = 1.0
            b[num] = 1.0
        
        print("A = \n{0}".format(A))
        print("det( A ) = {0:0.3f}".format(np.linalg.det(A)))
        print("b = \n{0}".format(b))
    
        x = np.linalg.solve(A, b)
        print("x = \n{0}".format(x))
        
        mx = 0
        for row in range(Ny):
            for col in range(Nx):
                num = cellNum(row, col)
                diff = abs(phi[row][col] - x[num])
                mx = max(mx, diff)
                phi[row][col] = x[num]
        print("mx = {0}".format(mx, "%.3f"))
        if mx < eps:
            break
    print("Result phi = \n{0}".format(phi))
