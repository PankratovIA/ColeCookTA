import numpy as np

WIDTH = 1.0
HEIGTH = 1.0

Nx = 3
Ny = 3
SIZE = Nx * Ny
       

dx = WIDTH / Nx
dy = HEIGTH / Ny

eps = 1e-9

phi = np.zeros((Ny+2, Nx+2))

def cellNum(row, col):
    assert 1<=row<=Ny
    assert 1<=col<=Nx
    return (row - 1) * Ny + col - 1
    
def createMatrix(Nx, Ny):
    #size = Nx * Ny
    print("SIZE = ", SIZE)
    A = 2* np.eye(SIZE)
    A = np.zeros((SIZE, SIZE))
    for row in range(1, Ny):
        for col in range(1, Nx):
            if not onBoundary(row, col):
                num = cellNum(row, col)
                A[num][num] = 10
    return A

def createBC(A):
    b = np.zeros((SIZE, 1))
    for col in range(1, Nx+1):
        # phi(x, 0) = 0
        num = cellNum(1, col)
        A[num] = np.zeros((1, SIZE))
        A[num][num] = 1.0
        b[num] = 0.0
            
        # phi(x, b) = 0
        num = cellNum(Nx, col)
        A[num] = np.zeros((1, SIZE))
        A[num][num] = 1.0
        b[num] = 0.0

    for row in range(1, Ny+1):
        # phi(0, y) = 1
        num = cellNum(row, 1)
        A[num] = np.zeros((1, SIZE))
        A[num][num] = 1.0
        b[num] = 1.0
           
        # phi(a, y) = 0
        num = cellNum(row, Nx)
        A[num] = np.zeros((1, SIZE))
        A[num][num] = 1.0
        b[num] = 0.0
    return (A, b)
    
def onBoundary(row, col):
    return (min(row, col) == 1) | (row == Ny) | (col == Nx)

if __name__ == "__main__":
    print("FDM Transonic steady flow")
    print("dx = {0}".format(dx, "%f0.3"))
    print("dy = {0}".format(dy, "%f0.3"))
    print("phi = \n{0}".format(phi))
    print("BC >>>")
    for i in range(Ny+2):
        phi[i][0] = 1.0
        phi[i][1] = 1.0
    print("phi = \n{0}".format(phi))
    print("BC <<<")
    
    iter = 0
    while 1:
        iter += 1
        print("Iteration {0}".format(iter))
        A = createMatrix(Nx, Ny)
        print("A = \n{0}".format(A))
        print("det( A ) = {0:0.3f}".format(np.linalg.det(A)))
        
        A, b = createBC(A)
        
        print("A = \n{0}".format(A))
        print("det( A ) = {0:0.3f}".format(np.linalg.det(A)))
        print("b = \n{0}".format(b))
    
        x = np.linalg.solve(A, b)
        print("x = \n{0}".format(x))
        
        mx = 0
        for row in range(1, Ny+1):
            for col in range(1, Nx+1):
                num = cellNum(row, col)
                diff = abs(phi[row][col] - x[num])
                mx = max(mx, diff)
                phi[row][col] = x[num]
        print("mx = {0}".format(mx, "%.3f"))
        if mx < eps:
            break
    print("Result phi = \n{0}".format(phi))
