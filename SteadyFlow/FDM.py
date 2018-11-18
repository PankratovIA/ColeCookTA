import numpy as np

width = 1.0
heigth = 1.0

Nx = 10
Ny = 10

dx = width / Nx
dy = heigth / Ny

phi = np.zeros((Ny, Nx))

if __name__ == "__main__":
    print("FDM Transonic steady flow")
    print("dx = {0}".format(dx, "%f0.3"))
    print("dy = {0}".format(dy, "%f0.3"))
    print("phi = \n{0}".format(phi))
    print("BC >>>")
    for i in range (Ny):
        phi[i][0] = 1.0
    print("phi = \n{0}".format(phi))
    print("BC <<<")