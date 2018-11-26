
from WRM import *
import matplotlib.pyplot as plt

M = 2

a = np.array([-0.00809222859605956, -1.0144816987212e-5])

a = np.array([-0.00809179468346368, 9.17361370121491e-6, 0.00213799041608868])

#a = np.array([0, 0])

def phi_calc(a, px, py):
    return phi(a).subs(x, px).subs(y, py)

if __name__ == "__main__":
    print("WRM Graph")
    
    Npoints = 21
    pX = np.linspace(0, WIDTH, Npoints)
    pY = np.linspace(dy, HEIGTH, Npoints)
    
    ans = []
    for curY in pY:
        row = []
        for curX in pX:
            row.append(phi_calc(a, curX, curY))
        ans.append(row)
    
    X, Y = np.meshgrid(pX, pY)    
    cs = plt.contour(X, Y, ans, 10)
    
    #plt.yticks(range(-5, 5, 1))
    #plt.xticks(range(-5, 5, 1))
    plt.grid(True)
    plt.clabel(cs)
    plt.show()
