
from WRM import *
import matplotlib.pyplot as plt

M = 2

a = np.array([-0.168580646613524, 4.48591275435969e-17])

a = np.array([-1.33451781e-01, 5.00999217e-17,  3.87752892e-02])

#a = np.array([ 0.07395811, -0.05404956, -0.00071295])

#a = np.array([0, 0])

def phi_calc(a, px, py):
    return phi(a).subs(x, px).subs(y, py)

def dphiY_calc(a, px, py):
    return phi(a).diff(y).subs(x, px).subs(y, py)

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
    
    
    ans = []
    for curX in pX:
        ans.append(dphiY_calc(a, curX, dy))
        
    cs = plt.plot(pX, ans)
    
    plt.grid(True)
    plt.show()


    
    
