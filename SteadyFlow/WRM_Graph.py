
from WRM import *
import matplotlib.pyplot as plt

M = 2

a = np.array([-0.168580646613524, 4.48591275435969e-17])

a = np.array([-1.33451781e-01, 5.00999217e-17,  3.87752892e-02])

a = np.array([0.435894451552019, -1.09126347681952,  0.600913856463807])

a =  np.array([ 0.56498618, 0.85399865, -9.80213127, 9.64475995])

#a = np.array([ 0.07395811, -0.05404956, -0.00071295])

#a = np.array([0, 0])

def phi_calc(a, px, py):
    return phi(a).subs(x, px).subs(y, py)

def dphiY_calc(a, px, py):
    return phi(a).diff(y).subs(x, px).subs(y, py)

if __name__ == "__main__":
    print("WRM Graph")
    
    Npoints = 21
    #pX = np.linspace(0, WIDTH, Npoints)
    pX = [i * WIDTH / (Npoints - 1) for i in range(Npoints)]
    #pY = np.linspace(dy, HEIGTH, Npoints)
    pY = [i * HEIGTH / (Npoints - 1) for i in range(Npoints)]
    
    print('---')
    ans = []
    for curY in pY:
        row = []
        for curX in pX:
            row.append(phi_calc(a, curX, curY))
        ans.append(row)
    
    print('---')
    X, Y = np.meshgrid(pX, pY)    
    #cs = plt.contour(X, Y, ans, 10)
    
    ##plt.yticks(range(-5, 5, 1))
    ##plt.xticks(range(-5, 5, 1))
    #plt.grid(True)
    #plt.clabel(cs)
    #plt.show()
    
    #plt.savefig("phi.png")
    
    ans = []
    for curX in pX:
        ans.append(dphiY_calc(a, curX, dy))
        
    cs = plt.plot(pX, ans)
    
    plt.grid(True)
    plt.show()
    plt.savefig("dphi.png")


    
    
