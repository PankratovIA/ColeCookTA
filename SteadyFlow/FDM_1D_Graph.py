from FDM_1D import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("FDM 1D Graph")
    
    N = 101
    x = np.linspace(0, WIDTH, N)
    print("x = ", x)
    y = np.array(list(map(analytical_solution, x)))
    par = np.array(list(map(parabola, x)))
    plt.plot(x, y)
    plt.plot(x, par, '--')
    
    phi = [ 0.250,  0.300,  0.375,  0.450,  0.525,  0.600,  0.675,  0.749,  0.824,  0.899,
  0.890,  0.881,  0.872,  0.863,  0.854,  0.845,  0.836,  0.827,  0.818,  0.809, 0.800]
  
    phi = [ 0.25000, 0.35000, 0.42890, 0.50780, 0.56909, 0.63037, 0.67438, 0.71839,
  0.74770, 0.77702, 0.80000]
  
     
    x = np.linspace(0, WIDTH, len(phi))
    print("x = ", x)
    plt.plot(x, phi, 'o')
    plt.grid()
    
    plt.ylim()
    print(plt.ylim())
    bottom, top = plt.ylim()  # return the current ylim
    plt.ylim((0.95*bottom, 1.05*top))   # set the ylim to bottom, top
    plt.show()
    
    
