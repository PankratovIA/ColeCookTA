from FDM_1D import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("FDM 1D Graph")
    
    N = 1001
    x = np.linspace(0, WIDTH, N)
    print("x = ", x)
    y = np.array(list(map(analytical_solution, x)))
    par = np.array(list(map(parabola, x)))
    plt.plot(x, y)
    plt.plot(x, par)
    plt.ylim()
    print(plt.ylim())
    bottom, top = plt.ylim()  # return the current ylim
    #plt.ylim((0.95*bottom, 1.05*top))   # set the ylim to bottom, top
    plt.show()
    
    