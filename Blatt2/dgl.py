import numpy as np
import matplotlib.pyplot as plt

def f(N, r, K):
    return N*r*(1-(N/K))

def pltN(r, h, K, N0, heun=False):
    N = []
    N.append(N0)
    for i in range(int(5/h)):
        N1 = N[i] + h*f(N[i], r, K)
        if(heun):
            N.append( N[i] + 0.5*h*(f(N[i], r, K) + f(N1, r, K)) )
        else:
            N.append(N1)
    label = "h="+str(h)+", K="+str(K)+", N0="+str(N0)+", Heun="+str(heun)
    plt.plot(np.arange(0, 5+h, h), N, label=(label))

def Aufgabe1():
    pltN(20, 0.01, 100, 5)
    pltN(20, 0.01, 100, 20)
    pltN(20, 0.01, 200, 5)
    pltN(20, 0.1, 100, 5)

def Aufgabe2():
    pltN(20, 0.1, 100, 5, True)

Aufgabe1()
Aufgabe2()
plt.legend()
plt.show()
