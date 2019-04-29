import numpy as np
import matplotlib.pyplot as plt

def f(p1, p2):
    return p1*(2-0.02*p2), -p2*(0.8-0.0002*p1)

def rk4(p1, p2, h):
    k1, l1 = f(p1, p2)
    k2, l2 = f(p1 + h/2*k1, p2 + h/2*l1)
    k3, l3 = f(p1 + h/2*k2, p2 + h/2*l2)
    k4, l4 = f(p1 + h*k3, p2 + h*l3)
    ks = k1 + 2*k2 + 2*k3 + k4
    ls = l1 + 2*l2 + 2*l3 + l4
    return p1 + h/6*(ks), p2 + h/6*(ls)

def rb(tmax, h, p10, p20):
    p1 = [0]*(int(tmax/h))
    p2 = [0]*(int(tmax/h))
    p1[0] = p10
    p2[0] = p20
    pro = 0
    mess = 0
    for i in range(1, int(tmax/h)):
        p1[i], p2[i] = rk4(p1[i-1], p2[i-1], h)
        pro = 100*(i*h/tmax)
        if (mess + 10 < pro):
            print(str(round(pro, 0))+"% Fortschritt")
            mess = pro
    plt.plot(np.arange(0, tmax, h), p1, label = "p1")
    plt.plot(np.arange(0, tmax, h), p2, label = "p2")
    plt.plot(np.arange(0, tmax, h), np.subtract(p1, p2), label = "p1-p2")
    plt.legend()
    plt.show()
    plt.plot(p1, p2)
    plt.show()

rb(500, 0.025, 1000, 20)