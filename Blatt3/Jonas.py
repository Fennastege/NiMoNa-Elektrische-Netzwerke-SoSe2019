import numpy as np
import matplotlib.pyplot as plt

tmax = 3
h = 0.01

def f(p1, p2):
    return [p1*(2-0.02*p2), -p2*(0.8-0.0002*p1)]

def rk4(p1, p2, h):
    k = []
    k.append(f(p1, p2))
    k.append(f(p1 + h/2*k[0][0], p2 + h/2*k[0][1]))
    k.append(f(p1 + h/2*k[1][0], p2 + h/2*k[1][1]))
    k.append(f(p1 + h*k[2][0], p2 + h*k[2][1]))
    ks = np.add(np.add(np.add(k[0], k[1]), k[2]), k[3])
    return p1 + h/6*(ks[0]), p2 + h/6*(ks[1])

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
    #plt.plot(np.arange(0, tmax, h), p1)
    #plt.plot(np.arange(0, tmax, h), p2)
    plt.plot(np.arange(0, tmax, h), np.subtract(p1, p2))
    #plt.plot(p1, p2)
    plt.show()

rb(200, 0.0005, 500, 400)