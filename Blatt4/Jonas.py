import numpy as np
import matplotlib.pyplot as plt

def o(t):
    N = len(t)
    re = 0
    im = 0
    for i in t:
        re += np.cos(i)/N
        im += np.sin(i)/N
    r = np.sqrt(re**2 + im**2)
    phi = np.arccos(re)
    return r, phi

def kuromoto(theta, args):
    """args = [w, K]"""
    r, phi = o(theta) #Ermittlung des Ordnungsparameters
    return args[0] + args[1]*r*np.sin(phi - theta)

def rk4(f, x, h, args):
    k1 = f(x, args)
    k2 = f(x + h/2*k1, args)
    k3 = f(x + h/2*k2, args)
    k4 = f(x + h*k3, args)
    ks = k1 + 2*k2 + 2*k3 + k4
    return x + h/6*(ks)

def plot(theta, h, T):
    plt.figure("phases")
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(np.sin(an), np.cos(an))
    plt.draw()
    plt.show(block=False)
    for i in range(1, int(T/h)):
        plt.title("t = " + str(round(h*i, 2)) + "s")
        t = plt.plot(np.sin(theta[i]), np.cos(theta[i]), "ro")
        r, phi = o(theta[i])
        t2 = plt.arrow(0, 0, r*np.sin(phi), r*np.cos(phi), head_width=0.05)
        plt.draw()
        plt.pause(0.01)
        t.pop(0).remove()
        t2.remove()

def main(N, K, T, h):
    w = 0.1*np.random.standard_cauchy(N) #Cauchy-Verteilung der Frequenzen
    theta = [np.random.random(N)*2*np.pi] #theta hat die für jedes theta[t] die 50 thetas
    pro = 0 #Prozentualer Fortschritt
    mess = 0 #Prozentualer Fortschritt der letzten Nachricht
    for i in range(1, int(T/h)): #For-Schleife durch alle Zeitschritte
        theta.append(rk4(kuromoto, theta[i-1], h, [w, K])) #neue Thetas bestimmen und an die Liste anhängen
        pro = 100*(i*h/T) #Fortschritt in Prozent ausrechnen
        if (mess + 1 < pro): #Entscheiden ob eine Neue Nachricht geschrieben werden soll
            print(str(round(pro, 0))+"% Fortschritt")
            mess = pro
    plot(theta, h, T)


main(50, 50, 0.5, 0.01)
