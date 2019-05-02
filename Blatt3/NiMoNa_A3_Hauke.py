# -*- coding: utf-8 -*-
# %% Preamble NiMoNa_A3 Hauke Bents
import numpy as np
import matplotlib.pyplot as plt
    
# %% Entwicklung des Raeuber-Beute-Modells mittels Runge-Kutta-4
def RB(b0, r0, t0, tE, h, e, g):
    # Anfangsfunktionen
    def fB(nB, nR):
        y = nB*(e[0] - g[0]*nR)
        return y

    def fR(nR, nB):
        y = -nR*(e[1] - g[1]*nB)
        return y
    
    # Berechnen der Arrays
    nb = [b0]
    nr = [r0]
    M = int(tE/h)
    i = 0
    while (i<M):
        k1r=fR(nr[i], nb[i])
        k1b=fB(nb[i], nr[i])
        
        k2r=fR(nr[i] + 0.5*h*k1r, nb[i] + 0.5*h*k1b)
        k2b=fB(nb[i] + 0.5*h*k1b, nr[i] + 0.5*h*k1r)
        
        k3r=fR(nr[i] + 0.5*h*k2r, nb[i] + 0.5*h*k2b)
        k3b=fB(nb[i] + 0.5*h*k2b, nr[i] + 0.5*h*k2r)
        
        k4r=fR(nr[i] + h*k3r, nb[i] + h*k3b)
        k4b=fB(nb[i] + h*k3b, nr[i] + h*k3r)
        
        nb.append( nb[i] + (h/6.)*(k1b + 2*k2b + 2*k3b + k4b))
        nr.append( nr[i] + (h/6.)*(k1r + 2*k2r + 2*k3r + k4r))
        i = i+1
    
    T = np.arange(t0, tE+h, h)
    
    return nb, nr, T
# %% Herumspielen
    
nb, nr, T = RB(100, 2, 0, 50, 0.025, [5, 0.8], [0.02, 0.021])
  
# Modell 1
beute01 = 1490
raeuber01




  
# Anfangswerte variieren
plt.xlabel("Zeit")
plt.ylabel("Population")
plt.plot(T, nr, color='red')
plt.plot(T, nb, color='lightgreen')
plt.grid()

   
        
