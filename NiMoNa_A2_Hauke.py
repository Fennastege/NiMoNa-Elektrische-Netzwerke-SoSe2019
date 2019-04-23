# -*- coding: utf-8 -*-

# %% Preamble ### NiMoNa_A2_1, Hauke Bents

import numpy as np
import matplotlib.pyplot as plt

# %% Vorarbeit

r = 20.
h = np.array([0.01, 0.01, 0.01, 0.1])
K = np.array([100., 100., 200., 100.])
N_0 = np.array([5., 20., 5., 5.])


t0 = 0
t_end = 5

# Definiere Funktion ( vgl. f(N,t) in Besprechung )

def func(n,k):
    y = r*n*(1-n/k)
    return y
# %% Modellieren und plotten (Euler-Verfahren)
name = np.array(["Euler_a.png", "Euler_b.png", "Euler_c.png", "Euler_d.png"])

# Eulerverfahren
j=3
N = np.array([ N_0[j] ])
M = t_end/h[j]

i=0
while (i<M):
    N=np.append(N, N[i] + h[j]*func( N[i], K[j] ))
    i=i+1

# Plotten
xlbl = "t"
ylbl = "N(t)"

lbl = "r=" + str(r) + ", h=" + str(h[j]) + ", K=" + str(K[j]) + ", $N_0$=" + str(N_0[j])

plt.xlabel(xlbl)
plt.ylabel(ylbl)

T = np.arange(t0, t_end + h[j], h[j])

plt.plot(T, N, linestyle='', marker = '.', label=lbl)
plt.legend()
plt.grid()
plt.savefig(name[j], dpi=400)
# %% Modellieren und Plotten (Heun-Verfahren)
name2 = np.array(["Heun_a.png", "Heun_b.png", "Heun_c.png", "Heun_d.png"])

# Definiere Heun-Trapez
def Trapez(n,k,H):
    y = 0.5*H*(func(n + H*func(n,k), k) + func(n,k))
    return y

# Heun-Verfahren
j=3
N = np.array([ N_0[j] ])
M = t_end/h[j]

i=0
while (i<M):
    N=np.append(N, N[i] + Trapez(N[i], K[j], h[j]) )
    i=i+1
       
# Plotten
xlbl = "t"
ylbl = "N(t)"

lbl = "r=" + str(r) + ", h=" + str(h[j]) + ", K=" + str(K[j]) + ", $N_0$=" + str(N_0[j])

plt.xlabel(xlbl)
plt.ylabel(ylbl)

T = np.arange(t0, t_end + h[j], h[j])

plt.plot(T, N, linestyle='', marker = '.', label=lbl)
plt.legend()
plt.grid()
plt.savefig(name2[j], dpi=400)