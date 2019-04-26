import numpy as np
import matplotlib.pyplot as plt
%matplotlib notebook
import math
import random

# Aufgabe 1
def Euler(r=0, h=0, K=0, N_0=0):
    N=[N_0]
    t=[0]
    
    i=0
    while i < (5/h):
        t.append(h+i*h)
        N.append(N[i]+h*r*N[i]*(1-N[i]/K))
        i=i+1
        
    plt.figure()
    plt.ylabel('N(t)')
    plt.xlabel('t')
    plt.title("Die explizite Euler-Lösung")
    plt.plot(t, N,'g')
    plt.grid(True)

Euler(20, 0.01, 100, 5)
Euler(20, 0.01, 100, 20)
Euler(20, 0.01, 200, 5)
Euler(20, 0.1, 100, 5)

# Aufgabe 2
def Heun(r=0, h=0, K=0, N_0=0):
    L=0
    N=[N_0]
    t=[0]
    
    i=0
    while i < (5/h):
        t.append(h+i*h)
        L=N[i]+h*r*N[i]*(1-N[i]/K)
        N.append(N[i]+0.5*h*r*(N[i]*(1-N[i]/K)+L*(1-L/K)))
        i=i+1
        
    plt.figure()
    plt.ylabel('N(t)')
    plt.xlabel('t')
    plt.title("Die Lösung des Heun-Verfahrens")
    plt.plot(t, N,'g')
    plt.grid(True)

Heun(20, 0.01, 100, 5)
Heun(20, 0.01, 100, 20)
Heun(20, 0.01, 200, 5)
Heun(20, 0.1, 100, 5)

# Es fällt auf, dass die Ergebnisse der beiden Verfahren bei kleinen Schrittweiten voneinander abweichen. Das Ergebnis des 
# Heun-Verfahrens ist nicht so gezackt wie das Ergebnis des expliziten Euler-Verfahrens. Der Graph nach ersterem Verfahren 
# liegt außerdem immer unter dem Graphen des anderen Verfahrens. 