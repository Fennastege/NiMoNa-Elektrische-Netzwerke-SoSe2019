import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd

epsilon1, epsilon2, gamma1, gamma2, h = 2.0, 0.8, 0.002, 0.0002, 0.025
#Startpopulation beider Parteien zum Zeitpunkt time=0:
p1_0, p2_0= 1000,20
time=[0]

#Array zur Beschreibung der Entwicklung der Population
p1=[p1_0]
p2=[p2_0]

def nxstep(f,array):
    """Errechnet ausgehend von einer Funktion f und einem Array, welches die Gesamtpopulation beschreibt,
    den nächsten Zustand der Population"""
    return array[-1]+h/6*(kn(1,f,array)+2*kn(2,f,array)+2*kn(3,f,array)+kn(4,f,array))

def kn(n,f,array):
    """beschreibt die Funktion kn, wobei 0<n<5 gilt."""
    if n==1:
        result= f(array[-1])
    elif n==2:
        result=f(array[-1]+h/2*kn(1,f,array))
    elif n==3:
        result= f(array[-1]+h/2*kn(2,f,array))
    else:
        result= f(array[-1]+h*kn(3,f,array))
    return result

def f1(aktuelle_pop):
    """Beschreibung der zeitl. Änderung der Population p1"""
    return aktuelle_pop*(epsilon1-gamma1*p2[-1])

def f2(aktuelle_pop):
    """Beschreibung der zeitl. Änderung der Population p2"""
    return -p2[-1]*(epsilon2-gamma2*p1[-1])

while time[-1]<50:
    p1_nxt= nxstep(f1,p1)
    p2_nxt= nxstep(f2,p2)
    p1.append(p1_nxt)
    p2.append(p2_nxt)
    time.append(time[-1]+h)

plt.plot(time,p1,'-', label='p1')
plt.plot(time,p2, '-', label='p2')
plt.legend()
plt.title("Der zeitliche Verlauf der zwei Populationen")
plt.grid()
plt.show()

plt.plot(p1,p2)
plt.title("Beide Populationen gegeneinander aufgetragen")
plt.grid()
plt.show()