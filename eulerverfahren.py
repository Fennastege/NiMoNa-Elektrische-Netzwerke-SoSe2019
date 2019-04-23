import numpy as np
import matplotlib.pyplot as plt

#dN/dt=r N (1- N/K)
# Variablen definieren
r=20.
K=100.
h=0.01
N0=5.

# Funktion f(x,t) (s. Aufzeichnungen) definieren
# #Euler-Verfahren: N_(j+1)=N_j+h*f(N_j,t_j)
def f(N):
    return r*N*(1-N/K)


#Arrays zum plotten definieren
x1,y1=[0],[N0]
i=0

#Punkte berechnen lassen
while x1[i]<=5:
    x1.append(i*h)
    y1.append(y1[i]+h*f(y1[i]))
    i+=1
i=0

#Plotten
plt.plot(x1,y1, ':', label=("r=20, h=0.01,   K=100, N_0=5"))
#plt.show()


# Variablen definieren
r=20.
K=100.
h=0.01
N0=20.

#Arrays zum plotten definieren
x2,y2=[0],[N0]

#Punkte berechnen lassen
while x2[i]<=5:
    x2.append(i*h)
    y2.append(y2[i]+h*f(y2[i]))
    i+=1
i=0

#Plotten
plt.plot(x2,y2, ':', label=("r=20, h=0.01,   K=100, N_0=20"))
#plt.show()


# Variablen definieren
r=20.
K=200.
h=0.01
N0=5.

#Arrays zum plotten definieren
x3,y3=[0],[N0]

#Punkte berechnen lassen
while x3[i]<=5:
    x3.append(i*h)
    y3.append(y3[i]+h*f(y3[i]))
    i+=1
i=0

#Plotten
plt.plot(x3,y3, ':', label=("r=20, h=0.01,   K=200, N_0=5"))
#plt.show()


# Variablen definieren
r=20.
K=100.
h=0.1
N0=5.

#Arrays zum plotten definieren
x4,y4=[0],[N0]

#Punkte berechnen lassen
while x4[i]<=5:
    x4.append(i*h)
    y4.append(y4[i]+h*f(y4[i]))
    i+=1
i=0

#Plotten
plt.plot(x4,y4, ':', label=("r=20, h=0.1,     K=100, N_0=5"))
plt.legend()
plt.xlabel('Zeit')
plt.ylabel('Populationsgröße')
plt.show()