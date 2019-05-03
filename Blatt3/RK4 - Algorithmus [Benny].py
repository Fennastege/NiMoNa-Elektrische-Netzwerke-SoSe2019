#!/usr/bin/env python
# coding: utf-8

# In[41]:


import matplotlib.pyplot as plt
# Runge-Kutta-Verfahren 4. Ordnung:

# Benötigte Variablen und Funktionen:
def f(x,t):
    # Hier die Funktion der DGL (x'(t) = f(x,t)) angeben:
    return -2*t*x
h = 0.025 # Schrittweite
AnzahlSchritte = 100
c=[1/6,1/3,1/3,1/6]
a=[0,1/2,1/2,1]
b=[[],[1/2],[0,1/2],[0,0,1]]
k=[0]*4
# Startwerte:
x = [1]
t = [0]

# Berechne oben festgelegte Anzahl an x- (und t-)Werten:
for step in range(0, AnzahlSchritte):
    # Berechne k-Array:
    for j in range(0,4):
        sum = 0
        for i in range(0,j-1):
            sum += b[j][i]*k[i]
        k[j] = f(x[step] + h*sum, t[step] + h*a[j])

    # Berechne neuen x-Wert:
    sum = 0
    for z in range(0,4):
        sum = sum + c[z]*k[z]
    x.append(x[step] + h * sum)
    
    # t-Update:
    t.append(t[step]+h)

# Plotten der Funktion:
plt.plot(t,x)

    

