#!/usr/bin/env python
# coding: utf-8

# In[23]:


import matplotlib.pyplot as plt
# 1. Explizites Euler-Verfahren:
r = 20
h = 0.1
K = 100
N = 5
t = 0
def N_ext(N_old):
    return r*N_old*(1-(N_old/K))
while t<=(5*h):
    plt.plot(t,N,'bo')
    N=N+h*N_ext(N)
    t=t+h    
# Und analog für die Werte aus a), b), c) - die Kurven sind allerdings ehrlich gesagt ziemlich unspannend.
# Wie man sieht, denkt das explizite Eulerverfahren hier nicht vorausschauend genug und kommt bei t=0.4 auf einen Wert größer als
# die Kapazitätsgrenze K, was - das sagt mir der gesunde Menschenverstand - nicht korrekt sein kann. Es wäre also sinnvoller, zu
# schauen, welche Steigung sich beim errechneten Wert ergeben würde und den Mittelwert beider Steigungen zu bilden, dann erhält 
# man realistischer die tatsächliche Steigung zwischen beiden Punkten. So geschehend im Heun-Verfahren:
t=0
N=5
plt.figure()
while t<=(5*h):
    plt.plot(t,N,'bo')
    N=N+(1/2)*h*(N_ext(N)+N_ext(N+h*N_ext(N)))
    t=t+h

