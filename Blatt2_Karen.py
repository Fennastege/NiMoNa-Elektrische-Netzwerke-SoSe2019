#!/usr/bin/env python
# coding: utf-8

# In[27]:


import numpy as np
import matplotlib.pyplot as plt


# In[28]:


def f(N, r, K):
    return N*r*(1-(N/K))


# In[63]:


def ploteuler(r, h, K, N0):
    heun=False
    N = []
    N.append(N0)
    for i in range(int(5/h)):
        N1 = N[i] + h*f(N[i], r, K)
        N.append(N1)
        
    label = "h="+str(h)+", K="+str(K)+", N0="+str(N0)
    plt.plot(np.arange(0, 5+h, h), N, label=(label))
   


# In[64]:


def plot1():
    ploteuler(20, 0.01, 100, 5)
    ploteuler(20, 0.01, 100, 20)
    ploteuler(20, 0.01, 200, 5)
    ploteuler(20, 0.1, 100, 5)


# In[65]:


def plotheun(r, h, K, N0):
    N = []
    N.append(N0)
    for i in range(int(5/h)):
        N1 = N[i] + h*f(N[i], r, K)
        N.append( N[i] + 0.5*h*(f(N[i], r, K) + f(N1, r, K)) )

    label2 = "h="+str(h)+", K="+str(K)+", N0="+str(N0)+"; Heun"    
    plt.plot(np.arange(0, 5+h, h), N, label=(label2))        
   
def plot2():
    plotheun(20, 0.1, 100, 5)        


# In[66]:


plot1()
plot2()
plt.legend()


# In[ ]:




