# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 18:14:34 2019

@author: Tim
"""

#Importiere mehrere Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
#für Ordnerstruktur
import time 
import os

#Adjazenz-Matrix laden
adjazenzmatrix = np.load("adjamat_nrw_1.npy")

#Anzahl der Oszillatoren (berechnen aus Größe der Adjazenzmatrix)
N = len(adjazenzmatrix[0])


def clusterLokal(adj,i):
    nachbarn = np.array([])
    #spalte = 0
    anzVerb = 0
    #Nachbarn bestimmen
    for spalte in range(0,len(adj)):
        if (adj[i,spalte] != 0):
            nachbarn = np.append(nachbarn,spalte)
    #Verbindungen der Nachbarn bestimmen
    for m,j,l in zip(*sp.find(adj)):
        if m in nachbarn and j in nachbarn:
            anzVerb = anzVerb + 0.5
    
    if len(nachbarn) < 2:
        return 0
    else:
        return 2 * anzVerb / (len(nachbarn) * (len(nachbarn) - 1))

def clusterGlobal(clusterLokal,adj):
    summe = 0
    for k in range(0,len(adj)):
        summe += clusterLokal(adj,k)
    return summe/len(adj)

print(clusterGlobal(clusterLokal,adjazenzmatrix))

addi = 0
for i,k,l in zip(*sp.find(adjazenzmatrix)):
    addi += 0.5

print(addi)