#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:50:40 2019

@author: henri
"""

#Importiere mehrere Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
#für Ordnerstruktur
import time 
import os


postmat = np.load('posmat_nrw_64.npy')
#adjamat1 = np.ndarray(shape=(N,N),0)

def init(connections):
    #Funktion die auf Basis von Verbidnungen Adjazenzmatrix initlaisiert
    #connections: Array/Tupel mit Verbidnungen zwischen Knoten
    #N: Anzahl der Knoten
    N = len(connections)
    adjamat1 = np.full((N,N), 0)
    
    i = 0
    for i in range(0,N):
        j = 0
        for j in range(0,len(connections[i][:])):
            adjamat1[i,(connections[i][j])-1] = 1
            adjamat1[(connections[i][j])-1,i] = 1
            j+=1
        i+=1
    return adjamat1

#Zeichen Funktion zur Kontrolle       
def zeichnen(adjamat1,postmat):     
    m=0
    j=0
    l=0
    #z gleich null für Nummer der Adjazenzmatrix
    #z=1 für Schaubild
    z=1
    plt.figure(figsize=(10,10))
    for m,j,l in zip(*sp.find(adjamat1)):
        plt.plot(postmat[[m,j],0],postmat[[m,j],1], "-", color="black",linewidth=1.05*l)
    for l in range(0,N):  
        
        plt.plot(postmat[l,0],postmat[l,1],'o',color='blue',marker='.')
        plt.text(postmat[l,0],postmat[l,1],str(z),size='large')
        z+=1
        
#Eingabe der Connections
Z1 = [3]
Z2 = [3,7,6]
Z3 = [1,4,5,6]
Z4 = [3]
Z5 = [11,3,10]
Z6 = [3,10,8,7,2]
Z7 = [2,8,6]
Z8 = [7,12,6,9]
Z9 = [10,12,8]
Z10 = [6,21,9,13,5]
Z11 = [5]
Z12 = [8,9,13]
Z13 = [12,10,20]
Z14 = [16,17,15]
Z15 = [16,14,17]
Z16 = [14,15,17]
Z17 = [16,18,14,15]
Z18 = [17,19,28]
Z19 = [18,23,20]
Z20 = [19,13,23,24]
Z21 = [22,25,10]
Z22 = [21]
Z23 = [24,20,19,28]
Z24 = [23,20,29,28]
Z25 = [21,26,27]
Z26 = [32,25]
Z27 = [25,31,33]
Z28 = [34,18,23,29,24]
Z29 = [31,30,24,28]
Z30 = [31,29,34]
Z31 = [30,29,27]
Z32 = [26]
Z33 = [27,51]
Z34 = [30,37,28]
Z35 = [50,36,53]
Z36 = [38,49,35,50]
Z37 = [38,40,39,34]
Z38 = [36,37]
Z39 = [40,37,41]
Z40 = [37,39,41,45]
Z41 = [39,44,43,40,42]
Z42 = [43,41]
Z43 = [42,41,44]
Z44 = [41,47,43]
Z45 = [40,46]
Z46 = [57,58,45]
Z47 = [44]
Z48 = [59]
Z49 = [53,54,36]
Z50 = [51,53,35,36]
Z51 = [52,33]
Z52 = [51]
Z53 = [50,55,49,35]
Z54 = [49,57,56]
Z55 = [56,53]
Z56 = [55,54]
Z57 = [46,54]
Z58 = [46,60,59]
Z59 = [48,61,60,62,58]
Z60 = [59,62,58,61]
Z61 = [59,62,60]
Z62 = [60,61,63,59]
Z63 = [62,64]
Z64 = [63]

Eingabe = [Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,Z21,Z22,Z23,Z24,Z25,Z26,Z27,Z28,Z29,Z30,Z31,Z32,Z33,Z34,Z35,Z36,Z37,Z38,Z39,Z40,Z41,Z42,Z43,Z44,Z45,Z46,Z47,Z48,Z49,Z50,Z51,Z52,Z53,Z54,Z55,Z56,Z57,Z58,Z59,Z60,Z61,Z62,Z63,Z64]

#Berechnung und Zeichnen der Adjamat
zeichnen(init(Eingabe),postmat)

#Nur Adjamat berechnen fürs abspeichern
adjamat1 = init(Eingabe)

#Aktivieren zum Abspeichern
#Nicht vergessen: Namen danach verändern
#np.save("adjamat_nrw_1.npy",adjamat1)