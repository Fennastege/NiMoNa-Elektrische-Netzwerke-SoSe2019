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
Z2 = []
Z3 = []
Z4 = []
Z5 = []
Z6 = []
Z7 = []
Z8 = []
Z9 = []
Z10 = []
Z11 = []
Z12 = []
Z13 = []
Z14 = []
Z15 = []
Z16 = []
Z17 = []
Z18 = []
Z19 = []
Z20 = []
Z21 = []
Z22 = []
Z23 = []
Z24 = []
Z25 = []
Z26 = []
Z27 = []
Z28 = []
Z29 = []
Z30 = []
Z31 = []
Z32 = []
Z33 = []
Z34 = []
Z35 = []
Z36 = []
Z37 = []
Z38 = []
Z39 = []
Z40 = []
Z41 = []
Z42 = []
Z43 = []
Z44 = []
Z45 = []
Z46 = []
Z47 = []
Z48 = []
Z49 = []
Z50 = []
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
