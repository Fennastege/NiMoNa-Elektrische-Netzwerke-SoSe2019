#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 23:30:44 2019

@author: henri
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
#für Ordnerstruktur
import time 
import os

#Matrix für N Oszillatoren
N=15
adjamat1 = np.ndarray(shape=(N,N))
adjamat1 = [[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[1,1,0,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,1,0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,1,1,0,1,1,0,1,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0,0,1,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,1],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]]

#Test ob Adjazenzmatrixlisten die richtige Länge haben
for i in range(N):
    x=0
    Laenge = len(adjamat1[i])
    if Laenge == N:
        print [i],'passt!'
        x=x+0
    else:
        print [i], 'falsch <---!!', 'Länge:', len(adjamat1[i])
        x=x+1



posmat1 = np.ndarray(shape=(N,2))
posmat1 = [[20,20],[70,20],[60,40],[55,65],[55,80],[50,85],[80,70],[70,85],[80,90],[95,90],[60,120],[100,130],[120,130],[145,135],[160,110]]

if x==0:
    np.save("adjamat2_h.npy",adjamat1)
    np.save("posmat2_h.npy",posmat1)
else:
    print 'Fehler in Matrix'