# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:35:38 2019

@author: Tim
"""

import numpy as np
import scipy.sparse as sp
import Kaskadenausfall as kask

adj = np.load("adjamat_nrw_1.npy")
leistung = np.load("synchronisierter_zustand_leistung_mit_leistung.npy")
theta = np.load("synchronisierter_zustand_theta.npy")

anzahl = 95
# modell_elektrisches_netzwerk(theta_anfang, adja, Leistung, i_Ausfallleitung, j_Ausfallleitung, Anzahl_ausfallgrenze_Leitungen = 5, technische_maximallast = 0.1, Kopplung = 30, 
# genauigkeit_output = 500, synchrogrenze = 0.999, synchronisiert_toleranz = 0.0001):
synchron = np.array([None]*anzahl)
ausgefallen = np.array([1]*anzahl)
zeit = np.array([0.0]*anzahl)
ordnungsding = np.array([0.0]*anzahl)


Kopplungsgrad = 30
maxlast = 0.001
i = 1

for k,l,m in zip(*sp.find(adj)):
    if k <= l:
        synchron[i],ausgefallen[i],zeit[i],ordnungsding[i] = kask.modell_elektrisches_netzwerk(theta, adj, leistung, i_Ausfallleitung = k, j_Ausfallleitung = l, technische_maximallast = maxlast, Kopplung = Kopplungsgrad)
        print(" Nummer:",i," : ",k,"-",l,"||", synchron[i],ausgefallen[i],zeit[i],ordnungsding[i])
        i += 1
