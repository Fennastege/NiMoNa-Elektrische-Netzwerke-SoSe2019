# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:35:38 2019

@author: Tim
"""

import numpy as np
import scipy.sparse as sp
import Kaskadenausfall as kask

adj = np.load("adjamat_nrw_1.npy")
leistung = np.load("synchronisierter_zustand_leistung.npy")
theta = np.load("synchronisierter_zustand_theta.npy")

anzahl = 95
# modell_elektrisches_netzwerk(theta_anfang, adja, Leistung, i_Ausfallleitung, j_Ausfallleitung, Anzahl_ausfallgrenze_Leitungen = 5, technische_maximallast = 0.1, Kopplung = 30,
# genauigkeit_output = 500, synchrogrenze = 0.999, synchronisiert_toleranz = 0.0001):
synchron = np.array([None]*anzahl)
ausgefallen = np.array([1]*anzahl)
zeit = np.array([0.0]*anzahl)
ordnungsding = np.array([0.0]*anzahl)
welche_rausgenommen = np.ndarray(shape=(anzahl,2))

Kopplungsgrad = 50
maxlast = 0.05
i = 0

#print(adj)
k = 0
l = 0

#Deepcopy

#synchron[i],ausgefallen[i],zeit[i],ordnungsding[i] = kask.modell_elektrisches_netzwerk(theta, adj, leistung, i_Ausfallleitung = k, j_Ausfallleitung = l, technische_maximallast = maxlast, Kopplung = Kopplungsgrad)
#welche_rausgenommen[i,:] = [k,l]
#print(" Nummer:",i," : ",k,"-",l,"||", synchron[i],ausgefallen[i],zeit[i],ordnungsding[i])
#i += 1

for k,l,m in zip(*sp.find(adj)):
    if i >= (anzahl-1):
        exit
    else:
        if k >= l:
            #print("vor",adj)
            synchron[i],ausgefallen[i],zeit[i],ordnungsding[i] = kask.modell_elektrisches_netzwerk(theta, adja = np.load("adjamat_nrw_1.npy"), Leistung = leistung, i_Ausfallleitung = k, j_Ausfallleitung = l, technische_maximallast = maxlast, Kopplung = Kopplungsgrad)
            #print("danach",adj)
            welche_rausgenommen[i,:] = [k,l]
            print(" Nummer:",i," : ",k,"-",l,"||", synchron[i],ausgefallen[i],zeit[i],ordnungsding[i])
            i += 1

#np.savetxt("ergebnisse.txt",[welche_rausgenommen[:,0],welche_rausgenommen[:,1],synchron,ausgefallen,zeit,ordnungsding],header="Gesamtergebnisse \n Verbindung | Synchronisiert | Wie viele Oszis ausgefallen | Dauer der Simulation (in s) | Ordnungsparameter am Ende")
