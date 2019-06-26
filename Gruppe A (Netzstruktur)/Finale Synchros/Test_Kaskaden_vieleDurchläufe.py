# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 01:12:23 2019

@author: Hauke
"""

import numpy as np
import scipy.sparse as sp
import Kaskadenausfall as kask

Ordnername = ["dezentral", "zentral", "flaechig", "random", "spinne"]

q = 4

Kopplungsgrad = 50
sync = 0.99
maxlast = 0.15


path = "adja_" + Ordnername[q]
adj = np.load(path + "/adja_start.npy")
leistung = np.load(path + "/synchronisierter_zustand_leistung.npy")
theta = np.load(path + "/synchronisierter_zustand_theta.npy")

print(path + "/adja_start.npy")


# Initiieren
anzahl = len(adj)
synchron = np.array([None]*anzahl)
ausgefallen = np.array([1]*anzahl)
zeit = np.array([0.0]*anzahl)
ordnungsding = np.array([0.0]*anzahl)
welche_rausgenommen = np.ndarray(shape=(anzahl,2))


i = 0
k = 0
l = 0


for k,l,m in zip(*sp.find(adj)):
    if i >= (anzahl-1):
        exit
    else:
        if k >= l:
            #print("vor",adj)
            synchron[i],ausgefallen[i],zeit[i],ordnungsding[i] = kask.modell_elektrisches_netzwerk(theta,adja = np.load(path + "/adja_start.npy"), Leistung = leistung, i_Ausfallleitung = k, j_Ausfallleitung = l,technische_maximallast = maxlast, Kopplung = Kopplungsgrad, synchrogrenze=sync)
            #print("danach",adj)
            welche_rausgenommen[i,:] = [k,l]
            print(" Nummer:",i," : ",k,"-",l,"||", synchron[i],ausgefallen[i],zeit[i],round(ordnungsding[i]),5)
            i += 1

np.savetxt(path + "/Kaskadenausfall_ergebnisse.txt",[welche_rausgenommen[:,0],
                welche_rausgenommen[:,1],synchron,ausgefallen,zeit,ordnungsding],
                header="Gesamtergebnisse \n Verbindung | Synchronisiert | Wie viele Oszis ausgefallen | Dauer der Simulation (in s) | Ordnungsparameter am Ende")
