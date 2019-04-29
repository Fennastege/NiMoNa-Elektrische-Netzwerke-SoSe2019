#Importiere mehrere (Standard-)Packages
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------
#Funktion für Räuber-Beute definieren
def rbmodell(beute_startwert, raeuber_startwert, ep1, gam1, ep2, gam2, tmax, h):
    
    #Zwischen-Funktionen aus Modell definieren
    def fkt1(p1_in, p2_in):
        wert = p1_in * (ep1 - gam1 * p2_in)
        return wert

    def fkt2(p1_in, p2_in):
        wert = (-1)*p2_in * (ep2 - gam2 * p1_in)
        return wert
    
    p1 = np.array([0.0]*int(tmax/h))
    p2 = np.array([0.0]*int(tmax/h))

    i = 0

    #Anfangsbedingungen
    p1[0] = beute_startwert
    p2[0] = raeuber_startwert

    #Berechnen der Werte
    while (i < int(tmax/h-1)):
        #Faktoren berechnen
        k1_1 = fkt1(p1[i], p2[i])
        k1_2 = fkt2(p1[i], p2[i])

        k2_1 = fkt1(p1[i] + h/2 * k1_1, p2[i] + h/2 * k1_2)
        k2_2 = fkt2(p1[i] + h/2 * k1_1, p2[i] + h/2 * k1_2)

        k3_1 = fkt1(p1[i] + h/2 * k2_1, p2[i] + h/2 * k2_2)
        k3_2 = fkt2(p1[i] + h/2 * k2_1, p2[i] + h/2 * k2_2)

        k4_1 = fkt1(p1[i] + h * k3_1, p2[i] + h * k3_2)
        k4_2 = fkt2(p1[i] + h * k3_1, p2[i] + h * k3_2)

        p1[i+1] = p1[i] + h/6 * (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1)
        p2[i+1] = p2[i] + h/6 * (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2)

        i = i + 1

    return p1, p2
    
#---------------------------------------------------------------
    
#Parameter für Berechnung festlegen

#Zeit-Dauer der Simulation
tmax = 20
#Zeitschritte
h = 0.025

#------------------------------------------------
#MODELL 01

#Beute
beute_Startwert = 1500
beute_Reproduktion =  2.0
beute_Sterberate = 0.02
#----------------------
#Räuber
raeuber_Startwert = 150
raeuber_Reproduktion = 0.8
raeuber_Sterberate = 0.0002

#----------------------------------------------
#Simulation laufen lassen
p1, p2 = rbmodell(beute_Startwert, raeuber_Startwert, beute_Reproduktion, beute_Sterberate, raeuber_Reproduktion, raeuber_Sterberate, tmax, h)

#Plotten
fig, axs = plt.subplots(1, 2, figsize=(15, 5))
axs[0].plot(np.arange(0,tmax,h), p1, label="Modell 01, Beute")
axs[0].plot(np.arange(0,tmax,h), p2, label="Modell 01, Räuber")
axs[1].plot(p1,p2)
plt.legend()
