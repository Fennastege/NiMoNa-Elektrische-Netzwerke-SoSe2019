# BLATT 02


#Importiere mehrere (Standard-)Packages
import numpy as np
import matplotlib.pyplot as plt

#Funktion definieren
def fkt(N_input, r_inp, K_inp):
    y = r_inp * N_input * (1 - N_input / K_inp)
    return y

#Gegebene Werte der Variablen:
r_alle = np.array([20]*4)
h_alle = np.array([0.01, 0.01, 0.01, 0.1])
K_alle = np.array([100, 100, 200, 100])
N0_alle = np.array([5, 20, 5, 5])
N = np.ndarray(shape=(4, int(5/h_alle[0])))
N_heun = np.ndarray(shape=(4, int(5/h_alle[0])))

#-------------------------------------------------------------
#AUFGABE 1

#for-Schleife für die Teilaufgaben
j = 0
for j in range(0,4):
    #Variablen übernehmen
    h = h_alle[j]
    K = K_alle[j]
    r = r_alle[j]
    #Initialisieren der Funktionswerte:
    N[j, 0] = N0_alle[j]
    #Beginnen des eigentlichen Euler-Verfahrens
    t = 0
    i = 0
    while(t < (5.0 - 2*h)):
        N[j, i+1] = N[j, i] + h * fkt(N[j, i], r, K)
        i = i+1
        t = t + h
        
#Plotten
for j in range(0,4):
    plt.plot(np.arange(0, 5, h_alle[j]), N[j][0:int(5/h_alle[j])], "-")

#Bennenung und Labeln
plt.xlabel("Zeit $t$")
plt.ylabel("Anzahl $N(t)$")


#-------------------------------------------------------------
#AUFGABE 2

#for-Schleife für die Teilaufgaben
j = 0
for j in range(0,4):
    #Variablen übernehmen
    h = h_alle[j]
    K = K_alle[j]
    #Initialisieren der Funktionswerte:
    #N = np.array([0.0]*(int(5/h)))
    N_heun[j, 0] = N0_alle[j]
    #Beginnen des eigentlichen Euler-Verfahrens
    t = 0
    i = 0
    while(t < (5.0 - 2*h)):
        #predictor-corrector Verfahren ()
        N_heun[j, i+1] = N_heun[j, i] + h *(1/2) * (fkt(N_heun[j, i], r, K) + fkt(N_heun[j, i] + h * fkt(N_heun[j, i], r, K), r, K))
        i = i+1
        t = t + h
        

#plotten
for j in range(0,4):
    plt.plot(np.arange(0, 5, h_alle[j]), N_heun[j][0:int(5/h_alle[j])], "-")

#Bennenung und Labeln
plt.xlabel("Zeit $t$")
plt.ylabel("Anzahl $N(t)$")

#-------------------------------------------------------------
#VERGLEICH

# j einstellen für Teilaufgabe, xlim für Funktionsaussschnitt
#j = 3
#plt.xlim(0,5)

#plt.plot(np.arange(0, 5, h_alle[j]), N[j][0:int(5/h_alle[j])], ".")
#plt.plot(np.arange(0, 5, h_alle[j]), N_heun[j][0:int(5/h_alle[j])], "--")

#Besonders im letzten Fall mit geringerer Diskretisierungsrate h fällt der größere Fehler des Euler-Verfahrens im Vergleich zum Heun-Verfahren auf.
#In den anderen Fällen gibt es Abweichungen zu Beginn (bis ca t = 0.5), die sich dann aber später annähern.
