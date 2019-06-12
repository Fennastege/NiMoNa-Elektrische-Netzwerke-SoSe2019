import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

def kuramoto(phasen, frequenz, args):
    """Kuramoto Gleichungen mit Reibungstherm. phasen sind die Phasen, frequenz die Frequenzen und args=[K, P],
    wobei K die Adjazenzmatrix ist und P die Leistungen.
    Der Rückgabewert sind die Differenzialgleichung für Phasen und Frequenz"""
    K, P = args
    summe = np.zeros(len(phasen))#Initialisieren des Summenvektors
    for m,j,l in zip(*sp.find(K)):#Für alle Einträge von K die nicht 0 sind
        summe[m] +=5.0*l*np.sin(phasen[j]-phasen[m]) #Auf die Summe addieren
    return frequenz, P - 0.1*frequenz + summe #Die Differentialgleichungen zurückgeben


def rk4(f, p1, p2, h, args=[]):
    """Runge-Kutta 4 Verfahren für die Zwei Differenzialgleichungen in f mit Vektoriellen x-Werten p1 und p2, 
    Schrittweite h und eventuell weitere für die Funktion wichtige Parameter als array in args übergeben"""
    k1, l1 = f(p1, p2, args) #Die k1 für beide Funktionen ermitteln (l1 ist das k1 der zweiten Funktion)
    k2, l2 = f(p1 + h/2*k1, p2 + h/2*l1, args)  #Die k2 für beide Funktionen ermitteln
    k3, l3 = f(p1 + h/2*k2, p2 + h/2*l2, args) #Die k3 für beide Funktionen ermitteln
    k4, l4 = f(p1 + h*k3, p2 + h*l3, args) #Die k4 für beide Funktionen ermitteln
    ks = k1 + 2*k2 + 2*k3 + k4 #Summe der ks für die erste Funktion
    ls = l1 + 2*l2 + 2*l3 + l4 ##Summe der ks für die erste Funktion
    return p1 + h/6*(ks), p2 + h/6*(ls) #Neue Werte nach dem Runge-Kutta 4 Vervahren zurückgeben


def o(p):
    """Bestimmt den Ordnungsparameter, wenn p der Vektor der Phasen ist"""
    N = len(p)
    re = 0
    im = 0
    for i in p:
        re += np.cos(i)/N
        im += np.sin(i)/N
    r = np.sqrt(re**2 + im**2)
    phi = np.arccos(re/r)
    return r, phi


def plot(phi, h, T, skip):
    """phi ist die Matrix der Phasen zu allen Zeitpunkten, h die Schrittweite, T die Simulationsdauer 
    und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    plt.figure("phases")
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(np.sin(an), np.cos(an)) #Kreis zeichnen
    plt.draw()
    plt.show(block=False)
    for i in range(1, int(T/h)):#Über die Zeitschritte iterieren
        if(i%(skip+1) == 0):#Nur plotten, wenn i ein vielfaches von skip ist
            plt.title("t = " + str(round(h*i, 2)) + "s")#Überschrift mit der Aktuellen Zeit erstellen
            t = plt.plot(np.sin(phi[i]), np.cos(phi[i]), "ro")#Punkte einzeichnen
            r, winkel = o(phi[i])#Ordnungsparameter bestimmen
            t2 = plt.arrow(0, 0, r*np.sin(winkel), r*np.cos(winkel), head_width=0.05)#Pfeil zeichnen
            plt.draw()#Zeichnen
            plt.pause(0.01)#Warten
            t.pop(0).remove()#Entfernen der Alten Punkte
            t2.remove()#Entfernen des alten Pheils


def netz(T, h, message, K, P, skip):
    """T ist die Simulationsdauer, h die Schrittweite, message ist abstand in dem der Fortschritt ausgegeben wird, K die Adjazenzmatrix, 
    P der Lesitungsvektor und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    phi = np.zeros(((int)(T/h), len(P))) #Phi hat die Form einer Matrix bei der der erste Indize den Zeitpunkt und der zweite den Oszilator angibt
    phipunkt = np.zeros(((int)(T/h), len(P))) #Phipunkt hat die gleiche Form
    phi[0] = 0*np.random.random(len(P))*2*np.pi #Anfangswinkel werden normalverteilt
    phipunkt[0] = 0*np.random.standard_cauchy(len(P)) #Anfangsfrequenzen sind normalverteilt
    pro = 0 #Prozentzahl des Aktuellen Berechnungsfortschritts
    mess = 0 #Prozentzahl, bei der die letzte Benachrichtigung geprinted wurde
    for i in range(1, int(T/h)):
        phi[i], phipunkt[i] = rk4(kuramoto, phi[i-1], phipunkt[i-1], h, [K, P]) #Berechnen der nächsten Werte
        pro = 100*(i*h/T) #Aktuellen Fortschritt errechnen
        if (mess + message < pro): #Überprüfen, ob eine neue Nachricht geprinted werden soll
            print(str(round(pro, 0))+"% Fortschritt")
            mess = pro
    plot(phi, h, T, skip)#Plotten der Ergebnisse


K = np.load('romAdj.npy')#Laden der Adjazenzmatrix
P=[]#Leistungen abwechselnd 1 und -1 setzen
for i in range((int)(len(K[0])/2)):
    P.append(1)
    P.append(-1)

netz(20, 0.01, 1, K, P, 10)#Hauptfunktion mit entsprechenden Werten aufrufen
