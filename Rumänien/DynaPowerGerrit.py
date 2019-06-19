import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

def kuramoto(x, K, P):
    """Kuramoto Gleichungen mit Reibungstherm. phasen sind die Phasen, frequenz die Frequenzen und args=[K, P],
    wobei K die Adjazenzmatrix ist und P die Leistungen.
    Der Rückgabewert sind die Differenzialgleichung für Phasen und Frequenz"""
    phasen, frequenz = x
    summe = np.zeros(len(phasen))#Initialisieren des Summenvektors
    for m,j,l in zip(*sp.find(K)):#Fuer alle Eintraege von K die nicht 0 sind
        summe[m] +=5.0*l*np.sin(phasen[j]-phasen[m]) #Auf die Summe addieren
    return np.array([frequenz, P - 0.1*frequenz + summe]) #Die Differentialgleichungen zurückgeben  


def rk4(f, x, h, K, P, i):
    """Runge-Kutta 4 Verfahren für Differenzialgleichungen in f mit Vektoriellen x-Werten,
    Schrittweite h und eventuell weitere für die Funktion wichtige Parameter als array in args übergeben"""
    k1 = f(x, K, P[i])
    k2 = f(x + h/2*k1, K, (P[i]+P[i+1])/2)  #da die Leistung nun zeitabhängig ist, muss diese auch zu verschiedenen Zeiten
    k3 = f(x + h/2*k2, K, (P[i]+P[i+1])/2)  #betrachtet werden (hier muss man u.a. einen halben Zeitschritt verschieben,
    k4 = f(x + h*k3, K, P[i+1])             #daher wurde Mittelwert zwischen beiden umliegenden Funktionswerten gebildet)
    return x + h/6*(k1 + 2*k2 + 2*k3 + k4) #Neue Werte nach dem Runge-Kutta 4 Vervahren zurückgeben


def o(p):
    """Bestimmt den Ordnungsparameter, wenn p der Vektor der Phasen ist"""
    N = len(p)
    re = 0
    im = 0
    for i in p:
        re += np.cos(i)/N
        im += np.sin(i)/N
    r = np.sqrt(re**2 + im**2)
    phi_c = np.arccos(re/r)
    phi_s = np.arcsin(im/r)
    return r, phi_c, phi_s


def plot(phi, h, T, skip, adjamatrix, posmatrix, power):
    """phi ist die Matrix der Phasen zu allen Zeitpunkten, h die Schrittweite, T die Simulationsdauer 
    und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    _, axs = plt.subplots(1, 2, figsize=(13, 6)) #erzeugt zwei Bilder nebeneinander
    an = np.linspace(0, 2*np.pi, 100)
    axs[1].plot(np.cos(an), np.sin(an)) #Kreis zeichnen
    for m,j,l in zip(*sp.find(adjamatrix)): #Verbindungslinien zeichnen
        axs[0].plot([posmatrix[m][0], posmatrix[j][0]], [posmatrix[m][1], posmatrix[j][1]], "-", color="black",linewidth=1.05*l)
    plt.draw()
    plt.show(block=False)
    for t in range(0, int(T/h)):#ueber die Zeitschritte iterieren
        if(t%(skip)== 0):#Nur plotten, wenn t ein vielfaches von skip ist
            axs[1].title.set_text("t = " + str(round(h*t, 2)) + "s")#ueberschrift mit der Aktuellen Zeit erstellen
            axs[0].title.set_text(r"$\sum_{i}P_i=$ "+ str(round(sum(power[t]))))
            r, winkel_c, winkel_s = o(phi[t])#Ordnungsparameter bestimmen
            pfeil = axs[1].arrow(0, 0, r*np.cos(winkel_c), r*np.sin(winkel_s), head_width=0.05)#Pfeil zeichnen
            points = []
            for i in range(len(phi[t])):
                if( (abs(np.cos(phi[t][i])-np.cos(winkel_c)) <0.1) and (abs(np.sin(phi[t][i])-np.sin(winkel_s)) <0.1) and (r > 0.8)):
                    p = axs[1].plot(np.cos(phi[t][i]), np.sin(phi[t][i]), "ro")#Punkt einzeichnen
                    points.append(p)
                    p = axs[0].plot(posmatrix[i][0], posmatrix[i][1], "ro")#Punkt einzeichnen
                    points.append(p)
                else:
                    p = axs[1].plot(np.cos(phi[t][i]), np.sin(phi[t][i]), "bo")#Punkt einzeichnen
                    points.append(p)
                    p = axs[0].plot(posmatrix[i][0], posmatrix[i][1], "bo")#Punkt einzeichnen
                    points.append(p)
            plt.draw()#Zeichnen
            plt.pause(0.01)#Warten
            for point in points:
                point.pop(0).remove()#Entfernen der Alten Punkte
            pfeil.remove()#Entfernen des alten Pfeils


def netz(T, h, K, skip, pos):
    """T ist die Simulationsdauer, h die Schrittweite, K die Adjazenzmatrix, 
    P der Lesitungsvektor und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    P=power(T, h, K)
    phi = np.zeros(((int)(T/h), K.shape[0])) #Phi hat die Form einer Matrix, bei der der erste Index den Zeitpunkt und der zweite den Oszilator angibt
    phipunkt = np.zeros(((int)(T/h), K.shape[0])) #Phipunkt hat die gleiche Form
    phi[0] = np.random.random(K.shape[0])*2*np.pi #Anfangswinkel werden normalverteilt
    phipunkt[0] = 0.1*np.random.standard_cauchy(K.shape[0]) #Anfangsfrequenzen sind normalverteilt
    pro = 0 #Prozentzahl des Aktuellen Berechnungsfortschritts
    mess = 0 #Prozentzahl, bei der die letzte Benachrichtigung geprinted wurde
    for i in range(1, int(T/h)):
        phi[i], phipunkt[i] = rk4(kuramoto, [phi[i-1], phipunkt[i-1]], h, K, P, i) #Berechnen der nächsten Werte
        pro = 100*(i*h/T) #Aktuellen Fortschritt errechnen
        if (mess + 1 < pro): #überprüfen, ob eine neue Nachricht geprinted werden soll
            print(str(round(pro, 0))+"% Fortschritt")
            r, _, _ = o(phi[i])#Ordnungsparameter bestimmen
            print (r)
            mess = pro
    plot(phi, h, T, skip, K, pos, P)#Plotten der Ergebnisse

def power(T, h, K):
    P = np.zeros(((int)(T/h)+1, K.shape[0]))
    for t in range (int(T/h)+1):
        P[t]=Leistungsverlauf(t, K)
    return np.array(P)

def Leistungsverlauf(Zeit, K):
    leistung=np.zeros(K.shape[0])   #hier wird die Leistung zunächst als konstant definiert,
    leistung[0]=-0.1                #sie kann jedoch auch als Funktion der Zeit definiert werden
    leistung[1]=-0.05
    leistung[2]=-0.1
    leistung[3]=-0.05
    leistung[4]=-0.35
    leistung[5]=-0.1
    leistung[6]=-0.2
    leistung[7]=-0.05
    leistung[8]=0.3
    leistung[9]=0.5
    leistung[10]=0.2
    return leistung

K = [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
     [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ]
K=sp.csr_matrix(K)
pos = [[1.5, 3], [3.5, 3], [4, 2.5], [4, 1.5], [3.5, 1], [1.5, 1], [1, 1.5], [1, 2.5], [2, 2], [3, 2], [2, 4]]

netz(100, 0.01, K, 50, pos)#Hauptfunktion mit entsprechenden Werten aufrufen


