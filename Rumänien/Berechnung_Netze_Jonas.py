import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import sys
import matplotlib.image as mpimg
import os
import matplotlib.patches as patch

def kuramoto(x, nr, K, P, Pn):
    """Kuramoto Gleichungen mit Reibungstherm.
    Der Rückgabewert sind die Differenzialgleichung für Phasen und Frequenz"""
    phasen, frequenz = x
    if nr==4: P = Pn
    elif not nr == 1: P = np.divide(np.add(P, Pn), 2) 
    summe = np.zeros(len(phasen))#Initialisieren des Summenvektors
    for m,j,l in zip(*sp.find(K)):#Fuer alle Eintraege von K die nicht 0 sind
        summe[m] +=50.0*l*np.sin(phasen[j]-phasen[m]) #Auf die Summe addieren
    return np.array([frequenz, P - 0.1*frequenz + summe]) #Die Differentialgleichungen zurückgeben  


def rk4(f, x, h, *args):
    """Runge-Kutta 4 Verfahren für Differenzialgleichungen in f mit Vektoriellen x-Werten,
    Schrittweite h und eventuell weitere für die Funktion wichtige Parameter"""
    k1 = f(x, 1, *args)
    k2 = f(x + h/2*k1, 2, *args)
    k3 = f(x + h/2*k2, 3, *args)
    k4 = f(x + h*k3, 4, *args)
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


def plotphi(phi, h, T, skip, adjamatrix, posmatrix, P):
    """phi ist die Matrix der Phasen zu allen Zeitpunkten, h die Schrittweite, T die Simulationsdauer 
    und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    _, axs = plt.subplots(1, 2, figsize=(13, 6)) #erzeugt zwei Bilder nebeneinander
    an = np.linspace(0, 2*np.pi, 100)
    axs[1].plot(np.sin(an), np.cos(an)) #Kreis zeichnen
    for m,j,l in zip(*sp.find(adjamatrix)): #Verbindungslinien zeichnen
        axs[0].plot([posmatrix[m][0], posmatrix[j][0]], [posmatrix[m][1], posmatrix[j][1]], "-", color="black",linewidth=1*l)
    axs[0].imshow(mpimg.imread('deutschland.png'), extent=[-3, 17, 0, 20])
    plt.draw()
    plt.show(block=False)
    for t in range(0, int(T/h)):#ueber die Zeitschritte iterieren
        if(t%(skip)== 0):#Nur plotten, wenn t ein vielfaches von skip ist
            plt.title("t = " + str(round(h*t, 2)) + "s")#ueberschrift mit der Aktuellen Zeit erstellen
            r, winkel_c, winkel_s = o(phi[t])#Ordnungsparameter bestimmen
            pfeil = axs[1].arrow(0, 0, r*np.sin(winkel_s), r*np.cos(winkel_c), head_width=0.05)#Pfeil zeichnen
            points = []
            recs = []
            for i in range(len(phi[t])):
                farbe = "g" if (abs(np.cos(phi[t][i])-np.cos(winkel_c)) <0.1) and (abs(np.sin(phi[t][i])-np.sin(winkel_s)) <0.1) and (r > 0.8) else "r"
                randfarbe = "y" if P[t][i]<0 else "b"
                p = axs[1].plot(np.sin(phi[t][i]), np.cos(phi[t][i]), farbe + "o")#Punkt in Kreis einzeichnen
                points.append(p)
                rec = axs[0].add_patch(patch.Rectangle((posmatrix[i][0], posmatrix[i][1]), 0.5*abs(P[t][i]), 0.5, facecolor=farbe, edgecolor=randfarbe, lw=1.5))#Punkt einzeichnen
                recs.append(rec)
            plt.draw()#Zeichnen
            plt.pause(0.01)#Warten
            for point in points:
                point.pop(0).remove()#Entfernen der Alten Punkte
            for rec in recs:
                rec.remove()
            pfeil.remove()#Entfernen des alten Pfeils


def netz(T, h, K, P, pos):
    """T ist die Simulationsdauer, h die Schrittweite, K die Adjazenzmatrix und P der Lesitungsvektor"""
    phi = np.zeros(((int)(T/h), len(P[0]))) #Phi hat die Form einer Matrix bei der der erste Indize den Zeitpunkt und der zweite den Oszilator angibt
    phipunkt = np.zeros(((int)(T/h), len(P[0]))) #Phipunkt hat die gleiche Form
    phi[0] = np.random.random(len(P[0]))*2*np.pi #Anfangswinkel werden normalverteilt
    phipunkt[0] = 0.1*np.random.standard_cauchy(len(P[0])) #Anfangsfrequenzen sind normalverteilt
    pro = 0 #Prozentzahl des Aktuellen Berechnungsfortschritts
    mess = 0 #Prozentzahl, bei der die letzte Benachrichtigung geprinted wurde
    for i in range(1, int(T/h)):
        phi[i], phipunkt[i] = rk4(kuramoto, [phi[i-1], phipunkt[i-1]], h, K, P[i], P[i+1]) #Berechnen der nächsten Werte
        pro = 100*(i*h/T) #Aktuellen Fortschritt errechnen
        if (mess + 1 < pro): #ueberpruefen, ob eine neue Nachricht geprinted werden soll
            print(str(round(pro, 0))+"% Fortschritt")
            r, _, _ = o(phi[i])#Ordnungsparameter bestimmen
            print (r)
            mess = pro
    return phi


def init(net): 
    p = lambda T, h, P, *args: [P]*(int(T/h + 1)) if(len(np.shape(P))==1) else P
    if(net=='rumänien'):
        ad = np.load('saves/romAdj.npy')#Laden der Adjazenzmatrix
        K=sp.csr_matrix(ad)
        P=np.zeros(K.shape[0])#Leistungen abwechselnd 1 und -1 setzen
        for i in range(K.shape[0]):
            if i%2==0:
                P[i]=-1.0
            else:
                P[i]=1.0
        pos = np.load('saves/romPos.npy')
    elif(net=='n11' or net=='n11_var'):
        K = [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
            [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ]
        P = [-0.1, -0.05, -0.1, -0.05, -0.35, -0.1, -0.2, -0.05, 0.3, 0.5, 0.2]
        pos = [[1.5, 3], [3.5, 3], [4, 2.5], [4, 1.5], [3.5, 1], [1.5, 1], [1, 1.5], [1, 2.5], [2, 2], [3, 2], [2, 4]]
        if(net=='n11_var'): p = lambda T, h, P: berechnep(T, h, P, lambda t, i, P, *a: abs(np.cos(t*i + i))*P[i] if not i==8 and not i==9 else P[i])
    elif(net=='n3'):
        K = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]
        P = [1, -0.3, -0.7]
        pos = [[1, 1], [1, 2], [2, 1.5]]
    elif(net=='nd'):
        K = adj([5], [4], [5, 6], [6, 7], [1, 3, 10], [3, 4, 7, 10], [4, 6, 8, 15], [7], [10], [5, 9, 11], [10, 13, 18], [13], [11, 12, 16], [15, 16], [7, 14], [13, 14, 17], [16], [11])
        P = [1, 1, -1, -1, -0.75, -0.75, 2, -0.75, 1, -1, -1, -1, -0.75, 0.5, 0.5, 0.5, 0.5, 1]
        pos = np.load("saves/dtpos.npy")
        p = lambda T, h, P, *args: berechnep(T, h, P, leistung)
    else: sys.exit('Unknown Net')

    return K, P, pos, p


def adj(*args):
    mat = []
    for arg in args:
        a = [0]*len(args)
        for i in arg:
            a[i-1]=1
        mat.append(a)
    return mat


def calc(net, T, h, skip=-1):
    if(skip==-1): skip=T
    K, P, pos, p = init(net)
    P = p(T, h, P)
    phi = netz(T, h, K, P, pos)#Hauptfunktion mit entsprechenden Werten aufrufen
    np.save('saves/' + net + '.npy', phi)#Speichern
    plotphi(phi, h, T, skip, K, pos, P)#Plotten der Ergebnisse


def plot(net, T, h, skip=-1):
    if(skip==-1): skip=T
    K, P, pos, p = init(net)
    P = p(T, h, P)
    phi = np.load('saves/' + net + '.npy')
    plotphi(phi, h, T, skip, K, pos, P)


def berechnep(T, h, P, f):
    a  = [P]
    for t in range(int(T/h)+1):
        b=[]
        for i in range(len(P)-1):
            b.append(f(t, i, P, h))
        a.append([*b, -sum(b)])
    return a

def leistung(t, i, P, h):
    p = P[i] + 0.3*P[i]*np.sin(t*(1 + 0.01*i))
    if(t*h>120 and (i==0 or i==1)): p*=2
    if(t*h>140 and (i==0 or i==1)): p/=6
    if(t*h>160 and (i==6)): p*=2
    if(t*h>190 and (i==13 or i==14 or i==15 or i==16)): p*=4
    if(t*h>230 and (i==6)): p/=4
    if(t*h>260 and (i==13 or i==14 or i==15 or i==16)): p/=4
    
    return p


plot('nd', 300, 0.01, 300)#Hauptfunktion mit entsprechenden Werten aufrufen
