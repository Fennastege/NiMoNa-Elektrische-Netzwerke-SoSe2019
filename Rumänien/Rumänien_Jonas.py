import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp

def kuramoto(x, args):
    """Kuramoto Gleichungen mit Reibungstherm. phasen sind die Phasen, frequenz die Frequenzen und args=[K, P],
    wobei K die Adjazenzmatrix ist und P die Leistungen.
    Der Rückgabewert sind die Differenzialgleichung für Phasen und Frequenz"""
    phasen, frequenz = x
    K, P = args
    summe = np.zeros(len(phasen))#Initialisieren des Summenvektors
    for m,j,l in zip(*sp.find(K)):#Fuer alle Eintraege von K die nicht 0 sind
        summe[m] +=5.0*l*np.sin(phasen[j]-phasen[m]) #Auf die Summe addieren
    return np.array([frequenz, P - 0.1*frequenz + summe]) #Die Differentialgleichungen zurückgeben  


def rk4(f, x, h, args=[]):
    """Runge-Kutta 4 Verfahren für Differenzialgleichungen in f mit Vektoriellen x-Werten,
    Schrittweite h und eventuell weitere für die Funktion wichtige Parameter als array in args übergeben"""
    k1 = f(x, args)
    k2 = f(x + h/2*k1, args)
    k3 = f(x + h/2*k2, args)
    k4 = f(x + h*k3, args)
    return x + h/6*(k1 + 2*k2 + 2*k3 + k4) #Neue Werte nach dem Runge-Kutta 4 Vervahren zurückgeben


def o(p):
    """Bestimmt den Ordnungsparameter, wenn p der Vektor der Phasen ist"""
    N = len(p)
    re = 0
    im = 0
    for i in p:
        re +=np.cos(i)/N
        im += np.sin(i)/N
    r = np.sqrt(re**2 + im**2)
    phi_c= np.arccos(re/r)
    phi_s= np.arcsin(im/r)
    return r, phi_c, phi_s


def plot(phi, h, T, skip):
    """phi ist die Matrix der Phasen zu allen Zeitpunkten, h die Schrittweite, T die Simulationsdauer 
    und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    plt.figure("phases")
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(np.sin(an), np.cos(an)) #Kreis zeichnen
    plt.draw()
    plt.show(block=False)
    for i in range(0, int(T/h)):#ueber die Zeitschritte iterieren
        if(i%(skip)== 0):#Nur plotten, wenn i ein vielfaches von skip ist
            plt.title("t = " + str(round(h*i, 2)) + "s")#ueberschrift mit der Aktuellen Zeit erstellen
            t = plt.plot(np.sin(phi[i]), np.cos(phi[i]), "ro")#Punkte einzeichnen
            r, winkel_c, winkel_s = o(phi[i])#Ordnungsparameter bestimmen
            t2 = plt.arrow(0, 0, r*np.sin(winkel_s), r*np.cos(winkel_c), head_width=0.05)#Pfeil zeichnen
            plt.draw()#Zeichnen
            plt.pause(0.01)#Warten
            t.pop(0).remove()#Entfernen der Alten Punkte
            t2.remove()#Entfernen des alten Pheils


def netz(T, h, message, K, P, skip):
    """T ist die Simulationsdauer, h die Schrittweite, message ist abstand in dem der Fortschritt ausgegeben wird, K die Adjazenzmatrix, 
    P der Lesitungsvektor und skip gibt an jeder wievielte Zeitschritt gezeichnet werden soll"""
    phi = np.zeros(((int)(T/h), len(P))) #Phi hat die Form einer Matrix bei der der erste Indize den Zeitpunkt und der zweite den Oszilator angibt
    phipunkt = np.zeros(((int)(T/h), len(P))) #Phipunkt hat die gleiche Form
    phi[0] = np.random.random(len(P))*2*np.pi #Anfangswinkel werden normalverteilt
    phipunkt[0] = 0.1*np.random.standard_cauchy(len(P)) #Anfangsfrequenzen sind normalverteilt
    pro = 0 #Prozentzahl des Aktuellen Berechnungsfortschritts
    mess = 0 #Prozentzahl, bei der die letzte Benachrichtigung geprinted wurde
    for i in range(1, int(T/h)):
        phi[i], phipunkt[i] = rk4(kuramoto, [phi[i-1], phipunkt[i-1]], h, [K, P]) #Berechnen der nächsten Werte
        pro = 100*(i*h/T) #Aktuellen Fortschritt errechnen
        if (mess + message < pro): #ueberpruefen, ob eine neue Nachricht geprinted werden soll
            print(str(round(pro, 0))+"% Fortschritt")
            r, _, _ = o(phi[i])#Ordnungsparameter bestimmen
            print (r)
            mess = pro
    plot(phi, h, T, skip)#Plotten der Ergebnisse


ad = np.load('romAdj.npy')#Laden der Adjazenzmatrix
K=sp.csr_matrix(ad)
P=np.zeros(K.shape[0])#Leistungen abwechselnd 1 und -1 setzen
sum_power=0.0
for i in range(K.shape[0]):
    if i%2==0:
        P[i]=-1.0
    else:
        P[i]=1.0
    sum_power=sum_power+P[i]
if sum_power!=0.0:
    print('sum over power does not equal 0!')


# K = [[0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1], [1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
#      [0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0],
#      [0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0], [1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
#      [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ]
# p = [-0.1, -0.05, -0.1, -0.05, -0.35, -0.1, -0.2, -0.05, 0.3, 0.5, 0.2]


K = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]
P = [1, -0.3, -0.7]

netz(200, 0.01, 1, K, P, 100)#Hauptfunktion mit entsprechenden Werten aufrufen
