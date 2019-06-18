
# # Grundmodell
# ## Elektrisches Netzwerk 

# Gruppe A
# Stand: 13.06.2019
# -------------------------------------------------------------------
# 
# Verändert in dieser Version (10.-13.6.):
#    - gezippte Adajazenzmatrix nicht in Kuarmoto sondern in Main berechnet und übergeben
#    - P-Werte an Oszillatoren rangeschrieben
#    - synchronisiert-Variable auf richtige Boolean geändert! (-> Farben-Abfrage für Plot in zeichnen() vereinfacht)
#    - iniLeistung() angepasst je nach N (gerade oder ungerade)
#    - Pfeil beschriftet (r)
#    - Synchronisationsmindestdifferenz abhängig von Schrittzahl
#    - Bilder werden in Ordnern pro durchlauf gespeichert (vlt noch schöner machen?)
#    - Anfangs-Verteilung von theta wird gespeichert
#    - Reihenfolge der Initialisierungswerte vertauscht (insb. N)
#    - Alle Anfangswerte speichern (ergänzt in Main)
#    - Kuramoto-Gleichung wieder mit + gesetzt (damit synchronisiert)
#    - bei Plotten von Verbindungslinie: Dicke mit l multipliziert (Wert der Adjazenzmatrix)
#    - Kopplung angepasst (vgl Paper: K > P0, damit Netz synchronisiert)

# -------------------------------------------------------------------

# ## Wichtige Infos:
# 
# ### theta[x,0]: Phase des Oszillators x (THETA)
# 
# ### theta[x,1]: Ableitung der Phase des Oszillators x (THETA PUNKT)

# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ## Teil 0 - Importieren von Packages und Festsetzen des Namens der Modulation

# In[1]:


#Importiere mehrere Packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
#für Ordnerstruktur
import time 
import os


#Nummer des Durchlaufs oder Name definieren (für Abspeicherung der Bilder)
bildname = "test_tim"
#Ordner erstellen für diesen Namen
#os.mkdir("video/"+bildname)


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ## Teil 1 - Initialisieren und Anfangswerte

# -------------------------------------------------------------------

#Feste Zahlen:

#-------------------------- 
#SIMULATIONSPARAMETER

#Anzahl der Maximal-Schritte der Simulation
anzahlschritte = 150000
#Nach wievielen Schritten ein Bild gespeichert werden soll:
output = 1000

#Schrittweite in der Zeit
dt = 0.01


#-------------------------- 
#MODELLPARAMETER

#Adjazenz-Matrix laden
adjazenzmatrix = np.load("adjamat2_h.npy")

#Positions-Matrix laden
position = np.load("posmat2_h.npy")

#Anzahl der Oszillatoren (berechnen aus Größe der Adjazenzmatrix)
N = len(adjazenzmatrix[0])


#Koppungsstärke (wird noch auf Adjazenzmatrix draufmultipliziert)
Kopplung = 30 #sollte deutlich größer als abs(P) sein, damit das Netz synchronisiert (vgl. Literatur)

#Faktor vor erster Ableitung
alpha = 0.1

#Leistung P(j) (kann auch wie Adjazenzmatrix eingelesen werden)
P = np.array([0.]*N) #np.load("Leistung.npy")



#------------------------------
#WEITERE GRUNDVARIABLEN
#(Nicht selbst einzustellen, wird direkt von System erstellt)


#Array für Thetas erstellen (jeweils mit zwei Einträgen pro Oszillator)
theta = np.ndarray(shape=(N,2))
theta_temp = np.ndarray(shape=(N,2)) #hilfsmatrix um synchronisation zu berechnen

#Synchronisations-Array
synchronisiert = np.array([False]*N)


# -------------------------------------------------------------------

#Initialisierungsfunktion definieren

def initialisieren():
    #Phasen zufällig verteilen und Änderung der Phasen auf Anfangswert 0 setzen
    theta[:,0] = np.random.random(N)*np.pi
    theta[:,1] = 0 # alternativ falls zufällige Anfangsverteilung gewünscht: np.random.standard_cauchy(N)
    
    
    #Synchronisierungsparameter zurücksetzen
    synchronisiert[:] = False
    #Theta_temp (für synchro) zurücksetzen
    theta_temp = theta
    
    #überprüfen ob Leistung der Vorraussetzung für das Modell entspricht
    summe = 0
    for i in range(0,N):
        summe = summe + P[i]
    if(summe != 0):
        print("Leistung überprüfen! ","Summe aller P:",summe)
    else:
        print("Erfolgreich initialisiert.")

# -------------------

#Funktion um Leistung zu verteilen (überflüssig wenn P-Array festgelegt wird!)
def iniLeistung():
    #Leistung festsetzen
    v = 1
    for i in range(0,N):
        P[i] = v
        v = v * -1
    if N%2 != 0:
        P[N-1] = 0
    print(P)
    
    #überprüfen ob Leistung der Vorraussetzung für das Modell entspricht
    summe = 0
    for i in range(0,N):
        summe = summe + P[i]
    if(summe != 0):
        print("Leistung überprüfen! ","Summe aller P:",summe)
    else:
        print("Leistung erfolgreich initialisiert.")


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ## Teil 2 - Formeln und Funktionen

# #### Inhalt:
# 1. Ordnungsparameter
# 2. Kuramoto-Funktion
# 3. Runge-Kutta-Verfahren
# 4. Zeichnungsfunktion

# -------------------------------------------------------------------


#1. Ordnungsparameter bestimmen
def ordnung(N,theta):
    #bertrachte Imaginärteil und Realteil des Parameters gesondert
    tempRe = 0
    tempIm = 0
    
    #Schleife um Summe der Theta zu berechnen
    for i in range(0,N):
        tempRe = tempRe + np.cos(theta[i,0])
        tempIm = tempIm + np.sin(theta[i,0])
        
    #Normieren um Mittelwert zu bekommen:
    tempRe = tempRe / N
    tempIm = tempIm / N
    
    #r und Phi berechnen
    r = (tempIm**2 + tempRe**2)**(1/2)
    
    #Fallunterscheidung um phi richtig zu berechnen
    if (tempIm >= 0):
            phi = np.arccos(tempRe)
    else:
            phi = 2 * np.pi - np.arccos(tempRe)
    
    return r,phi


# -------------------


#Kuramoto-Funktion (auf Basis von Fennas Vorlage)
def kuramoto(theta,power, ziped_adja):
    
    #Zwischenwerte auf 0 setzen
    sum_kop= np.array([0.0]*N)
    res = np.ndarray(shape=(N,2))
    m=0
    j=0
    l=0
    
    #Schleife um die Summen für alle Werte zu bekommen (sum_kop[m] ist Summe der Kopplung für Oszillator m)
    for m,j,l in ziped_adja:
        sum_kop[m] += Kopplung * l * np.sin(theta[j,0] - theta[m,0])
        
    #Funktionswerte ermitteln (in 2d)
    res[:,0]= theta[:,1] #Ableitung von Theta ist Theta-Punkt
    res[:,1]= power - alpha*theta[:,1] + sum_kop #Ableitung von Theta-Punkt entspricht Kuramoto-Gleichung

    return res


# -----------------


#Runge-Kutta-4 Verfahren
def rungekutta4(funktion, x, h, power, ziped_adja):
    #Faktoren berechnen
    k1 = funktion(x, power, ziped_adja)    
    k2 = funktion(x + h/2 * k1, power, ziped_adja)   
    k3 = funktion(x + h/2 * k2, power, ziped_adja)
    k4 = funktion(x + h * k3, power, ziped_adja)
    
    #neuen Funktionswert berechnen (mehrdimensional)
    return x + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)


# -----------------


#Hilfsfunktion zum zeichnen
def zeichne(theta, r, phi, synchronisiert, adjazenzmatrix, path, nummern = True):
    
    #Aufbau des Schaubilds
    fig, axs = plt.subplots(1, 2, figsize=(13, 6)) #erzeugt zwei Bilder nebeneinander
    axs[0].plot() #erzeugt plot im ersten Bild
    axs[1].plot() #erzeugt plot im zweiten Bild

    plt.figure("phases")

    
    #Kreis plotten
    an = np.linspace(0,2*np.pi,100)
    axs[1].plot(np.sin(an),np.cos(an), color="black")

    #Ordnungsparameter
    axs[1].arrow(0,0,r*np.cos(phi),r*np.sin(phi),head_width=0.05, color="black", label="K = "+str(Kopplung))
    axs[1].text(0,0,str(round(r*100,1))+"%")
    axs[1].text(0.5,1,"K = "+str(Kopplung)+", a = "+str(alpha))

    #Verbindungslinien plotten
    m=0
    j=0
    l=0
    for m,j,l in zip(*sp.find(adjazenzmatrix)):
        axs[0].plot(position[[m,j],0],position[[m,j],1], "-", color="black",linewidth=1.05*l)


    #Schleife um alle Punkte (Oszillatoren) zu plotten
    for l in range(0,N):
        #Synchronisierung checken
        color = "red" if synchronisiert[l] else "blue"
        axs[0].plot(position[l,0],position[l,1],"o", color=color)
        if nummern == True:
            axs[0].text(position[l,0],position[l,1],str(int(P[l])), size="large")
        axs[1].plot(np.cos(theta[l,0]),np.sin(theta[l,0]),"o", color=color)
            
    #zeichnen und Abspeichern
    plt.draw()
    fileName = path+"/"+str(i)+".png"
    fig.savefig(fileName,dpi=300)
    plt.show(block=False)

    plt.pause(0.01)


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ## Teil 3 - Hauptzelle

# -------------------------------------------------------------------

#Einfügen falls Leistung zurückgesetzt werden soll:
#iniLeistung()

#---------------------------
#VORBEREITEN

#dann alles andere zurücksetzen
initialisieren()



#Um Bilder in Ordnern zu speichern:
tm = time.gmtime()
path = "video/"+bildname+"/"+str(tm.tm_mon)+"-"+str(tm.tm_mday)+"-"+str(tm.tm_hour+2)+str(tm.tm_min)+str(tm.tm_sec)
os.mkdir(path)
#falls benötigt: Anfangswerte werden gespeichert
np.save(path+"/theta_start.npy",theta)
np.save(path+"/P_start.npy",P)
np.save(path+"/adja_start.npy",adjazenzmatrix)
np.save(path+"/pos_start.npy",position)


#zurücksetzen für Prozentanzeige
zaehler = 0





#dann berechnen
for i in range(0,anzahlschritte):


    #------------------
    #RECHNEN
    
    #um zeit zu verringern: berechne zip-funktion:
    zipd_adja = zip(*sp.find(adjazenzmatrix))
    
    #berechne nun die neuen Theta
    theta = rungekutta4(kuramoto,theta,dt, P,zipd_adja)
    

    if (i%(output/10) == 0): #um Zwischenstand anzuzeigen
        print(zaehler,"%")
        zaehler += 10
    if (i%output == 0): #um Bilder zu plotten
        zaehler = 0
        
        #Ordnungsparameter berechnen
        r, phi = ordnung(N,theta)
        print("Synchro:", round(r,5))
        
        #Synchronisation ermitteln
        for o in range(0,N):
            #überprüfen ob 
            if (abs(theta_temp[o,1]- theta[o,1]) <= (0.00001*output) and r>0.8): #(die synchronisierungsbedingung sollte man sich vlt nochmal anschauen)
                synchronisiert[o] = True
            else:
                synchronisiert[o] = False

        #theta_temp auf theta setzen für den nächsten Vergleich
        theta_temp = theta
        
        #---------------
        #PLOTTEN
        
        #Zeichnen
        zeichne(theta, r, phi, synchronisiert, adjazenzmatrix, path)