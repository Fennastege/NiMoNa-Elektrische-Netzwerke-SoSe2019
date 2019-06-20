
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

# Update 20.Juni:
#   - Synchronisationsabfrage gebastelt
#   - Berechnung der Gesamtlänge eingebaut
#   -

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
bildname = "test_nrw_adja1"
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
output = 25

#Schrittweite in der Zeit
dt = 0.01

#Synchronisationsgrenze (zwischen 0 und 1 - Abbruchbedingung)
synchrogrenze = 0.999

#Toleranz der Synchronisation
synchronisiert_toleranz = 0.0001

#WAS WIRD UNTERSUCHT: (Kaskaden oder synchronisierung)
kaskadenuntersuchung = True

#-------------------------- 
#MODELLPARAMETER

#Adjazenz-Matrix laden
adjazenzmatrix = np.load("adjamat_nrw_1.npy")

#Positions-Matrix laden
position = np.load("posmat_nrw_64.npy")

#Anzahl der Oszillatoren (berechnen aus Größe der Adjazenzmatrix)
N = len(adjazenzmatrix[0])


#Koppungsstärke (wird noch auf Adjazenzmatrix draufmultipliziert)
Kopplung = 30 #sollte deutlich größer als abs(P) sein, damit das Netz synchronisiert (vgl. Literatur)

#Faktor vor erster Ableitung
alpha = 0.1

#Leistung P(j) (kann auch wie Adjazenzmatrix eingelesen werden)
P = np.array([0.]*N) #np.load("Leistung.npy")

#Kaskadenausfälle: Obergrenze für theta_j - theta_i (Relativ zur Kopplung an der Leitung)
max_Last_Faktor = 1/Kopplung * 0.1

#Kaskadenausfälle: Farben zur Darstellung der Last auf jeder Leitung
heatmap_colors = ["black", "black", "black", "black", "yellow", "orange", "orangered", "red", "red"]


adjazenzmatrix[0,2] = 0
adjazenzmatrix[2,0] = 0

#------------------------------
#WEITERE GRUNDVARIABLEN
#(Nicht selbst einzustellen, wird direkt von System erstellt)


#Array für Thetas erstellen (jeweils mit zwei Einträgen pro Oszillator)
theta = np.ndarray(shape=(N,2))
#falls Kaskadenausfall untersucht werden soll:

#Synchronisations-Array
if kaskadenuntersuchung:
    synchronisiert = np.array([True]*N)
    theta = np.load("synchronisierter_zustand_theta.npy")
    P = np.load("synchronisierter_zustand_leistung.npy")
else:
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
    v = 0.1
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
# 5. Distanzermittlung
# 6. Synchronisationsabfrage
# 7. a last_quotient
#    b lastTest
# 8. ClusterKoeffizient
#   a lokal
#   b global
        
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
    r = np.sqrt(tempIm**2 + tempRe**2)
    
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
        #Berechne Eintrag im Farbenarray über lastTest() und plotte Verbindungslinien:
        lastQuotient = last_quotient(theta, m, j, l)
        last_index = int((len(heatmap_colors)-1)*lastQuotient)
        if np.abs(lastQuotient) <= 1:
            axs[0].plot(position[[m,j],0],position[[m,j],1], "-", color=heatmap_colors[last_index],linewidth=1.05*l)
        else:
            axs[0].plot(position[[m,j],0],position[[m,j],1], "--", color="grey",linewidth=1.05*l)

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
    #5. Distanzberechnung
def distanz(adjmat, posmat):
    
    Einheit = 24/100 # Skala der Karte in km
    N=len(posmat)
    
    Sum = 0
    for i in range(N):
        for j in range(i):
            
            if (adjmat[i,j]==0): #falls Adjazenzmatrix null (kein Kabel)
                a = 0
            else:
                a = np.sqrt((posmat[j,0] - posmat[i,0])**2 + (posmat[j,1] - posmat[i,1])**2)
            Sum = Sum + a
            
    Laenge = Sum * Einheit
    
    print("Gesamtlänge:",round(Laenge,1)," km") #Gesamtlaenge der benötigten Leitung in km
    
    return Laenge
    
# -------------------------------------------------------------------
    #6. Synchronisationsabfrage

def synchronisierungsabfrage(theta_jetzt,theta_alt, theta_alt_alt, synchrobedingung):
    #Differenz berechnen
    delta_theta = abs(theta_jetzt[:,0]-theta_alt[:,0])
    delta_theta_alt = abs(theta_alt[:,0]-theta_alt_alt[:,0])
    N = len(theta[:,0])
    sync = np.array([False]*N)
    
    for o in range(0,N):
            #überprüfen ob 
            if (abs(delta_theta_alt[o] - delta_theta[o]) <= synchrobedingung):
                sync[o] = True
    
    return sync

# -------------------------------------------------------------------
    #7. Kaskadenausfälle
    
    # a 
    
#Gibt Quotient zwischen Last und maximaler Last der Leitung 
#zwischen m-ten und j-ten Oszillator aus, auf eine Nachkommastelle gerundet. (l: Kopplungsgrad der Leitung)
def last_quotient(theta, m, j, l):
    LastQuotient = (np.sin(np.abs(theta[m][0] - theta[j][0])) / (max_Last_Faktor * np.abs(l) * Kopplung))
    return LastQuotient
    
    # b

#Überprüft die Last aller Leitungen und passt ggf. die Adjazenzmatrix an.
def lastTest(theta, adj):
    a = 0
    b = 0
    c = 0
    for a,b,c in zip(*sp.find(adj)):
        if last_quotient(theta, a, b, c) > 1:
            adj[a,b] = 0
            print("GAU:",a,b)
    return adj


# -------------------------------------------------------------------
    #7. ClusterKoeffizient

def clusterLokal(adj,i):
    nachbarn = np.array([])
    #spalte = 0
    anzVerb = 0
    #Nachbarn bestimmen
    for spalte in range(0,len(adj)):
        if (adj[i,spalte] != 0):
            nachbarn = np.append(nachbarn,spalte)
    #Verbindungen der Nachbarn bestimmen
    for m,j,l in zip(*sp.find(adj)):
        if m in nachbarn and j in nachbarn:
            anzVerb = anzVerb + 0.5
    
    if len(nachbarn) < 2:
        return 0
    else:
        return 2 * anzVerb / (len(nachbarn) * (len(nachbarn) - 1))

def clusterGlobal(clusterLokal,adj):
    summe = 0
    for k in range(0,len(adj)):
        summe += clusterLokal(adj,k)
    return summe/len(adj)


# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ## Teil 3 - Hauptzelle

# -------------------------------------------------------------------

#Einfügen falls Leistung zurückgesetzt werden soll:
#iniLeistung()

#---------------------------
#VORBEREITEN

#dann alles andere zurücksetzen
if kaskadenuntersuchung == False:
    initialisieren()
    #iniLeistung()
    #Länge berechnen und abspeichern
    laenge = distanz(adjazenzmatrix,position)
    cluster = clusterGlobal(clusterLokal,adjazenzmatrix)
    
    text = "synchro"
else:
    text = "kaskaden"
#Um Bilder in Ordnern zu speichern:
tm = time.gmtime()
path = "video/"+bildname+"/"+str(tm.tm_mon)+"-"+str(tm.tm_mday)+"-"+str(tm.tm_hour+2)+str(tm.tm_min)+str(tm.tm_sec)+"_"+text
os.mkdir(path)
#falls benötigt: Anfangswerte werden gespeichert
np.save(path+"/theta_start.npy",theta)
np.save(path+"/P_start.npy",P)
np.save(path+"/adja_start.npy",adjazenzmatrix)
np.save(path+"/pos_start.npy",position)


#zurücksetzen für Prozentanzeige
zaehler = 0
#Zurücksetzen für Abbruchbedingung
gespeichert = False

#------------------
#Bild zuerst einmal Plotten

#Ordnungsparameter berechnen
r, phi = ordnung(N,theta)
print("Synchro am Anfang:", round(r,5))

i = -1
#Zeichnen
zeichne(theta, r, phi, synchronisiert, adjazenzmatrix, path, nummern = False)


#------------------
#dann berechnen
for i in range(0,anzahlschritte):


    #------------------
    #RECHNEN
    
    #um zeit zu verringern: berechne zip-funktion:
    zipd_adja = zip(*sp.find(adjazenzmatrix))
    
    if kaskadenuntersuchung:
        #Überprüfe, ob Leitungen ausfallen:
        adjazenzmatrix = lastTest(theta, adjazenzmatrix)
    
    #berechne nun die neuen Theta
    theta = rungekutta4(kuramoto,theta,dt, P,zipd_adja)
    


    #------------------
    #IN BESTIMMTEN FÄLLEN PASSIERT ETWAS:

    if (i%(output/10) == 0): #um Zwischenstand anzuzeigen und r und phi zwischendurch zu berechnen
        print(zaehler,"%")
        zaehler += 10
        
        if kaskadenuntersuchung == False:
            r, phi = ordnung(N,theta)
        

        
    if (i%output == (output-2)):
        theta_temp_2 = theta
    elif (i%output == (output-1)):
        theta_temp_1 = theta
    elif (i%output == 0 and i != 0): #um Bilder zu plotten
        zaehler = 0
        
        #Ordnungsparameter berechnen
        r, phi = ordnung(N,theta)
        print("Synchro:", round(r,5))
        
        #---------------------------------
        #SYNCHRONISIERUNG
        #überprüfen mit Toleranzgrenze (hinterste Zahl in der Funktion)
        synchronisiert = synchronisierungsabfrage(theta,theta_temp_1,theta_temp_2,synchronisiert_toleranz)
        
        #---------------------------------
        #PLOTTEN
        
        #Zeichnen
        zeichne(theta, r, phi, synchronisiert, adjazenzmatrix, path, nummern = False)
        
        #---------------------------------
        #ABBRUCH
        #zum Abspeichern des Synchronisierten Zustandes
        if synchronisiert.all() and kaskadenuntersuchung == False: 
            np.save(path+"/synchronisierter_zustand_theta.npy",theta)
            np.save(path+"/synchronisierter_zustand_leistung.npy",P)
            gespeichert = True
    #Abbruchbedingung, falls genug synchronisiert
    if gespeichert and r >= synchrogrenze and kaskadenuntersuchung == False:
        break 
            

#GESAMTERGEBNISSE SPEICHERN
if kaskadenuntersuchung == False:
    np.savetxt(path+"/ergebnisse.txt",[laenge,i*dt, cluster],header="Gesamtergebnisse \n \n 1. Länge \n 2. Synchronisationsdauer \n 3. Globaler Clusterkoeffizient", comments="")