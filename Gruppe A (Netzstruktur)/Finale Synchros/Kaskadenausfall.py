
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
import scipy.sparse as sp#für Ordnerstruktur


#für Ordnerstruktur

def modell_elektrisches_netzwerk(theta_anfang, adja, Leistung, i_Ausfallleitung, j_Ausfallleitung, Anzahl_ausfallgrenze_Leitungen = 5, technische_maximallast = 0.1, Kopplung = 30, genauigkeit_output = 500, synchrogrenze = 0.999, synchronisiert_toleranz = 0.0001):
    '''
    argumente:
        theta_anfang: synchronisierte Anfangsbedingung
        adja: Adjazenzmatrix
        Leistung: P (Matrix)
        Kopplung: kopplungsstärke
        technische_maximallast: Ausfallgrenze eines Kabels
        i_Ausfallleitung, j_Ausfallleitung: Oszillatoren zwischen denen Leitung ausfällt
        Anzahl_ausfallgrenze_Leitungen: ab wann abgebrochen wird und Kaskadenausfall zurückgegeben wird


        genauigkeit_output: wie oft theoretisch Bild/Synchronisation berechnet wird
        synchrogrenze: wann abgebrochen wird, weil netz synchronisiert ist
        synchronisiert_toleranz: Ab wann Oszi-Differenz als synchronisiert angenommen wird

    outputs:
        1. synchronisiert - True oder False
        2. anzahl_ausgefallener_Leitungen
        3. Zeit bis zur Synchronisation oder Kaskadenausfall
        4. Ordnungsparameter r
    '''

    #print(adja)

    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    # ## Teil 1 - Initialisieren und Anfangswerte

    # -------------------------------------------------------------------

    #Feste Zahlen:

    #--------------------------
    #SIMULATIONSPARAMETER

    #Anzahl der Maximal-Schritte der Simulation
    anzahlschritte = 50000
    #Nach wievielen Schritten ein Bild gespeichert werden soll:
    output = genauigkeit_output

    #Schrittweite in der Zeit
    dt = 0.01

    #Synchronisationsgrenze (zwischen 0 und 1 - Abbruchbedingung)
    synchrogrenze = synchrogrenze

    #Toleranz der Synchronisation
    synchronisiert_toleranz = 0.0001


    #--------------------------
    #MODELLPARAMETER

    #Adjazenz-Matrix laden
    adjazenzmatrix = adja


    #Leitung fällt aus!!!
    adjazenzmatrix[i_Ausfallleitung,j_Ausfallleitung] = 0.0
    adjazenzmatrix[j_Ausfallleitung,i_Ausfallleitung] = 0.0

    #Anzahl der Oszillatoren (berechnen aus Größe der Adjazenzmatrix)
    N = len(adjazenzmatrix[0])


    #Koppungsstärke (wird noch auf Adjazenzmatrix draufmultipliziert)
    Kopplung = Kopplung #sollte deutlich größer als abs(P) sein, damit das Netz synchronisiert (vgl. Literatur)

    #Faktor vor erster Ableitung
    alpha = 0.1

    #Leistung P(j) (kann auch wie Adjazenzmatrix eingelesen werden)
    P = Leistung #np.load("Leistung.npy")

    #Kaskadenausfälle: Obergrenze für theta_j - theta_i (Relativ zur Kopplung an der Leitung)
    max_Last_Faktor = 1/Kopplung * technische_maximallast



    #------------------------------
    #WEITERE GRUNDVARIABLEN
    #(Nicht selbst einzustellen, wird direkt von System erstellt)


    #Array für Thetas erstellen (jeweils mit zwei Einträgen pro Oszillator)
    theta = theta_anfang

    synchronisiert = np.array([True]*N)


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
    def last_quotient(theta, m, j, l, max_Last_Faktor):
        LastQuotient = (np.sin(np.abs(theta[m,0] - theta[j,0])) / (max_Last_Faktor * np.abs(l) * Kopplung))
        return LastQuotient

        # b

    #Überprüft die Last aller Leitungen und passt ggf. die Adjazenzmatrix an.
    def lastTest(theta, adj_einzeln, max_Last_Faktor):
        a = 0
        b = 0
        c = 0
        veraendert = 0
        #print("davor",adj_einzeln)
        for a,b,c in zip(*sp.find(adj_einzeln)):
            if last_quotient(theta, a, b, c, max_Last_Faktor) > 1:
                #print(last_quotient(theta, a, b, c, max_Last_Faktor))
                adj_einzeln[a,b] = 0
                #print("schleife",a,b)
                veraendert += 0.5
                #print("ver",veraendert)
                #print("danach",adjazenzmatrix)

        return adj_einzeln, veraendert


    # -------------------------------------------------------------------
    # -------------------------------------------------------------------

    # ## Teil 3 - Hauptzelle

    # -------------------------------------------------------------------

    #zurücksetzen für Prozentanzeige
    zaehler = 0

    anzahl_ausgefallen = 1
    #------------------
    #dann berechnen
    for i in range(0,anzahlschritte):

        #------------------
        #RECHNEN



        adjazenzmatrix, neu_ausgefallen = lastTest(theta, adjazenzmatrix, max_Last_Faktor)

        anzahl_ausgefallen += neu_ausgefallen

        #um zeit zu verringern: berechne zip-funktion:
        zipd_adja = zip(*sp.find(adjazenzmatrix))
        #print("zip", adjazenzmatrix)
        #berechne nun die neuen Theta
        theta = rungekutta4(kuramoto,theta,dt, P,zipd_adja)


        #------------------
        #IN BESTIMMTEN FÄLLEN PASSIERT ETWAS
        if (i%output == (output-2)):
            theta_temp_2 = theta
        elif (i%output == (output-1)):
            theta_temp_1 = theta
        elif (i%output == 0 and i < output*3):
            r, phi = ordnung(N,theta)
            print(zaehler, anzahl_ausgefallen, round(r,5),np.sum(synchronisiert == True))
        elif (i%output == 0 and i >= output*3): #um Bilder zu plotten

            #Ordnungsparameter berechnen
            r, phi = ordnung(N,theta)

            #---------------------------------
            #SYNCHRONISIERUNG
            #überprüfen mit Toleranzgrenze (hinterste Zahl in der Funktion)
            synchronisiert = synchronisierungsabfrage(theta,theta_temp_1,theta_temp_2,synchronisiert_toleranz)

            print(zaehler, anzahl_ausgefallen, round(r,5),np.sum(synchronisiert == True))
            zaehler += 1


            if  r >= synchrogrenze and synchronisiert.all():
                return  True, anzahl_ausgefallen, i*dt, r

            if anzahl_ausgefallen >=  Anzahl_ausfallgrenze_Leitungen:
                return False, anzahl_ausgefallen, i*dt, r

            if zaehler == 25:
                print("timeout")
                return None, anzahl_ausgefallen, i*dt, r

    return None, anzahl_ausgefallen, i*dt, r
