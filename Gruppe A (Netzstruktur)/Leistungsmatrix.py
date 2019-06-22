#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:42:06 2019

@author: henri
"""
#packages importieren
import numpy as np
import matplotlib.pyplot as plt
import random

#Funktion zur Zufälligen Verteilung der Leistung
def Leistung(N):
    t = []
    N=N
    i = 0
    w = 0
    
    #Leistung wird zwischen +-a verteilt; taucht unten als e wieder auf
    a=0.5
    
    while i <=N:
        z=0
        my_number = 0
        x = random.random()
        if x < 0.5:
            my_number = 1
        else:
            my_number = -1
        z = random.random()*my_number*a
        
        #Runden der Zahlen auf 3 Nachkommastellen
        z = z*1000
        z = round(z)
        z = z/1000
        t.append(z)
        i=i+1
        w = w + z
    t[1] = t[1]- w
    
    #Überschussleistung zur Kontrolle ausgeben
    print w
    
    #Matrix zur Kontrolle ausgeben
    #print t
    
    return t,w,a
#Erhöhen, um Anzahl der Durchläufe zur Matrixfindung zu erhöhen
b = 500

p=0    

while p < b:
    #Programm: Mit N=63 sollte eine 64 Array herauskommen   
    N=63
    x,o,e = Leistung(N)
    
    #Abweichung von Leistung herausrechnen, damit Pgesamt=0
    x[1] = x[1]- o
    
    #Neuen Wert driekt wieder runden
    x[1] = (x[1]*1000)
    x[1] = round(x[1])
    x[1] = x[1]/1000
    
    #Veränderte Matrix zur kontrolle ausgeben
    #print x
    
    #Abspeichern, nur wenn Abweichung nicht zu groß und veränderter Wert nicht über 1
    #Eventuell k erhöhen, falls nichts gefunden wird
    k = 0.1
    
    #Kleiner als e, um große Unterschiede zu vermeiden
    if (o**2)<k**2 and x[1]**2<e**2:
    #Name ändern, falls Matrix nicht überschrieben werden soll
        np.save('Leistung2',x)
        p=b+1
        print x
        print 'Gefunden!'
        
    else:
        print 'Zu hoher Rest!',p
        p=p+1
        
print 'Fertig'
