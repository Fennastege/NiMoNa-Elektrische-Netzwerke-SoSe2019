import matplotlib.pyplot as plt
import numpy as np

def zufall_omega(N=50):                     #generiert eine Cauchy-Verteilung der Kreisfrequenzen Omega in Form eines arrays
    return(np.random.standard_cauchy(N))
def zufall_theta(N=50):                     #generiert eine zufällige Verteilung von Phasen Theta zwischen 0 und 2Pi
    return(np.random.random(N)*2*np.pi)     #in beiden Fällen ist die Anzahl der Oszillatoren N, standardmäßig 50
def ordnungsparameter(theta):               #berechnet den ordnungsparameter für ein gegebenes Theta aus dem Realteil bzw. 
    N=len(theta)                            #Imaginärteil von 1/N mal die Summe aller e**(i*Theta)
    re=im=0
    for j in theta:
        re=re+np.cos(j)/N
        im=im+np.sin(j)/N
    return(np.sqrt(re**2+im**2), np.arccos(re)) #gibt am Ende r und Phi aus
def f(K, omega, theta):                     #berechnet die erste Ableitung nach der Zeit für arrays aus gegebenen Omega und Theta
    a=[]                                    #bei gegebenem Kopplungsgrad K
    N=len(theta)
    r, phi=ordnungsparameter(theta)
    i=0
    while i<N:
        a.append(omega[i]+K*r*np.sin(phi-theta[i]))    #der i-te Eintrag entspricht der Ableitung des i-ten Oszillators
        i+=1
    return np.array(a)                      #np.array erlaubt es, später über den einfachen Befehl c*a alle Elemente mit einem
                                            #Skalar zu multiplizieren
def rk4(K, omega, theta, dt):               #Berechnung der neuen Phase Theta im nächsten Zeitschritt als array
    k1=f(K, omega, theta)
    k2=f(K, omega, theta+dt/2*k1)
    k3=f(K, omega, theta+dt/2*k2)
    k4=f(K, omega, theta+dt*k3)
    theta=theta+dt/6*(k1+2*k2+2*k3+k4)
    return theta

def ergebnis(K, N_t, dt, omega, theta):     #Schleife zum Berechnen und Plotten von Theta in jedem Zeitschritt
    i=0
    plt.figure("phases")
    an = np.linspace(0, 2*np.pi, 100)
    plt.plot(np.cos(an), np.sin(an), "b")        #Plotten des Kreises in der Farbe blau
    while i<(N_t/dt):
        r, phi=ordnungsparameter(theta)     #Berechnung des Ordnungsparameters zum Plotten des durchschnittlichen Phasenpfeils
        plt.title("t = " + str(round(dt*i, 2)) + "s")
        a=plt.plot(np.cos(theta), np.sin(theta), "ro")  #Plotten von Punkten auf dem Kreis
        p=plt.arrow(0, 0, r*np.cos(phi), r*np.sin(phi), head_width=0.05) #Plotten des Pfeils
        plt.draw()
        plt.show(block=False) 
        plt.pause(0.01)                     #Pausieren des Plottens
        a.pop(0).remove()                   #Entfernen alter Punkte und Pfeile nachdem diese geplottet und angezeigt wurden
        p.remove()                          #Entfernen alter Pfeile  
        theta=rk4(K, omega, theta, dt)      #Berechnung von Theta für den nächsten Zeitschritt
        i+=1                                #Zählvariable
    return

ergebnis(50, 100, 0.01, zufall_omega(), zufall_theta())