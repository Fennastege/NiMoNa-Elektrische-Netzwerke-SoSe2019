import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import matplotlib.patches as patches

###Butcher Tabelle
#[a1, a2, a3, a4,], [b21, b31, b32, b41, b42, b43], [c1, c2, c3, c4]
table=[[0.,0.,0.,0.,0.],["a-Koeffizienten", 0., 0.5, 0.5, 1],["b2-Koeffizienten",0.5],["b3",0., 1/2], ["b4", 0.,0.,1.]]
clist=["c-Koeffizienten",1/6,1/3,1/3,1/6]

##Butcher Tabelle
#RK4 Koeffizienten
#### falls f explizit von t abhängt, muss für die Faktoren k[i] noch die Zeit t+a[i]*h eingesetzt werden
	
def rk4(f,x,h,args):
	k1= f(x, args)
	k2= f(x+h*(k1*table[2][1]),args)
	k3= f(x+h*(k1*table[3][1]+k2*table[3][2]),args)
	k4= f(x+h*(k1*table[4][1]+k2*table[4][2]+k3*table[4][3]),args)
	ksum= clist[1]*k1+clist[2]*k2+clist[3]*k3+clist[4]*k4
	return x + h*ksum
	
###############################################

def kuramoto(x, args):
    """Kuramoto Gleichungen mit Reibungstherm. phasen sind die Phasen, frequenz die Frequenzen und args=[K, P],
    wobei K die Adjazenzmatrix ist und P die Leistungen.
    Der Rückgabewert sind die Differenzialgleichung für Phasen und Frequenz"""
    phasen, frequenz = x
    K, P = args
    K=sp.csr_matrix(K)
    summe = np.zeros(len(phasen))#Initialisieren des Summenvektors
    for m,j,l in zip(*sp.find(K)):#Fuer alle Eintraege von K die nicht 0 sind
        summe[m] +=5.0*l*np.sin(phasen[j]-phasen[m]) #Auf die Summe addieren
    return np.array([frequenz, P - 0.1*frequenz + summe]) #Die Differentialgleichungen zurückgeben  

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


def netz(phi, phipunkt, pos, T, h, message, K, P):
	"""phi sind die phasen, phipunkt die frequenzen,T ist die Simulationsdauer, h die Schrittweite, message ist abstand in dem der Fortschritt ausgegeben wird, K die Adjazenzmatrix, 
	P der Lesitungsvektor """
	pro = 0 #Prozentzahl des Aktuellen Berechnungsfortschritts
	mess = 0 #Prozentzahl, bei der die letzte Benachrichtigung geprinted wurde
	i=0
	while (i < int(T/h)):
		r, phi_c, phi_s = o(phi)#Ordnungsparameter bestimmen
		phi, phipunkt = rk4(kuramoto, [phi, phipunkt], h, [K, P]) #Berechnen der nächsten Werte
		pro = 100*(i*h/T) #Aktuellen Fortschritt errechnen
		
		if (mess + message < pro): #ueberpruefen, ob eine neue Nachricht geprinted werden soll
			print(str(round(pro, 0))+"% Fortschritt")
			r,_, _ = o(phi)#Ordnungsparameter bestimmen
			print (r)
			mess=pro
		dynamisch(phi, phi_c, phi_s, pos, r, i,i*h)
		i+=1
		
def Plot(pos, adj, synchro):
	#pos=np.load("romPos.npy")
	#adj=np.load("romAdj.npy")
	
	length=len(position)
	#======================================================================================Verbindungen plotten
	adjrow=0
	while(adjrow<length):
		adjcol=0
		while(adjcol<length):
			if(adj[adjrow][adjcol]==1):
				axs[0].plot([pos[adjrow][0],pos[adjcol][0]],[pos[adjrow][1],pos[adjcol][1]],
						 color="black",linewidth=1) #also plot [x1,x2][y1,y2]
			adjcol+=1
		adjrow+=1
	an=np.linspace(0,2*np.pi,100)
	axs[1].plot(np.sin(an),np.cos(an))		 
	#plt.show()
	#axs[1].add_patch(pat.Rectangle((-10,-10),15,15,facecolor="black"))

def dynamisch(phase,mitphas_c,mitphas_s, pos, ordnungs, index,time):  #diese Funktion in gekoppelt implementieren	
	synchro=False
	points=[]
	if(ordnungs>0.8):
		synchro=True
	
	if(index%50==0):
		plt.title("t="+ str(round(time,2))+"s")
		length=len(pos) ###linker plot
		lengthindex=0
		while(lengthindex<length):#Positionen plotten
			if(synchro==True and abs(np.cos(phase[lengthindex])-np.cos(mitphas_c))<0.1 and abs(np.sin(phase[lengthindex])-np.sin(mitphas_s))<0.1):
				axs[0].plot(pos[lengthindex][0],pos[lengthindex][1], "o", markersize=7,color="red")
				p=axs[1].plot(np.cos(phase[lengthindex]),np.sin(phase[lengthindex]),"o", markersize=7,color="red")
				points.append(p)
			else:
				axs[0].plot(pos[lengthindex][0],pos[lengthindex][1], "o",markersize=7,color="blue")
				p=axs[1].plot(np.cos(phase[lengthindex]),np.sin(phase[lengthindex]),"o", markersize=7,color="blue")
				points.append(p)
			lengthindex+=1
		
		#rechter plot
		arrow=axs[1].arrow(0,0, ordnungs*np.cos(mitphas_c), ordnungs*np.sin(mitphas_s),color = list(np.random.choice(range(256), size=3)/256),linewidth=4, head_width=0.02)
		#dots=axs[1].plot(np.cos(phase),np.sin(phase),"ro")
		
		
		
		plt.draw()
		#plt.show()
		plt.pause(0.03)
		for point in points:
			point.pop(0).remove() ###Entfernen alter Punkte der Population dots
		arrow.remove()



#pos=np.load("sim1pos.npy")
#K=np.load("sim1adj.npy")

K = [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0],
      [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ]
P = [-0.1, -0.05, -0.1, -0.05, -0.35, -0.1, -0.2, -0.05, 0.3, 0.5, 0.2]
position=[[1.5,3],[3.5,3],[4,2.5],[4,1.5],[3.5,1],[1.5,1],[1,1.5],[1,2.5],[2,2],[3,2],[2,4]]

#K = [[0, 1, 1], [1, 0, 0], [1, 0, 0]]
#P = [1, -0.3, -0.7]

#Zufällige Anfangsbedingungen für die Oszillatoren theta_i
theta = np.random.random(len(P))*2*np.pi

#Zufällige Anfangsbedingungen der omega_i nach Cauchy-Verteilung
w=0.1*np.random.standard_cauchy(len(P))

#position=[[1,1],[5,1],[3,3]]
fig, axs=plt.subplots(1,2,figsize=(13,6))

print (K)
Plot(position,K,False)
netz(theta, w, position, 100, 0.01, 1, K, P)

		
