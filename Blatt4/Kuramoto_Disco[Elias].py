import numpy as np		# ist so'n zahlen und rechenpaket
import matplotlib.pyplot as plt   # ist plotstuff

##Anfangswerte
N=50 #Anzahl der Oszillatoren
K=50 #Kopplungsstärke (später Matrix)
Nt= 1 #(Länge des Zeitintervalls)
h=0.01 #Schrittweite der Zeit
t0=0 #Startwert der Zeit

###50 Oszillatoren mit Cauchyverteilten Phasen und Eigenfrequenzen generiere
theta=np.random.random(N)*2*np.pi
omega=np.random.standard_cauchy(N)
#print(theta, omega)			 			 
	
def ordng(phase): 
	L=len(phase)	#Ordnungsparameter
	
	re = 0
	im = 0
	i=0
	
	while (i<=(L-1)):
		re = re + np.cos(phase[i])/L
		im = im + np.sin(phase[i])/L
		
		i+=1	
		
	
	r = np.sqrt(re**2 + im**2)
	phi = np.arccos(re)
	return r,phi, L

###Die zu lösende dgl mit ihren Parametern:
def kuramoto(phase, Eigenfrq, arg): #args: K
	r, phi, L = ordng(phase)
	### args ist ein Vektor der weiteren Variablen, der hier mit allgemeiner Länge definiert werden und nur später aufgerufen werden muss. ###
	dtheta=[]
	i=0
	while (i<= (L-1)):
		dtheta.append(Eigenfrq[i]+arg*r*np.sin(phi-phase[i]))    		
		i+=1

	return np.array(dtheta)

#===================================================================================================================
#====================RK4============RK4========================RK4====================RK4===========================
###Butcher Tabelle
#[a1, a2, a3, a4,], [b21, b31, b32, b41, b42, b43], [c1, c2, c3, c4]
table=[[0.,0.,0.,0.,0.],["a-Koeffizienten", 0., 0.5, 0.5, 1],["b2-Koeffizienten",0.5],["b3",0., 1/2], ["b4", 0.,0.,1.]]
clist=["c-Koeffizienten",1/6,1/3,1/3,1/6]

##Butcher Tabelle
#RK4 Koeffizienten
#### falls f explizit von t abhängt, muss für die Faktoren k[i] noch die Zeit t+a[i]*h eingesetzt werden
	
def rk4(f,x,h,y,args):
	k1= f(x,y, args)
	k2= f(x+h*(k1*table[2][1]),y,args)
	k3= f(x+h*(k1*table[3][1]+k2*table[3][2]),y,args)
	k4= f(x+h*(k1*table[4][1]+k2*table[4][2]+k3*table[4][3]),y,args)
	ksum= clist[1]*k1+clist[2]*k2+clist[3]*k3+clist[4]*k4
	return x + h*ksum
	
#====================================================================================================================
#====================================================================================================================
	
	
def Gekoppelt(t, theta, omega, args): #args: K,Nt,h
	i=1
	while (i <= (args[1]/args[2])):
		
		r, phi, L = ordng(theta)
		t += args[2]
		theta=rk4(kuramoto, theta, h, omega, args[0])	
		plot(theta,phi,r,i)
		i += 1	 
#color = list(np.random.choice(range(256), size=3)/256)
		
def plot(phase,mitphas, ordnungs, index):  #diese Funktion in gekoppelt implementieren	
	if(index%1==0):
		arrow=plt.arrow(0,0, ordnungs*np.cos(mitphas), ordnungs*np.sin(mitphas),color = list(np.random.choice(range(256), size=3)/256), head_width=0.05)
		dots=plt.plot(np.cos(phase),np.sin(phase),"go")
		plt.draw()
		plt.pause(0.05)
		dots.pop(0).remove() ###Entfernen alter Punkte der Population dots
		arrow.remove()
		




plt.figure("KURAMOTOINCAPS")
an=np.linspace(0,2*np.pi,100)
plt.plot(np.sin(an),np.cos(an))		 
plt.draw()
plt.show(block=False)
Gekoppelt(t0, theta, omega,[K,Nt,h])
			 
			 
			 
			 
			 
			 
			 
	