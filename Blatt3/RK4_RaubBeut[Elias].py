import numpy as np		# ist so'n zahlen und rechenpaket
import matplotlib.pyplot as plt   # ist plotstuff

###Die zu lösende dgl mit ihren Parametern
def f(x, y, a, b):
	return x*(a-b*y)
	
raub=1; beut=1#x und y, die Räuber Beute Parameter
e1 = 2.0; e2 = -0.8; y1 = 0.02; y2 = -0.0002


###

### RK Verfahren

##Anfangswerte
x=1
y=1
t=0
dt=100 ## Zeit über die gelöst wird

###Butcher Tabelle
#[a1, a2, a3, a4,], [b21, b31, b32, b41, b42, b43], [c1, c2, c3, c4]
table=[[0.,0.,0.,0.,0.],["a-Koeffizienten", 0., 0.5, 0.5, 1],["b2-Koeffizienten",0.5],["b3",0., 1/2], ["b4", 0.,0.,1.]]
clist=["c-Koeffizienten",1/6,1/3,1/3,1/6]
##Butcher Tabelle


# Der Schrittabstand des Verfahrens.
print('Enter the steps in time, h')
h= 0.025; sumx=0; sumy=0

##Variablen und Zählindex definieren
i=0
### Listen definieren
time=[]; xlist=[]; ylist=[]

#RK4 Koeffizienten
	#### falls f explizit von t abhängt, muss für die Faktoren k[i] noch die Zeit t+a[i]*h eingesetzt werden
	
def rk(x,y,h, e1,e2,y1,y2):
	k1= f(x,y,e1,y1)
	l1= f(y,x,e2,y2)
	k2= f(x+h*k1*table[2][1], y+h*l1*table[2][1], e1, y1)
	l2= f(y+h*l1*table[2][1], x+h*k1*table[2][1], e2, y2)
	k3= f(x+h*(k1*table[3][1]+k2*table[3][2]), y+h*(l1*table[3][1]+l2*table[3][2]), e1, y1)
	l3= f(y+h*(l1*table[3][1]+l2*table[3][2]), x+h*(k1*table[3][1]+k2*table[3][2]), e2, y2)
	k4= f(x+h*(k1*table[4][1]+k2*table[4][2]+k3*table[4][3]), y+h*(l1*table[4][1]+l2*table[4][2]+l3*table[4][3]), e1,y1)
	l4= f(y+h*(l1*table[4][1]+l2*table[4][2]+l3*table[4][3]), x+h*(k1*table[4][1]+k2*table[4][2]+k3*table[4][3]), e2, y2)
	ksum= clist[1]*k1+clist[2]*k2+clist[3]*k3+clist[4]*k4
	lsum= clist[1]*l1+clist[2]*l2+clist[3]*l3+clist[4]*l4
	return x + h*ksum, y + h*lsum

	

while (t< dt):
	t=t+h
	xlist.append(x)
	ylist.append(y)
	time.append(t)
	x,y=rk(x,y,h, e1,e2,y1,y2)
	#######Parameter und Schleife für Räuber
	#### falls f explizit von t abhängt, muss für die Faktoren k[i] noch die Zeit t+a[i]*h eingesetzt werden
	#a=e1
	#b=y1
	#x=raub
	#y=beut
	
	
	
	#######Parameter und Schleife für Räuber

	
	
	
plt.xlabel("Zeit t")
plt.ylabel("rot: Räuber, blau: Beute")
plt.plot(time,xlist,'b')
plt.plot(time,ylist,'r')
plt.show()
plt.savefig("python.eps")

#with open('xdata.txt', 'w') as xfile:
#	for i in xrange(0, len(xlist)):
#		xfile.write('{:.3f} {:.3f} \n'.format(time[i], xlist[i]))


#def k1(x,y,a,b):
#	return f(x,y,a,b)
#def k2(x,y,a,b,h):
#	return f(x+h*(k1*table[2][1]), y, a, b)
#def k3(x,y,a,b,h):
#	return f(x+h*(k1*table[2][1]+k2*table[3][1]), y, a, b)
#def k4(x,y,a,b,h):
#	f(x+h*(k1*table[2][1]+k2*table[3][1]+k3*table[4][1]), y, a, b)

#def k(x,y,a,b,h, t2,t3,t4 ):
#	f(x+h*(k1*table[2][1]*t2 +k2*table[3][1]*t3+k3*table[4][1]*t4), y, a, b)

