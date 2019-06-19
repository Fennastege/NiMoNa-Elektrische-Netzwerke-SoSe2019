import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img

#Die Matrix wird auch als dtpos.npy abgespeichert. Damit diese funktioniert
#muss allerdings auch die Skalierung der Karte angepasst werden

pos=[[3.3,15.3],        #01
     [9.7,16],          #02
     [7.5,15.5],        #03
     [12.5,13.5],       #04
     [4.8,13.9],        #05
     [7.3,12],          #06
     [11,11.5],         #07
     [13.3,10.3],       #08
     [1.5,10.5],        #09
     [4.9,11.6],        #10
     [4.3,8.4],         #11
     [10.2,4.5],        #12
     [5.7,5.5],         #13
     [10,7.7],          #14
     [10,9.2],          #15
     [8,7.7],           #16
     [8,9.2],           #17
     [5.7,3.3]]         #18

np.save("dtpos.npy",pos)

k=0
while (k<len(pos)):
    plt.plot(pos[k][0],pos[k][1],"o",color="white")
    k=k+1
plt.grid()
plt.imshow(img.imread('deutschland.png'), extent=[-3, 17, 0, 20])
