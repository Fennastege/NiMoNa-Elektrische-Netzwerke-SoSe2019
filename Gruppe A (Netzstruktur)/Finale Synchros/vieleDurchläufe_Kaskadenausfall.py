'''
Idee des Programms:
-
'''
import numpy as np
# -----------------------------------------------------------------------------
def TestKaskadenVieleAdj(maxlast, Ordnername):
    import scipy.sparse as sp
    import Kaskadenausfall as kask

    Kopplungsgrad = 50
    sync = 0.99


    path = "adja_" + Ordnername
    adj = np.load(path + "/adja_start.npy")
    leistung = np.load(path + "/synchronisierter_zustand_leistung.npy")
    theta = np.load(path + "/synchronisierter_zustand_theta.npy")

    print(path + "/adja_start.npy")


    # Initiieren
    anzahl = int(np.sum(adj==1)/2)
    synchron = np.array([None]*anzahl) # Synchro nach dem Leitung herausgenommen wurde
    ausgefallen = np.array([1]*anzahl)
    zeit = np.array([0.0]*anzahl)
    ordnungsding = np.array([0.0]*anzahl)
    welche_rausgenommen = np.ndarray(shape=(anzahl,2))


    i = 0
    k = 0
    l = 0


    #
    anzahlTrue = 0
    anzahlFalse = 0
    anzahlNone = 0

    for k,l,m in zip(*sp.find(adj)):
        if i >= (anzahl-1):
            exit
        else:
            if k >= l:
                synchron[i],ausgefallen[i],zeit[i],ordnungsding[i] = kask.modell_elektrisches_netzwerk(theta,adja = np.load(path + "/adja_start.npy"), Leistung = leistung, i_Ausfallleitung = k, j_Ausfallleitung = l,technische_maximallast = maxlast, Kopplung = Kopplungsgrad, synchrogrenze=sync)
                #welche_rausgenommen[i,:] = [k,l]
                #print(" Nummer:",i," : ",k,"-",l,"||", synchron[i],ausgefallen[i],zeit[i],round(ordnungsding[i]),5)
                i += 1

    numbers = np.array([np.sum(synchron==True), np.sum(synchron==False), np.sum(synchron==None)])
    return numbers

# -----------------------------------------------------------------------------


Ordnername = ["dezentral", "zentral", "flaechig", "random", "spinne"]
TestMaxLast = np.arange(0.11, 0.11, 0.01)
q=0
#for q in range(0,len(Ordnername)):
for last in TestMaxLast:
    kaskaden = TestKaskadenVieleAdj(last, Ordnername[q])
    np.savetxt(Ordnername[q] + "_Kaskadenergebnisse.txt", kaskaden, fmt='%1.0f', delimiter='\t', header="True   False   None")
