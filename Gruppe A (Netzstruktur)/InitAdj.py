#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 13:03:54 2019

@author: haukebents
"""
import numpy as np


def AdjInit(fileName, connections, output=False):
    #    fileName: Dateiname unter dem abgespeichert werden soll 
    # connections: Array/Tupel mit Verbindungen zwischen Knoten
    #      output: Soll die Funktion gleich ein Array ausgeben, das verwendet werden kann?
    
    N = len(connections)
    Adj = np.full((N, N), 0)
    
    i = 0
    for i in range(0, N):
        j = 0
        for j in range(0, len(connections[i][:])):
            Adj[i][(connections[i][j])-1] = 1
            Adj[(connections[i][j])-1][i] = 1
            j+=1
        i+=1
        
    np.save(fileName, Adj)
    
    if output:
        return Adj
    
