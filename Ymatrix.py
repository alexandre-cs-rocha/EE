import numpy as np
import pandas as pd
import scipy as sp
from pathlib import Path
from dss import DSS as dss_engine

def reducao_kron(Yij: complex, Ykk: complex):
    yprim = np.zeros((3, 3), dtype=np.complex128)
    for linha in range(3):
        for coluna in range(3):
            if linha == coluna:
                yprim[linha, coluna] = (Yij) - ((-Yij)**2) / (Ykk)
            else:
                yprim[linha, coluna] =  - ((-Yij)**2) / (Ykk)
    return yprim

def Ymatrix(DSSCircuit):
    baseY1 = 0.11575495704526467
    baseY2 = 1.0401257396449703
    #Linha 1
    DSSCircuit.SetActiveElement('line1')
    Zvetor = (DSSCircuit.Lines.Rmatrix + DSSCircuit.Lines.Xmatrix*1j) * DSSCircuit.Lines.Length
    Zmatrix = np.zeros((3, 3), dtype=np.complex128)
    i = 0
    for p in range(3):
        for m in range(3):
            Zmatrix[p, m] = Zvetor[i]
            i += 1
    Yl1 = np.linalg.inv(Zmatrix)
    Ypriml11 = np.concatenate([Yl1, -Yl1])
    Ypriml12 = np.concatenate([-Yl1, Yl1])
    Ypriml1 = np.concatenate([Ypriml11, Ypriml12], axis=1) / baseY1
    
    #Linha 2
    DSSCircuit.SetActiveElement('line2')
    Zvetor = (DSSCircuit.Lines.Rmatrix + DSSCircuit.Lines.Xmatrix*1j) * DSSCircuit.Lines.Length
    Zmatrix = np.zeros((3, 3), dtype=np.complex128)
    i = 0
    for p in range(3):
        for m in range(3):
            Zmatrix[p, m] = Zvetor[i]
            i += 1
    Yl2 = np.linalg.inv(Zmatrix)
    Ypriml21 = np.concatenate([Yl2, -Yl2])
    Ypriml22 = np.concatenate([-Yl2, Yl2])
    Ypriml2 = np.concatenate([Ypriml21, Ypriml22], axis=1) / baseY2

    #Carga
    Yij = (DSSCircuit.Loads.kW - DSSCircuit.Loads.kvar*1j)*1000 / ((DSSCircuit.Loads.kV*1000)**2)
    Yprimc = reducao_kron(Yij, 3*Yij) / baseY2

    #Transformador
    xl = (((4.16*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.06j
    rl = (((4.16*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.005
    xh = (((12.47*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.06j
    rh = (((12.47*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.005
    z1 = (rl+xl)
    z2 = (rh+xh)
    yt = z2/(z1*z2)
    Yprimt = np.zeros((6, 6), np.complex128)
    for i in range(6):
        for j in range(6):
            if i == j:
                Yprimt[i, j] = yt
            elif i == j-3 or j == i-3:
                Yprimt[i, j] = yt
                   
    #Transformador por redução de kron
    '''Ypp = reducao_kron(0.1042837451-0.6257024898j, 0.3128512353-1.877107489j)
    Ysp = reducao_kron(-0.3126005532+1.875603319j, -0.9378016595+5.626809957j)
    Yps = reducao_kron(-0.3126005532+1.875603319j, -0.9378016595+5.626809957j)
    Yss = reducao_kron(0.9378016595+5.626809957j, 2.811150648-16.86690458j)

    Yprimt = np.concatenate([Ypp, Ysp, Yps, Yss])
    Yprimt1 = np.concatenate([Ypp, Ysp])
    Yprimt2 = np.concatenate([Yps, Yss])
    Yprimt = np.concatenate([Yprimt1, Yprimt2], axis=1)'''
    
    #Vsource
    Yvector = [ [343.2684315-1237.766042j , 31.32731086+9.998440971j  , 31.32731086+9.998440971j],
                [31.32731086+9.998440971j  , 343.2684315-1237.766042j , 31.32731086+9.998440971j],
                [31.32731086+9.998440971j  , 31.32731086+9.998440971j  , 343.2684315-1237.766042j]]
    Yprimv = np.zeros((3, 3), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            Yprimv[i, j] = Yvector[i][j] / baseY1

    Ybus = np.zeros((12, 12), np.complex128)
    for i in range(DSSCircuit.NumBuses*3):
        for j in range(DSSCircuit.NumBuses*3):
            if i < 3 and j < 3:
                Ybus[i, j] += Yprimv[i, j]
            
            if i < 6 and j < 6:
                Ybus[i, j] += Ypriml1[i, j]
            
            if 9 > i > 2 and 9 > j > 2:
                Ybus[i, j] += Yprimt[i-3, j-3]
                 
            if i >= 6 and j >= 6:
                Ybus[i, j] += Ypriml2[i-6, j-6]
            
            if i >= 9 and j >= 9:
                Ybus[i, j] += Yprimc[i-9, j-9]
                
    return Ybus
