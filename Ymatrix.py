import numpy as np
import pandas as pd

def reducao_kron(Yij: complex, Ykk: complex):
    yprim = np.zeros((3, 3), dtype=np.complex128)
    for linha in range(3):
        for coluna in range(3):
            if linha == coluna:
                yprim[linha, coluna] = (Yij) - ((-Yij)**2) / (Ykk)
            else:
                yprim[linha, coluna] =  - ((-Yij)**2) / (Ykk)
    return yprim

def KronRed(matrix, neutral):
    # Matrix: deve ser uma matriz NxN
    # Neutral: deve ser um número inteiro

    # Crie a nova matriz com o mesmo tamanho da matriz original
    n = len(matrix)
    newmatrix = np.zeros((n, n), dtype=np.complex128)

    for r in range(n):  # Loop pelas linhas
        if r != neutral:  # Pule a linha do neutro
            for c in range(n):  # Loop pelas colunas
                if c != neutral:  # Pule a coluna do neutro
                    newmatrix[r, c] = matrix[r, c] - (matrix[r, neutral] * matrix[neutral, c]) / matrix[neutral, neutral]

    # Agora reduza a linha e coluna específica
    # A nova matriz se tornou uma matriz (N-1)x(N-1)
    newmatrix = np.delete(newmatrix, neutral, axis=0)
    newmatrix = np.delete(newmatrix, neutral, axis=1)

    return newmatrix

def Ymatrix(DSSCircuit, baseva):
    baseY1 = baseva / (7.199557856794634*1000)**2
    baseY2 = baseva / (2.4017771198288433*1000)**2
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
    Ypriml1 = np.concatenate([Ypriml11, Ypriml12], axis=1)
    
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
    Ypriml2 = np.concatenate([Ypriml21, Ypriml22], axis=1)

    #Carga
    Yij = (DSSCircuit.Loads.kW - DSSCircuit.Loads.kvar*1j)*1000 / ((DSSCircuit.Loads.kV*1000)**2)
    Yprimc4 = np.array([[0.3120377219-0.1511267663j, 0, 0, -0.3120377219+0.1511267663j], 
              [0, 0.3120377219-0.1511267663j, 0, -0.3120377219+0.1511267663j],
              [0, 0, 0.3120377219-0.1511267663j, -0.3120377219+0.1511267663j], 
              [-0.3120377219+0.1511267663j, -0.3120377219+0.1511267663, -0.3120377219+0.1511267663, 0.9361141018-0.4533807521j]], dtype=np.complex128)
    Yprimc = KronRed(Yprimc4, 3)
    #Yprimc = reducao_kron(Yij, 3*Yij)

    #Transformador
    xl = (((2.4017771198288433*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.06j
    rl = (((2.4017771198288433*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.005
    xh = (((7.199557856794634*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.06j
    rh = (((7.199557856794634*1000)**2) / (DSSCircuit.Transformers.kVA*1000))*0.005
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
            Yprimv[i, j] = -Yvector[i][j]
            
    #Yprimv = Ypriml1
    
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
    
    Ybarra = pd.DataFrame(Ybus)
    Ybarra.to_csv('Ybarra.csv', index=False)
    return Ybus
