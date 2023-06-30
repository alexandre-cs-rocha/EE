import numpy as np

def Calcula_Jacobiana(vet_estados: np.array, num_med: int) -> np.array:
    H = np.zeros((num_med, len(vet_estados)))
    mag1 = vet_estados[1]
    mag2 = vet_estados[2]
    magt = mag1*mag2
    anglediff = -vet_estados[0]
    H[0, 0] = -magt*g*np.sin(anglediff)+magt*b*np.cos(anglediff)
    H[0, 1] = 2*mag1*g-mag2*g*np.cos(anglediff)-mag2*b*np.sin(anglediff)
    H[0, 2] = -mag1*g*np.cos(anglediff)-mag1*b*np.sin(anglediff)
    H[1, 0] = magt*g*np.sin(-anglediff)-magt*b*np.cos(-anglediff)
    H[1, 1] = -mag2*g*np.cos(-anglediff)-mag2*b*np.sin(-anglediff)
    H[1, 2] = 2*mag2*g-mag1*g*np.cos(-anglediff)-mag1*b*np.sin(-anglediff)
    H[2, 0] = magt*b*np.sin(anglediff)+magt*g*np.cos(anglediff)
    H[2, 1] = -2*mag1*(b+bsh)+mag2*b*np.cos(anglediff)-mag2*g*np.sin(anglediff)
    H[2, 2] = mag1*b*np.cos(anglediff)-mag1*g*np.sin(anglediff)
    H[3, 0] = -magt*b*np.sin(-anglediff)-magt*g*np.cos(-anglediff)
    H[3, 1] = mag2*b*np.cos(-anglediff)-mag2*g*np.sin(-anglediff)
    H[3, 2] = -2*mag2*(b+bsh)+mag1*b*np.cos(-anglediff)-mag1*g*np.sin(-anglediff)
    H[4, 0] = 0
    H[4, 1] = 1
    H[4, 2] = 0
    return H

def Calcula_Residuo(vet_med: np.array, vet_estados: np.array) -> np.array:
    mag1 = vet_estados[1]
    mag2 = vet_estados[2]
    magt = mag1*mag2
    anglediff = -vet_estados[0]

    residuos = []
    #Calcula as estimativas das potências ativas
    for medida, mag in zip(vet_med[:2], (mag1, mag2)):
        estimativa = (mag**2)*g-magt*g*np.cos(anglediff)-magt*b*np.sin(anglediff)
        residuos.append(medida - estimativa)
        anglediff = -anglediff
        
    #Calcula as estimativas das potências reativas
    for medida, mag in zip(vet_med[2:4], (mag1, mag2)):
        estimativa = -(mag**2)*(b+bsh)+magt*b*np.cos(anglediff)-magt*g*np.sin(anglediff)
        residuos.append(medida - estimativa)
        anglediff = -anglediff
        
    #Calcula a estimativa da Tensão
    residuos.append(vet_med[4] - vet_estados[1])
        
    return np.array(residuos)

def Calcula_Pesos(num_med: int, des_pad: np.array) -> np.array:
    matriz_pesos = np.zeros((num_med, num_med))

    for i, dp in enumerate(des_pad):
        matriz_pesos[i, i] = 1/(dp**2)
        
    return matriz_pesos

def EE(vet_med: np.array, vet_estados: np.array, matriz_pesos: np.array, erro_max: float) -> np.array:
    delx = 1
    while(np.max(delx) > erro_max):
        jacobiana = Calcula_Jacobiana(vet_estados, len(vet_med))
        
        residuo = Calcula_Residuo(vet_med, vet_estados)

        #Calcula a matriz ganho
        matriz_ganho = np.dot(np.dot(jacobiana.T, matriz_pesos), jacobiana)

        #Calcula o outro lado da Equação normal
        seinao = np.dot(np.dot(jacobiana.T, matriz_pesos), residuo)

        delx = np.dot(np.linalg.inv(matriz_ganho), seinao)
        
        #Atualiza o vetor de estados
        vet_estados[0] += delx[0]
        vet_estados[1:] += delx[1:]
    
    return vet_estados

#Inicializa os valores de modelagem da linha dados na questão
y = 1 / (0.0062 + 0.036j)
g = np.real(y)
b = np.imag(y)
bsh = 0.0104/2

#Inicializa as medidas dadas na questão
vet_med = np.array([1.711, -1.614, 0.2245, -0.1566, 1.016])

vet_estados = np.array([0, 1, 1], dtype=float)

des_pad = np.array([1 for _ in range(5)])

W = Calcula_Pesos(len(vet_med), des_pad)

vet_estados = EE(vet_med, vet_estados, W, 10**-3)
print(vet_estados)