import numpy as np
import random
from pathlib import Path
from dss import DSS as dss_engine

def reshape(matriz) -> np.array:
    #Transforma as matrizes 1x1 e 2x2 em uma matriz 3x3 de acordo com suas fases
    matriz2 = np.zeros((3, 3))
    fase1 = DSSCircuit.ActiveCktElement.NodeOrder[0]-1
    #Diferencia as matrizes pelo seu tamanho total
    if len(matriz) == 1:
        matriz2[fase1][fase1] = matriz[0]
        return matriz2

    fase2 = DSSCircuit.ActiveCktElement.NodeOrder[1]-1
    if len(matriz) == 4:
        matriz2[fase1][fase1] = matriz[0]
        matriz2[fase1][fase2] = matriz[1]
        matriz2[fase2][fase1] = matriz[2]
        matriz2[fase2][fase2] = matriz[3]
    else:
        matriz2 = matriz.reshape(3, 3)
    
    return matriz2

def Ymatrix() -> tuple:
    #Separa a matriz de admitâncias em uma matriz de Condutâncias e Reatâncias
    Ymatrix = 1/((DSSCircuit.Lines.Rmatrix + DSSCircuit.Lines.Xmatrix*1j)*DSSCircuit.Lines.Length)
    Gs = reshape(np.real(Ymatrix))
    Bs = reshape(np.imag(Ymatrix))
    Bsh = reshape((DSSCircuit.Lines.Cmatrix*DSSCircuit.Lines.Length*60*2*np.pi*10**-9)/2)
    
    return Gs, Bs, Bsh

def indexar_barras() -> dict:
    #Designa indíces às barras
    barras = {}
    for i, barra in enumerate(DSSCircuit.AllBusNames):
        barras[barra] = i
        
    return barras

def pegar_fases() -> np.array:
    return DSSCircuit.ActiveCktElement.NodeOrder - 1

def InitializeDSS() -> tuple:
    DSSObj = dss_engine
    flag = DSSObj.Start(0)
    if flag:
        print('OpenDSS COM Interface initialized succesfully.')
        
    else:
        print('OpenDSS COMInterface failed to start.')
        
    #Set up the interface variables - Comunication OpenDSS with Python
    DSSText = DSSObj.Text
    DSSCircuit = DSSObj.ActiveCircuit
    DSSMonitors = DSSCircuit.Monitors
            
    return DSSCircuit, DSSText, DSSObj, DSSMonitors

def desvio_padrao(vet_med: np.array) -> np.array:
    pr = 0.002
    dp = (pr * vet_med)/(3*1000)
    return dp

def iniciar_medidores(num_med_pot: int, num_med_tensao: int) -> None:
    #Define aleatóriamente as linhas para medir o fluxo de potência
    linhas = []
    DSSCircuit.Lines.First
    for _ in range(DSSCircuit.Lines.Count):
        if not (DSSCircuit.Lines.IsSwitch):
            linhas.append(DSSCircuit.Lines.Name)
            DSSCircuit.Lines.Next

    elementos = random.sample(linhas, num_med_pot)
    for i, elem in enumerate(elementos):
        DSSText.Command = f'New Monitor.{i} element=Line.{elem}, terminal=1, mode=1'
        
    #Define aleatóriamente os elementos para medir as tensões na barra em que eles estão conectados
    elementos = random.sample(DSSCircuit.AllElementNames, num_med_tensao)
    for i, elem in enumerate(elementos):
        DSSText.Command = f'New Monitor.{num_med_pot + i} element={elem}, terminal=1, mode=32'
    
def medidas() -> tuple:
    #Amostra e salva os valores dos medidores no sistema
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()
    
    vet_med = np.array([])
    num_medidas = 0
    DSSMonitors.First
    for _ in range(DSSMonitors.Count):
        DSSCircuit.SetActiveElement(DSSMonitors.Element)
        fases = pegar_fases()
        fases = set(fases)
        fases.discard(-1)
        matriz_medidas = DSSMonitors.AsMatrix()[0][2:]
        
        #Arruma as medidas de potência
        if DSSMonitors.Mode == 1:
            medidas = np.full([6], np.nan)
            
            for i, fase in enumerate(fases):
                medidas[fase*2] = matriz_medidas[i*2]
                medidas[(fase*2)+1] = matriz_medidas[(i*2)+1]
                num_medidas += 2
                
        #Arruma as medidas de tensão em pu
        elif DSSMonitors.Mode == 32:
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
            base = DSSCircuit.Buses.kVBase*1000
            medidas = np.full([3], np.nan)
            
            for fase, medida in zip(fases, matriz_medidas):
                medidas[fase] = medida / base
                num_medidas += 1
                
        vet_med = np.append(vet_med, medidas)
        DSSMonitors.Next
    
    return vet_med, num_medidas

def Calcula_pesos(num_medidas: int) -> np.array:
    matriz_pesos = np.zeros((num_medidas, num_medidas))
    for i in range(num_medidas):
        if DSSMonitors.Mode == 1:
            matriz_pesos[i, i] = 1 / 1
        elif DSSMonitors.Mode == 32:
            matriz_pesos[i, i] = 1 / 1
    
    return matriz_pesos
  
def Calcula_residuo(vet_med: np.array, vet_estados: np.array, num_med_pot: int) -> np.array:
    vetor_residuos = []
    
    DSSMonitors.First
    for i in range(DSSMonitors.Count):
        DSSCircuit.SetActiveElement(DSSMonitors.Element)
        DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
        index_barra = barras[DSSCircuit.Buses.Name]
        
        #Salva os valores das tensões estimadas da barra atual
        tensao_estimada = vet_estados[DSSCircuit.NumBuses*3:DSSCircuit.NumBuses*3 + 3]
        ang_estimado = vet_estados[(index_barra*3):(index_barra*3)+3]
        
        #Acha a ordem das fases na barra da medida, 0 é a fase 'a', 1 é a fase 'b' e 2 é a fase 'c'
        fases = pegar_fases()
        fases = set(fases)
        fases.discard(-1)
        for fase in fases:  
            #Separa os medidores por modo de medição, modo 1 ele mede potências e modo 32 ele mede tensões
            if DSSMonitors.Mode == 32:
                tensao = vet_med[(num_med_pot*6) + ((i-num_med_pot)*3) + fase]
                vetor_residuos.append(tensao - tensao_estimada[fase])

            elif DSSMonitors.Mode == 1:
                DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[1])
                index_barra = barras[DSSCircuit.Buses.Name]
                
                Gs, Bs, Bsh = Ymatrix()

                pot_ativa_estimada = 0
                for m in range(len(fases)):
                    tensao_estimada2 = vet_estados[DSSCircuit.NumBuses*3 + (index_barra*3) + m]
                    ang_estimado2 = vet_estados[(index_barra*3) + m]
                    #Calcula a potencia com base nas tensões e ângulos estimados
                    parte1 = tensao_estimada[fase]*tensao_estimada[m]*(Gs[fase][m]*np.cos(ang_estimado[fase]-ang_estimado[m])+(Bs[fase][m]+Bsh[fase][m])*np.sin(ang_estimado[fase]-ang_estimado[m]))
                    parte2 = tensao_estimada2*tensao_estimada[fase]*(Gs[fase][m]*np.cos(ang_estimado[fase]-ang_estimado2) + (Bs[fase][m]*np.sin(ang_estimado[fase]-ang_estimado2)))
                    pot_ativa_estimada += parte1 - parte2
                potencia = vet_med[i*6 + fase*2]
                vetor_residuos.append(potencia - pot_ativa_estimada)
                
                pot_reativa_estimada = 0
                for m in range(len(fases)):
                    tensao_estimada2 = vet_estados[DSSCircuit.NumBuses*3 + (index_barra*3) + m]
                    ang_estimado2 = vet_estados[(index_barra*3) + m]
                    #Calcula a potencia com base nas tensões e ângulos estimados
                    parte1 = tensao_estimada[fase]*tensao_estimada[m]*(Gs[fase][m]*np.sin(ang_estimado[fase]-ang_estimado[m])-(Bs[fase][m]+Bsh[fase][m])*np.cos(ang_estimado[fase]-ang_estimado[m]))
                    parte2 = tensao_estimada2*tensao_estimada[fase]*(Gs[fase][m]*np.sin(ang_estimado[fase]-ang_estimado2) - (Bs[fase][m]*np.cos(ang_estimado[fase]-ang_estimado2)))
                    pot_reativa_estimada += parte1 - parte2
                potencia = vet_med[i*6 + fase*2 + 1]
                vetor_residuos.append(potencia - pot_reativa_estimada)
                
        DSSMonitors.Next
        
    return np.array(vetor_residuos)

def Derivadas_tensao(jacobiana: np.array, index_barra: int, fases: set, i: int, medida: int) -> int:
    for fase in fases:
        delta = 0
        if i == DSSMonitors.idx-1:
            delta = 1

        j = (DSSCircuit.NumBuses*3) + (index_barra*3) + fase

        jacobiana[medida][j] = delta
        medida += 1
    
    return medida
   
def Derivadas_potencia(jacobiana: np.array, index_barra: int, fases: set, i: int, medida: int) -> int:
    for fase in fases:
        Gs, Bs, Bsh = Ymatrix()
        
        tensao_estimada = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
        ang_estimado = vet_estados[(index_barra*3) + fase]
        
        #Derivada do fluxo de Potência ativa com relação a tensão na barra inicial
        for m in fases:
            if m == fase:
                continue
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = tensao_estimada*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2)+(Bs[fase][m]+Bsh[fase][m])*np.sin(ang_estimado-ang_estimado2))
            print(delta)
            jacobiana[medida][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta

        print(jacobiana)
                
        #Derivada do fluxo de Potência ativa com relação ao ângulo na barra inicial
        for m in fases:
            if m == fase:
                continue
            tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2) - (Bs[fase][m]+Bsh[fase][m])*np.cos(ang_estimado-ang_estimado2))
            jacobiana[medida][(index_barra*3) + m] = delta
        
        medida += 1
            
        #Derivada do fluxo de Potência reativa com relação a tensão na barra inicial
        for m in fases:
            if m == fase:
                continue
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = -tensao_estimada*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2)-(Bs[fase][m])*np.cos(ang_estimado-ang_estimado2))
            jacobiana[medida][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta
            
        #Derivada do fluxo de Potência reativa com relação ao ângulo na barra inicial
        for m in fases:
            if m == fase:
                continue
            tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = -tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2) + (Bs[fase][m]+Bsh[fase][m])*np.sin(ang_estimado-ang_estimado2))
            jacobiana[medida][(index_barra*3) + m] = delta
            
        medida -= 1

        DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[1])
        index_barra = barras[DSSCircuit.Buses.Name]

        #Derivada do fluxo de Potência ativa com relação a tensão na barra final
        for m in fases:
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = -tensao_estimada*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2)+Bs[fase][m]*np.sin(ang_estimado-ang_estimado2))
            jacobiana[medida][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta
            
        #Derivada do fluxo de Potência ativa com relação ao ângulo na barra final
        for m in fases:
            tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = -tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2) - Bs[fase][m]*np.cos(ang_estimado-ang_estimado2))
            jacobiana[medida][(index_barra*3) + m] = delta
        
        medida += 1
            
        #Derivada do fluxo de Potência reativa com relação a tensão na barra final
        for m in fases:
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = -tensao_estimada*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2)-Bs[fase][m]*np.cos(ang_estimado-ang_estimado2))
            jacobiana[medida][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta
            
        #Derivada do fluxo de Potência reativa com relação ao ângulo na barra final
        for m in fases:
            tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
            ang_estimado2 = vet_estados[(index_barra*3) + m]
            delta = tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2) + Bs[fase][m]*np.sin(ang_estimado-ang_estimado2))
            jacobiana[medida][(index_barra*3) + m] = delta
            
        medida += 1
    
    return medida
        
def Calcula_Jacobiana(vet_estados: np.array, num_medidas: int) -> np.array:
    jacobiana = np.zeros((num_medidas, len(vet_estados)))
    
    DSSMonitors.First
    medida = 0
    for i in range(DSSMonitors.Count):
        DSSCircuit.SetActiveElement(DSSMonitors.Element)
        DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
        index_barra = barras[DSSCircuit.Buses.Name]
        
        fases = pegar_fases()
        fases = set(fases)
        fases.discard(-1)
        
        if DSSMonitors.Mode == 32:
            #Calcular a derivada das tensoes medidas com relação as outras variáveis de estado
            medida = Derivadas_tensao(jacobiana, index_barra, fases, i, medida)
            
        elif DSSMonitors.Mode == 1:
            #Calcular a derivada dos fluxos de potencias medidas com relação as outras variáveis de estado
            medida = Derivadas_potencia(jacobiana, index_barra, fases, i, medida)
                                 
        DSSMonitors.Next
        
    return jacobiana

#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
MasterFile = CurrentFolder / '13Bus' / 'IEEE13Nodeckt.dss'

#Inicializar variaveis
num_med_tensao = 0
num_med_pot = 1

DSSCircuit, DSSText, DSSObj, DSSMonitors = InitializeDSS()
DSSObj = dss_engine
DSSText = DSSObj.Text
DSSCircuit = DSSObj.ActiveCircuit
DSSMonitors = DSSCircuit.Monitors

DSSText.Command = f'Compile {MasterFile}'

iniciar_medidores(num_med_pot, num_med_tensao)

DSSText.Command = 'Solve'

barras = indexar_barras()

vet_med, num_medidas = medidas()

#Inicializar o vetor de estados com perfil de tensão neutro
vet_estados = np.zeros(DSSCircuit.NumBuses*6)
for i in range(DSSCircuit.NumBuses*3, DSSCircuit.NumBuses*6):
    vet_estados[i] = 1

matriz_pesos = Calcula_pesos(num_medidas)

jacobiana = Calcula_Jacobiana(vet_estados, num_medidas)
print(jacobiana)

residuos = Calcula_residuo(vet_med, vet_estados, num_med_pot)

#Calcula a matriz ganho
matriz_ganho = np.dot(np.dot(jacobiana.T, matriz_pesos), jacobiana)

#Calcula o outro lado da Equação normal
seinao = np.dot(np.dot(jacobiana.T, matriz_pesos), residuos)

delx = np.dot(np.linalg.inv(matriz_ganho), seinao)
