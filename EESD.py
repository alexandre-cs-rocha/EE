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

def pegar_fase() -> int:
    return DSSCircuit.ActiveCktElement.NodeOrder[0] - 1

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
    elementos = random.sample(DSSCircuit.Lines.AllNames, num_med_pot)
    for i, elem in enumerate(elementos):
        DSSText.Command = f'New Monitor.{i} element=Line.{elem}, terminal=1, mode=1'
        
    #Define aleatóriamente os elementos para medir as tensões na barra em que eles estão conectados
    elementos = random.sample(DSSCircuit.AllElementNames, num_med_tensao)
    for i, elem in enumerate(elementos):
        DSSText.Command = f'New Monitor.{num_med_pot + i} element={elem}, terminal=1, mode=32'
    
def medidas(num_med_pot: int) -> np.array:
    #Amostra e salva os valores dos medidores no sistema
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()
    
    vet_med = np.zeros(DSSMonitors.Count)
    DSSMonitors.First
    for i in range(DSSMonitors.Count):
        medida = DSSMonitors.Channel(1)
        #Transforma as medidas de tensão em pu
        if i >= num_med_pot:
            DSSCircuit.SetActiveElement(DSSMonitors.Element)
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
            medida = medida / (DSSCircuit.Buses.kVBase*1000)
        vet_med[i] = medida
        DSSMonitors.Next
    
    return vet_med
   
def Calcula_residuo(vet_med: np.array, vet_estados: np.array, num_med_pot: int, num_med_tensao: int) -> np.array:
    vetor_residuos = []
    
    DSSMonitors.First
    for i in range(DSSMonitors.Count):
        DSSCircuit.SetActiveElement(DSSMonitors.Element)
        DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
        #Acha a fase da medida, 0 é a fase 'a', 1 é a fase 'b' e 2 é a fase 'c'
        fase = pegar_fase()
        index_barra = barras[DSSCircuit.Buses.Name]
        #Salva os valores das tensões estimadas da barra atual
        tensao_estimada = vet_estados[DSSCircuit.NumBuses*3:DSSCircuit.NumBuses*3 + 3]
        ang_estimado = vet_estados[(index_barra*3):(index_barra*3)+3]
        
        #Separa os medidores por modo de medição, modo 1 ele mede potências e modo 32 ele mede tensões
        if DSSMonitors.Mode == 1:
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[1])
            index_barra = barras[DSSCircuit.Buses.Name]
            
            Gs, Bs, Bsh = Ymatrix()

            potencia_estimada = 0
            for m in range(DSSCircuit.Lines.Phases):
                tensao_estimada2 = vet_estados[DSSCircuit.NumBuses*3 + (index_barra*3) + m]
                ang_estimado2 = vet_estados[(index_barra*3) + m]
                #Calcula a potencia com base nas tensões e ângulos estimados
                parte1 = tensao_estimada[fase]*tensao_estimada[m]*(Gs[fase][m]*np.cos(ang_estimado[fase]-ang_estimado[m])+(Bs[fase][m]+Bsh[fase][m])*np.sin(ang_estimado[fase]-ang_estimado[m]))
                parte2 = tensao_estimada2*tensao_estimada[fase]*(Gs[fase][m]*np.cos(ang_estimado[fase]-ang_estimado2) + (Bs[fase][m]*np.sin(ang_estimado[fase]-ang_estimado2)))
                potencia_estimada += parte1 - parte2
            potencia = vet_med[i]
            vetor_residuos.append(potencia - potencia_estimada)
                
        elif DSSMonitors.Mode == 32:
            tensao = vet_med[i]
            vetor_residuos.append(tensao - tensao_estimada[fase])
        DSSMonitors.Next
        
    return np.array(vetor_residuos)
            
def Calcula_Jacobiana(vet_med: np.array, vet_estados: np.array, num_med_pot: int, num_med_tensao: int) -> np.array:
    num_medidas = num_med_tensao + num_med_pot
    jacobiana = np.zeros((num_medidas, len(vet_estados)))
    
    DSSMonitors.First
    delta = 0
    for i in range(len(vet_med)):
        DSSCircuit.SetActiveElement(DSSMonitors.Element)
        DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
        
        fase = pegar_fase()
        index_barra = barras[DSSCircuit.Buses.Name]
        
        #Calcular a derivada das tensoes medidas com relação as outras variáveis de estado
        if DSSMonitors.Mode == 32:
            delta = 0
            if i == DSSMonitors.idx-1:
                delta = 1
        
            j = (DSSCircuit.NumBuses*3) + (index_barra*3) + fase
        
            jacobiana[i][j] = delta

        #Calcular a derivada dos fluxos de potencias ativas medidas com relação as outras variáveis de estado
        elif DSSMonitors.Mode == 1:
            Gs, Bs, Bsh = Ymatrix()
            
            tensao_estimada = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
            ang_estimado = vet_estados[(index_barra*3) + fase]
            num_fases = DSSCircuit.Lines.Phases
            
            #Derivada com relação a tensão na barra inicial
            for m in range(num_fases):
                if m == fase:
                    continue
                ang_estimado2 = vet_estados[(index_barra*3) + m]
                delta = tensao_estimada*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2)+(Bs[fase][m]+Bsh[fase][m])*np.sin(ang_estimado-ang_estimado2))
                jacobiana[i][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta
            
            #Derivada com relação ao ângulo na barra inicial
            for m in range(num_fases):
                if m == fase:
                    continue
                tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
                ang_estimado2 = vet_estados[(index_barra*3) + m]
                delta = tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2) - (Bs[fase][m]+Bsh[fase][m])*np.cos(ang_estimado-ang_estimado2))
                jacobiana[i][(index_barra*3) + m] = delta
            
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[1])
            index_barra = barras[DSSCircuit.Buses.Name]

            #Derivada com relação a tensão na barra final
            for m in range(num_fases):
                ang_estimado2 = vet_estados[(index_barra*3) + m]
                delta = -tensao_estimada*(Gs[fase][m]*np.cos(ang_estimado-ang_estimado2)+Bs[fase][m]*np.sin(ang_estimado-ang_estimado2))
                jacobiana[i][DSSCircuit.NumBuses*3 + (index_barra*3) + m] = delta
            
            #Derivada com relação ao ângulo na barra final
            for m in range(num_fases):
                tensao_estimada2 = vet_estados[(DSSCircuit.NumBuses*3) + (index_barra*3) + fase]
                ang_estimado2 = vet_estados[(index_barra*3) + m]
                delta = -tensao_estimada*tensao_estimada2*(Gs[fase][m]*np.sin(ang_estimado-ang_estimado2) - Bs[fase][m]*np.cos(ang_estimado-ang_estimado2))
                jacobiana[i][(index_barra*3) + m] = delta
                 
        DSSMonitors.Next
        
    return jacobiana

#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
print(CurrentFolder)
MasterFile = CurrentFolder / '13Bus' / 'IEEE13Nodeckt.dss'

#Inicializar variaveis
num_med_tensao = 2
num_med_pot = 2
num_medidas = num_med_pot + num_med_tensao

DSSCircuit, DSSText, DSSObj, DSSMonitors = InitializeDSS()
DSSObj = dss_engine
DSSText = DSSObj.Text
DSSCircuit = DSSObj.ActiveCircuit
DSSMonitors = DSSCircuit.Monitors

DSSText.Command = f'Compile {MasterFile}'

iniciar_medidores(num_med_pot, num_med_tensao)

DSSText.Command = 'Solve'



barras = indexar_barras()

vet_med = medidas(num_med_pot)
print(vet_med)

#Inicializar o vetor de estados com perfil de tensão neutro
vet_estados = np.zeros(DSSCircuit.NumBuses*6)
for i in range(DSSCircuit.NumBuses*3, DSSCircuit.NumBuses*6):
    vet_estados[i] = 1

jacobiana = Calcula_Jacobiana(vet_med, vet_estados, num_med_pot, num_med_tensao)

residuos = Calcula_residuo(vet_med, vet_estados, num_med_pot, num_med_tensao)

print(jacobiana)