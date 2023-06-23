import numpy as np
import random
from pathlib import Path
from dss import DSS as dss_engine

def indexar_barras() -> dict:
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
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()
    
    vet_med = np.zeros(DSSMonitors.Count)
    DSSMonitors.First
    for i in range(DSSMonitors.Count):
        medida = DSSMonitors.Channel(1)
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
        tensao_estimada = vet_estados[DSSCircuit.NumBuses + (index_barra*3) + fase]
        ang_estimado = vet_estados[(index_barra*3) + fase]
        if DSSMonitors.Mode == 1:
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[1])
            index_barra = barras[DSSCircuit.Buses.Name]
            Ymatrix = 1/((DSSCircuit.Lines.Rmatrix + DSSCircuit.Lines.Xmatrix*1j)*DSSCircuit.Lines.Length)
            pot = 0
            for j in range(DSSCircuit.Lines.Phases):
                tensao_estimada2 = vet_estados[DSSCircuit.NumBuses + (index_barra*3) + fase]
                ang_estimado2 = vet_estados[(index_barra*3) + fase]
                pot += tensao_estimada*tensao_estimada2
                
        elif DSSMonitors.Mode == 32:
            tensao = vet_med[i]
            vetor_residuos.append(tensao - tensao_estimada)
        DSSMonitors.Next
        
    return np.array(vetor_residuos)
            
def Calcula_Jacobiana(vet_med: np.array, vet_estados: np.array, num_med_pot: int, num_med_tensao: int) -> np.array:
    num_medidas = num_med_tensao + num_med_pot
    jacobiana = np.zeros((num_medidas, len(vet_estados)))
    
    DSSMonitors.First
    for _ in range(DSSMonitors.Count):
        if DSSMonitors.Mode == 32:
            DSSCircuit.SetActiveElement(DSSMonitors.Element)
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
            continue

        elif DSSMonitors.Mode == 1:
            DSSCircuit.SetActiveElement(DSSMonitors.Element)
            DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[0])
            
            tensao1 = DSSCircuit.ActiveBus.puVoltages[0]
            
        DSSMonitors.Next
        
    return jacobiana

#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
MasterFile = CurrentFolder / '13Bus' / 'IEEE13Nodeckt.dss'

#Inicializar variaveis
num_med_tensao = 5
num_med_pot = 5
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

#Inicializar o vetor de estados com perfil de tensão neutro
vet_estados = np.zeros(DSSCircuit.NumBuses*6)
for i in range(DSSCircuit.NumBuses*3):
    vet_estados[i] = 1

#des_pad = desvio_padrao((vet_med * 0.002)/300)

jacobiana = Calcula_Jacobiana(vet_med, vet_estados, num_med_pot, num_med_tensao)

residuos = Calcula_residuo(vet_med, vet_estados, num_med_pot, num_med_tensao)

DSSMonitors.First
DSSCircuit.SetActiveElement(DSSMonitors.Element)
print(1/(DSSCircuit.Lines.Cmatrix*DSSCircuit.Lines.Length))
print(1/(DSSCircuit.Lines.Cmatrix*DSSCircuit.Lines.Length*1j))