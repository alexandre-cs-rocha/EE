import numpy as np
import pandas as pd
import scipy as sp
from pathlib import Path
from dss import DSS as dss_engine
from Derivadas import *

def achar_index_barra(barras: pd.DataFrame, barra: int) -> int:
    #Retorna o index da barra do monitor ativo
    DSSCircuit.SetActiveElement(DSSMonitors.Element)
    
    #Verifica se é o objeto Vsource.source e retorna a primeira barra se for
    if DSSCircuit.ActiveCktElement.Name == 'Vsource.source':
        return 0

    DSSCircuit.SetActiveBus(DSSCircuit.ActiveCktElement.BusNames[barra])
    nome = DSSCircuit.Buses.Name
    
    return barras.index[barras['nome_barra'] == nome].to_list()[0]

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

def organizar_nodes() -> dict:
    nodes = {}
    for i, node in enumerate(DSSCircuit.YNodeOrder):
        nodes[node.lower()] = i
    
    return nodes

def indexar_barras() -> pd.DataFrame:
    #Designa indíces às barras
    nomes = []
    bases = []
    for barra in DSSCircuit.AllBusNames:
        #if barra.isdigit(): è possível que o sourcebus e o reg não entrem para a EE
        DSSCircuit.SetActiveBus(barra)
        base = DSSCircuit.Buses.kVBase
        nomes.append(barra)
        bases.append(base)
    
    idx = [i for i in range(len(nomes))]
    barras = pd.DataFrame(columns=['nome_barra', 'Bases', 'Inj_pot_at', 'Inj_pot_rat', 'Flux_pot_at', 'Flux_pot_rat', 'Tensao'],
                          index=idx)
    barras['nome_barra'] = nomes
    barras.loc[idx, 'Bases'] = bases

    return barras

def pegar_fases() -> np.array:
    fases = DSSCircuit.ActiveCktElement.NodeOrder - 1
    fases = set(fases)
    fases.discard(-1)
    
    return fases

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

def iniciar_medidores(elem_inj_pot: list, elem_flux_pot: list, elem_tensao: list) -> None:
    #Incia medidores para Injeção de Potência
    for i, elem in enumerate(elem_inj_pot):
        DSSText.Command = f'New Monitor.pqi{i} element={elem}, terminal=1, mode=1'
    
    #Incia medidores para Fluxo de Potência
    for i, elem in enumerate(elem_flux_pot):
        DSSText.Command = f'New Monitor.pqij{i} element={elem}, terminal=1, mode=1'
        
    #Incia medidores para Módulo de Tensão
    for i, elem in enumerate(elem_tensao):
        DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=1, mode=32'

def medidas() -> pd.DataFrame:
    barras = indexar_barras()
    
    #Amostra e salva os valores dos medidores no sistema
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()

    num_medidas = 0
    baseva = 2 * 10**6
    DSSMonitors.First
    for _ in range(DSSMonitors.Count):
        index_barra = achar_index_barra(barras, 0)
        fases = pegar_fases()
        matriz_medidas = DSSMonitors.AsMatrix()[0][2:]
        
        if 'pqij' in DSSMonitors.Name:
            if type(barras['Flux_pot_at'][index_barra]) != list and type(barras['Flux_pot_rat'][index_barra]) != list:
                barras['Flux_pot_at'][index_barra] = []
                barras['Flux_pot_rat'][index_barra] = []
                
            elemento = DSSMonitors.Element
            DSSCircuit.ActiveCktElement.BusNames[1]
            medidas_at = np.full([3], np.NaN)
            medidas_rat = np.full([3], np.NaN)
            
            for i, fase in enumerate(fases):
                medidas_at[fase] = matriz_medidas[i*2]*1000 / baseva
                medidas_rat[fase] = matriz_medidas[i*2+1]*1000 / baseva
                num_medidas += 2
                
            barras['Flux_pot_at'][index_barra].append((elemento, medidas_at))
            barras['Flux_pot_rat'][index_barra].append((elemento, medidas_rat))
        
        elif 'pqi' in DSSMonitors.Name:
            medidas_at = np.full([3], np.NaN)
            medidas_rat = np.full([3], np.NaN)
            
            for i, fase in enumerate(fases):
                medidas_at[fase] = matriz_medidas[i*2]
                medidas_rat[fase] = matriz_medidas[i*2+1]
                num_medidas += 2
            
            '''if DSSMonitors.Element == 'vsource.source':
                medidas_at = -medidas_at
                medidas_rat = -medidas_rat'''
                
            barras['Inj_pot_at'][index_barra] = medidas_at*1000/baseva
            barras['Inj_pot_rat'][index_barra] = medidas_rat*1000/baseva
              
        elif 'v' in DSSMonitors.Name:
            if type(barras['Tensao'][index_barra]) != np.ndarray:
                medidas = np.full([3], np.NaN)

                for i, fase in enumerate(fases):
                    medidas[fase] = matriz_medidas[i]
                    num_medidas += 1

                basekv = DSSCircuit.Buses.kVBase
                barras['Tensao'][index_barra] = medidas / (basekv*1000)
        
        DSSMonitors.Next
    
    return barras, num_medidas

def Calcula_pesos(num_medidas: int) -> np.array:
    matriz_pesos = np.zeros((num_medidas, num_medidas))
    for i in range(num_medidas):
        if DSSMonitors.Mode == 1:
            matriz_pesos[i, i] = 1 / 1
        elif DSSMonitors.Mode == 32:
            matriz_pesos[i, i] = 1 / 1
    
    return matriz_pesos

def Residuo_inj_pot_at(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra: int,
                       num_buses: int, barras: pd.DataFrame) -> int:
    
    barra1 = barras['nome_barra'][index_barra]
    basekv = barras['Bases'][index_barra]
    baseY = (2*10**6) / ((basekv*1000)**2)
    i = 0
    for fase in fases:
        inj_pot_est = 0
        no1 = nodes[barra1+f'.{fase+1}']
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]
        for index_barra2 in range(len(barras['nome_barra'])):
            barra2 = barras['nome_barra'][index_barra2]
            for m in range(3):
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2]
                if Yij != 0:
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    inj_pot_est += tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2)+Bs*np.sin(ang_estimado-ang_estimado2))
                    
        inj_pot_med = barras['Inj_pot_at'][index_barra][fase]
        inj_pot_est = tensao_estimada*inj_pot_est
        vetor_residuos.append(inj_pot_med - inj_pot_est)
    
    residuo_atual += 1
        
    return residuo_atual
    
def Residuo_inj_pot_rat(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra: int,
                        num_buses: int, barras: pd.DataFrame) -> int:
    
    basekv = barras['Bases'][index_barra]
    baseY = (2*10**6) / ((basekv*1000)**2)
    barra1 = barras['nome_barra'][index_barra]
    for fase in fases:
        inj_pot_est = 0
        no1 = nodes[barra1+f'.{fase+1}']
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]
        for index_barra2 in range(len(barras['nome_barra'])):
            barra2 = barras['nome_barra'][index_barra2]
            for m in range(3):
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2]
                if Yij != 0:
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    inj_pot_est += tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
        inj_pot_med = barras['Inj_pot_at'][index_barra][fase]
        inj_pot_est = tensao_estimada*inj_pot_est
        vetor_residuos.append(inj_pot_med - inj_pot_est)
    
    residuo_atual += 1
        
    return residuo_atual

def Residuo_tensao(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra: int) -> int:
        tensao_estimada = vet_estados[(DSSCircuit.NumBuses+index_barra)*3:(DSSCircuit.NumBuses+index_barra)*3 + 3]
        
        for fase in fases:
            tensao = barras['Tensao'][index_barra][fase]
            vetor_residuos.append(tensao - tensao_estimada[fase])
            residuo_atual += 1
            
        return residuo_atual

def Residuo_fluxo_pot_at(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra1: int, elemento: str) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = (2*10**6) / ((basekv*1000)**2)
    num_buses = DSSCircuit.NumBuses
    Bshmatrix = np.zeros((3, 3))
    if 'line' in elemento:
        DSSCircuit.SetActiveElement(elemento)
        Gsmatrix, Bsmatrix, Bshmatrix = Shuntmatrix(DSSCircuit, fases)
    barra2 = DSSCircuit.ActiveCktElement.BusNames[1]
    index_barra2 = barras[barras['nome_barra'] == barra2].index.values[0]
    
    tensao_estimada = vet_estados[(DSSCircuit.NumBuses+index_barra1)*3:(DSSCircuit.NumBuses+index_barra1)*3 + 3]
    ang_estimado = vet_estados[(index_barra1*3):(index_barra1*3) + 3]
    
    for fase in fases:
        no1 = nodes[barra1+f'.{fase+1}']

        pot_ativa_estimada = 0
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2]
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m]
            tensao_estimada2 = vet_estados[DSSCircuit.NumBuses*3 + (index_barra2*3) + m]
            ang_estimado2 = vet_estados[(index_barra2*3) + m]
            #Calcula a potencia com base nas tensões e ângulos estimados
            parte1 = tensao_estimada[fase]*tensao_estimada[m]*(Gs*np.cos(ang_estimado[fase]-ang_estimado[m])+(Bs+Bsh)*np.sin(ang_estimado[fase]-ang_estimado[m]))
            parte2 = tensao_estimada2*tensao_estimada[fase]*(Gs*np.cos(ang_estimado[fase]-ang_estimado2) + (Bs*np.sin(ang_estimado[fase]-ang_estimado2)))
            pot_ativa_estimada += (parte1 - parte2)
        potencia_at = barras['Flux_pot_at'][index_barra1][0][1][fase]
        vetor_residuos.append(potencia_at - pot_ativa_estimada)
        residuo_atual += 1
        
    return residuo_atual
  
def Residuo_fluxo_pot_rat(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra1: int, elemento: str) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = (2*10**6) / ((basekv*1000)**2)
    num_buses = DSSCircuit.NumBuses
    Bshmatrix = np.zeros((3, 3))
    if 'line' in elemento:
        DSSCircuit.SetActiveElement(elemento)
        Gsmatrix, Bsmatrix, Bshmatrix = Shuntmatrix(DSSCircuit, fases)
    barra2 = DSSCircuit.ActiveCktElement.BusNames[1]
    index_barra2 = barras[barras['nome_barra'] == barra2].index.values[0]
    
    tensao_estimada = vet_estados[(DSSCircuit.NumBuses+index_barra1)*3:(DSSCircuit.NumBuses+index_barra1)*3 + 3]
    ang_estimado = vet_estados[(index_barra1*3):(index_barra1*3) + 3]
    
    for fase in fases:
        no1 = nodes[barra1+f'.{fase+1}']

        pot_reativa_estimada = 0
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2]
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m]
            tensao_estimada2 = vet_estados[DSSCircuit.NumBuses*3 + (index_barra2*3) + m]
            ang_estimado2 = vet_estados[(index_barra2*3) + m]
            #Calcula a potencia com base nas tensões e ângulos estimados
            parte1 = tensao_estimada[fase]*tensao_estimada[m]*(Gs*np.sin(ang_estimado[fase]-ang_estimado[m])-(Bs+Bsh)*np.cos(ang_estimado[fase]-ang_estimado[m]))
            parte2 = tensao_estimada2*tensao_estimada[fase]*(Gs*np.sin(ang_estimado[fase]-ang_estimado2) - (Bs*np.cos(ang_estimado[fase]-ang_estimado2)))
            pot_reativa_estimada += (parte1 - parte2)
        potencia_rat = barras['Flux_pot_rat'][index_barra1][0][1][fase]
        vetor_residuos.append(potencia_rat - pot_reativa_estimada)
        residuo_atual += 1
        
    return residuo_atual

def Calcula_residuo(vet_estados: np.array) -> np.array:
    vetor_residuos = []
    ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
    vet_estados_aux = np.concatenate((ang_ref, vet_estados))
    
    residuo_atual = 0
    for idx, medida in enumerate(barras['Inj_pot_at']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_inj_pot_at(vetor_residuos, vet_estados_aux, fases, residuo_atual, idx, DSSCircuit.NumBuses, barras)

    for idx, medida in enumerate(barras['Inj_pot_rat']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_inj_pot_rat(vetor_residuos, vet_estados_aux, fases, residuo_atual, idx, DSSCircuit.NumBuses, barras)
            
    for idx1, medidas in enumerate(barras['Flux_pot_at']):
        if type(medidas) == list:
            for medida in medidas:
                elemento = medida[0]
                fases = np.where((np.isnan(medida[1]) == False))[0]
                residuo_atual = Residuo_fluxo_pot_at(vetor_residuos, vet_estados_aux, fases, residuo_atual, idx1, elemento)
            
    for idx1, medidas in enumerate(barras['Flux_pot_rat']):
        if type(medidas) == list:
            for medida in medidas:
                elemento = medida[0]
                fases = np.where((np.isnan(medida[1]) == False))[0]
                residuo_atual = Residuo_fluxo_pot_rat(vetor_residuos, vet_estados_aux, fases, residuo_atual, idx1, elemento)
            
    for idx, medida in enumerate(barras['Tensao']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_tensao(vetor_residuos, vet_estados_aux, fases, residuo_atual, idx)
        
    return np.array(vetor_residuos)

def Calcula_Jacobiana(barras: pd.DataFrame, vet_estados: np.array, num_medidas: int) -> np.array:
    jacobiana = np.zeros((num_medidas, len(vet_estados)))
    ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
    vet_estados_aux = np.concatenate((ang_ref, vet_estados))
    
    medida_atual = 0
    for idx, medida in enumerate(barras['Inj_pot_at']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            medida_atual = Derivadas_inj_pot_at(jacobiana, fases, medida_atual, idx, DSSCircuit.NumBuses, vet_estados_aux, barras,
                                                nodes, medida, Ybus)
    
    for idx, medida in enumerate(barras['Inj_pot_rat']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            medida_atual = Derivadas_inj_pot_rat(jacobiana, fases, medida_atual, idx, DSSCircuit.NumBuses, vet_estados_aux, barras,
                                                nodes, medida, Ybus)
            
    for idx1, medida in enumerate(barras['Flux_pot_at']):
        if type(medida) == list:
            elemento = medida[0][0]
            fases = np.where((np.isnan(medida[0][1]) == False))[0]
            medida_atual = Derivadas_fluxo_pot_at(jacobiana, fases, medida_atual, idx1, elemento, barras, nodes, vet_estados_aux,
                                                  DSSCircuit, Ybus)
            
    for idx1, medida in enumerate(barras['Flux_pot_rat']):
        if type(medida) == list:
            elemento = medida[0][0]
            fases = np.where((np.isnan(medida[0][1]) == False))[0]
            medida_atual = Derivadas_fluxo_pot_rat(jacobiana, fases, medida_atual, idx1, elemento, barras, nodes, vet_estados_aux,
                                                  DSSCircuit, Ybus)
            
    for idx, medida in enumerate(barras['Tensao']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            medida_atual = Derivadas_tensao(jacobiana, fases, medida_atual, idx, DSSCircuit.NumBuses)
        
    return jacobiana

def EE(barras: pd.DataFrame, vet_estados: np.array, matriz_pesos: np.array, erro_max: float, lim_iter: int) -> np.array:
    k = 0
    delx = 1
    while(np.max(np.abs(delx)) > erro_max):
        if k >= lim_iter:
            break
        
        jacobiana = Calcula_Jacobiana(barras, vet_estados, num_medidas)
        '''for i, medida in enumerate(jacobiana):
            print(i+1, medida)'''
        residuo = Calcula_residuo(vet_estados)
        
        #Calcula a matriz ganho
        matriz_ganho = np.dot(np.dot(jacobiana.T, matriz_pesos), jacobiana)
        
        #Calcula o outro lado da Equação normal
        seinao = np.dot(np.dot(jacobiana.T, matriz_pesos), residuo)

        delx = np.dot(np.linalg.inv(matriz_ganho), seinao)
        
        #Atualiza o vetor de estados
        vet_estados += delx
        
        k += 1
        
    print(k)
    return vet_estados


#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
MasterFile = CurrentFolder / '4Bus-YY-Bal' / '4Bus-YY-Bal.DSS'

DSSCircuit, DSSText, DSSObj, DSSMonitors = InitializeDSS()
DSSObj = dss_engine
DSSText = DSSObj.Text
DSSCircuit = DSSObj.ActiveCircuit
DSSMonitors = DSSCircuit.Monitors

DSSText.Command = f'Compile {MasterFile}'

elem_inj_pot = ['Vsource.source', 'Load.load1']
elem_flux_pot = []
elem_tensao = ['Line.line1', 'Line.line2', 'Transformer.t1', 'Load.load1']

iniciar_medidores(elem_inj_pot, elem_flux_pot, elem_tensao)

DSSText.Command = 'Solve'

barras, num_medidas = medidas()

nodes = organizar_nodes()

Ybus = sp.sparse.csc_matrix(DSSObj.YMatrix.GetCompressedYMatrix())

#Inicializar o vetor de estados com perfil de tensão neutro
vet_estados = np.zeros(len(barras)*6 - 3)
for i in range(len(barras)*3 - 3, len(barras)*6 - 3):
    vet_estados[i] = 1

matriz_pesos = Calcula_pesos(num_medidas)

vet_estados = EE(barras, vet_estados, matriz_pesos, 10**-3, 3)

ang = np.array([])
tensoes = np.array([])
for barra in DSSCircuit.AllBusNames:
    DSSCircuit.SetActiveBus(barra)
    ang = np.concatenate([ang, DSSCircuit.Buses.puVmagAngle[1::2]])
    tensoes = np.concatenate([tensoes, DSSCircuit.Buses.puVmagAngle[::2]])

gabarito = np.concatenate([ang, tensoes])
print(gabarito)
print(np.concatenate((np.array([0, -2*np.pi/3, 2*np.pi/3]), vet_estados)))
