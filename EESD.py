import numpy as np
import pandas as pd
import scipy as sp
from pathlib import Path
from dss import DSS as dss_engine
from Derivadas import *
from Residuos import *
from Ymatrix import *

def achar_index_barra(barras: pd.DataFrame, barra: int) -> int:
    #Retorna o index da barra do monitor ativo
    DSSCircuit.SetActiveElement(DSSMonitors.Element)
    
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
    geracao = []
    for barra in DSSCircuit.AllBusNames:
        #if barra.isdigit(): è possível que o sourcebus e o reg não entrem para a EE
        DSSCircuit.SetActiveBus(barra)
        g = False
        #print(DSSCircuit.Buses.AllPCEatBus)
        for elem in DSSCircuit.Buses.AllPCEatBus:
            if 'Vsource' in elem:
                g = True
        #Base é em fase-neutro
        base = DSSCircuit.Buses.kVBase
        nomes.append(barra)
        bases.append(base)
        geracao.append(g)

    nomes = np.concatenate([nomes[1:], [nomes[0]]])
    bases = np.concatenate([bases[1:], [bases[0]]])
    geracao = np.concatenate([geracao[1:], [geracao[0]]])

    idx = [i for i in range(len(nomes))]
    barras = pd.DataFrame(columns=['nome_barra', 'Bases', 'Fases', 'Inj_pot_at', 'Inj_pot_rat', 'Flux_pot_at', 'Flux_pot_rat', 'Tensao', 'Geracao'],
                          index=idx)
    barras['nome_barra'] = nomes
    barras.loc[idx, 'Bases'] = bases
    barras.loc[idx, 'Geracao'] = geracao

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

def iniciar_medidores(elem_inj_pot: list, elem_flux_pot: list, elem_tensao: list) -> None:
    for i, barra in enumerate(DSSCircuit.AllBusNames):
        DSSCircuit.SetActiveBus(barra)
        for j, elem in enumerate(DSSCircuit.Buses.AllPCEatBus):
            if 'Load' in elem or 'Generator' in elem or 'Vsource' in elem:
                DSSText.Command = f'New Monitor.pqi{i}{j} element={elem}, terminal=1, mode=1, ppolar=no'
                
        elem = DSSCircuit.Buses.AllPDEatBus[0]
        if elem != 'None':
            DSSCircuit.SetActiveElement(elem)
            if DSSCircuit.ActiveCktElement.BusNames[0].split('.')[0] == barra:
                DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=1, mode=32'
                
            elif DSSCircuit.ActiveCktElement.BusNames[1].split('.')[0] == barra:
                DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=2, mode=32'
                
            else:
                print('Deu errado')

    '''#Incia medidores para Injeção de Potência
    for i, elem in enumerate(elem_inj_pot):
        DSSText.Command = f'New Monitor.pqi{i} element={elem}, terminal=1, mode=1'
    
    #Incia medidores para Fluxo de Potência
    for i, elem in enumerate(elem_flux_pot):
        DSSText.Command = f'New Monitor.pqij{i} element={elem}, terminal=1, mode=1'
        
    #Incia medidores para Módulo de Tensão
    for i, elem in enumerate(elem_tensao):
        DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=1, mode=32'''

def medidas(baseva: int) -> pd.DataFrame: 
    barras = indexar_barras()

    num_medidas = 0
    for idx in range(len(DSSCircuit.AllBusNames)):
        barras['Inj_pot_at'][idx] = np.array([0, 0, 0])
        barras['Inj_pot_rat'][idx] = np.array([0, 0, 0])
        num_medidas += 6
                    
    #Amostra e salva os valores dos medidores no sistema
    DSSMonitors.SampleAll()
    DSSMonitors.SaveAll()

    DSSMonitors.First
    for _ in range(DSSMonitors.Count):
        barra = DSSMonitors.Terminal - 1
        index_barra = achar_index_barra(barras, barra)
        fases = pegar_fases()
        barras['Fases'][index_barra] = fases
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
            medidas_at = np.zeros(3)
            medidas_rat = np.zeros(3)
            
            for i, fase in enumerate(fases):
                medidas_at[fase] = matriz_medidas[i*2]
                medidas_rat[fase] = matriz_medidas[i*2+1]
                
            barras['Inj_pot_at'][index_barra] = -medidas_at*1000 /baseva
            barras['Inj_pot_rat'][index_barra] = -medidas_rat*1000 /baseva
              
        elif 'v' in DSSMonitors.Name:
            if type(barras['Tensao'][index_barra]) != np.ndarray:
                medidas = np.zeros(3)

                for i, fase in enumerate(fases):
                    medidas[fase] = matriz_medidas[i]
                    
                basekv = DSSCircuit.Buses.kVBase
                barras['Tensao'][index_barra] = medidas / (basekv*1000)
                num_medidas += 3
        
        DSSMonitors.Next
        
    return barras, num_medidas

def Conserta_Ybus(Ybus):
    DSSCircuit.Transformers.First
    for _ in range(DSSCircuit.Transformers.Count):
        trafo = DSSCircuit.Transformers.Name
        DSSCircuit.SetActiveElement(trafo)
        num_phases = DSSCircuit.ActiveCktElement.NumPhases
        barras_conectadas = DSSCircuit.ActiveCktElement.BusNames
        DSSCircuit.SetActiveBus(barras_conectadas[0])
        basekv1 = DSSCircuit.Buses.kVBase
        DSSCircuit.SetActiveBus(barras_conectadas[1])
        basekv2 = DSSCircuit.Buses.kVBase
        if '.' in barras_conectadas[0] or '.' in barras_conectadas[1]:
            barras_conectadas[0] = barras_conectadas[0].split('.')[0]
            barras_conectadas[1] = barras_conectadas[1].split('.')[0]
        no1 = nodes[f"{barras_conectadas[0]}.{1}"]
        no2 = nodes[f"{barras_conectadas[1]}.{1}"]
        
        if basekv1 > basekv2:
            n = basekv1 / basekv2
            Ybus[no1:no1+num_phases, no2:no2+num_phases] = (Ybus[no1:no1+num_phases, no2:no2+num_phases])/n
            Ybus[no2:no2+num_phases, no1:no1+num_phases] = (Ybus[no2:no2+num_phases, no1:no1+num_phases])*n
        else:
            n = basekv2 / basekv1
            Ybus[no1:no1+num_phases, no2:no2+num_phases] = (Ybus[no1:no1+num_phases, no2:no2+num_phases])*n
            Ybus[no2:no2+num_phases, no1:no1+num_phases] = (Ybus[no2:no2+num_phases, no1:no1+num_phases])/n
            
        DSSCircuit.Transformers.Next

    DSSCircuit.Loads.First
    if DSSCircuit.Loads.IsDelta: #Caso carga em delta
        pass
    else: #Caso carga em estrela
        DSSMonitors.First
        for _ in range(DSSMonitors.Count):
            elemento = DSSMonitors.Element
            
            if 'load' in elemento:
                DSSCircuit.SetActiveElement(elemento)
                Yij = (DSSCircuit.Loads.kW - DSSCircuit.Loads.kvar*1j)*1000 / ((DSSCircuit.Loads.kV*1000)**2)
                barra_correspondente = DSSCircuit.ActiveCktElement.BusNames[0]

                for k in DSSCircuit.Buses.Nodes:
                    no1 = nodes[f"{barra_correspondente}.{k}"]
                    Ybus[no1, no1] -= Yij

                DSSCircuit.Loads.Next

            if DSSMonitors.Next == None: #Critério de parada
                break
        
    DSSCircuit.SetActiveElement('Vsource.source')
    Yprim = DSSCircuit.ActiveCktElement.Yprim
    real = Yprim[::2]
    imag = Yprim[1::2]*1j
    Yprim = real+imag
    Yprim = np.reshape(Yprim, (6, 6))
    Ybus[:3, :3] -= Yprim[:3, :3]

    return Ybus

def Calcula_pesos(barras: pd.DataFrame, num_medidas: int) -> np.array:
    dp = []
    for medidas in barras['Inj_pot_at']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    continue
                dp.append((medida * 0.01) / (3 * 100))
    
    for medidas in barras['Inj_pot_rat']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    continue
                dp.append((medida * 0.01) / (3 * 100))
                
    for medidas in barras['Flux_pot_at']:
        if type(medidas) == list:
            for medida in medidas:
                dp.append((medida * 0.01) / (3 * 100))
                
    for medidas in barras['Flux_pot_rat']:
        if type(medidas) == list:
            for medida in medidas:
                dp.append((medida * 0.01) / (3 * 100))
    
    for medidas in barras['Tensao']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    continue
                dp.append((medida * 0.002) / (3 * 100))
    
    dp = np.array(dp)**-2
    aux = []
    for x in dp:
        if x > 10**10:
            x = 10**10
        aux.append(x)
         
    matriz_pesos = np.diag(aux)
    
    return matriz_pesos

def Calcula_residuo(vet_estados: np.array, baseva: int) -> np.array:
    vetor_residuos = []
    # ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
    # vet_estados_aux = np.concatenate((ang_ref, vet_estados))
    
    residuo_atual = 0
    for idx, medida in enumerate(barras['Inj_pot_at']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_inj_pot_at(vetor_residuos, vet_estados, residuo_atual, idx, DSSCircuit.NumBuses,
                                               barras, baseva, nodes, Ybus)

    for idx, medida in enumerate(barras['Inj_pot_rat']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_inj_pot_rat(vetor_residuos, vet_estados, residuo_atual, idx, DSSCircuit.NumBuses,
                                                barras, baseva, nodes, Ybus)
            
    for idx1, medidas in enumerate(barras['Flux_pot_at']):
        if type(medidas) == list:
            for medida in medidas:
                elemento = medida[0]
                fases = np.where((np.isnan(medida[1]) == False))[0]
                residuo_atual = Residuo_fluxo_pot_at(vetor_residuos, vet_estados, fases, residuo_atual, idx1, elemento,
                                                     baseva, barras, DSSCircuit, nodes, Ybus)
            
    for idx1, medidas in enumerate(barras['Flux_pot_rat']):
        if type(medidas) == list:
            for medida in medidas:
                elemento = medida[0]
                fases = np.where((np.isnan(medida[1]) == False))[0]
                residuo_atual = Residuo_fluxo_pot_rat(vetor_residuos, vet_estados, fases, residuo_atual, idx1, elemento,
                                                      baseva, barras, DSSCircuit, nodes, Ybus)
            
    for idx, medida in enumerate(barras['Tensao']):
        if type(medida) == np.ndarray:
            fases = np.where((np.isnan(medida) == False))[0]
            residuo_atual = Residuo_tensao(vetor_residuos, vet_estados, fases, residuo_atual, idx, barras, DSSCircuit.NumBuses)
        
    return np.array(vetor_residuos)

def Calcula_Jacobiana(barras: pd.DataFrame, vet_estados: np.array, num_medidas: int, baseva: int) -> np.array:
    jacobiana = np.zeros((num_medidas, len(vet_estados)))
    # ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
    # vet_estados_aux = np.concatenate((ang_ref, vet_estados))
    
    medida_atual = 0
    for idx, medida in enumerate(barras['Inj_pot_at']):
        if type(medida) == np.ndarray:
            medida_atual = Derivadas_inj_pot_at(jacobiana, medida_atual, idx, DSSCircuit.NumBuses, vet_estados, barras,
                                                nodes, medida, Ybus, baseva)
    
    for idx, medida in enumerate(barras['Inj_pot_rat']):
        if type(medida) == np.ndarray:
            medida_atual = Derivadas_inj_pot_rat(jacobiana, medida_atual, idx, DSSCircuit.NumBuses, vet_estados, barras,
                                                nodes, medida, Ybus, baseva)
            
    for idx1, medida in enumerate(barras['Flux_pot_at']):
        if type(medida) == list:
            elemento = medida[0][0]
            fases = np.where((np.isnan(medida[0][1]) == False))[0]
            medida_atual = Derivadas_fluxo_pot_at(jacobiana, fases, medida_atual, idx1, elemento, barras, nodes, vet_estados,
                                                  DSSCircuit, Ybus, baseva)
            
    for idx1, medida in enumerate(barras['Flux_pot_rat']):
        if type(medida) == list:
            elemento = medida[0][0]
            fases = np.where((np.isnan(medida[0][1]) == False))[0]
            medida_atual = Derivadas_fluxo_pot_rat(jacobiana, fases, medida_atual, idx1, elemento, barras, nodes, vet_estados,
                                                  DSSCircuit, Ybus, baseva)
            
    for idx, medida in enumerate(barras['Tensao']):
        if type(medida) == np.ndarray:
            medida_atual = Derivadas_tensao(jacobiana, barras, medida_atual, idx, DSSCircuit.NumBuses)
        
    return jacobiana

def EE(barras: pd.DataFrame, vet_estados: np.array, matriz_pesos: np.array, baseva:int, erro_max: float, lim_iter: int) -> np.array:
    k = 0
    delx = 1
    while(np.max(np.abs(delx)) > erro_max and lim_iter > k):

        residuo = Calcula_residuo(vet_estados, baseva)
        print(residuo)

        jacobiana = Calcula_Jacobiana(barras, vet_estados, num_medidas, baseva)

        #Calcula a matriz ganho
        matriz_ganho = jacobiana.T @ matriz_pesos @ jacobiana
        
        #Calcula o outro lado da Equação normal
        seinao = jacobiana.T @ matriz_pesos @ residuo

        delx = np.linalg.inv(matriz_ganho) @ seinao
        
        #Atualiza o vetor de estados
        vet_estados += delx
        
        k += 1
        print(k)
    return vet_estados

#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
MasterFile = CurrentFolder / '6Bus' / 'Caso_teste.dss'

DSSCircuit, DSSText, DSSObj, DSSMonitors = InitializeDSS()
DSSObj = dss_engine
DSSText = DSSObj.Text
DSSCircuit = DSSObj.ActiveCircuit
DSSMonitors = DSSCircuit.Monitors

DSSText.Command = 'Clear'
DSSText.Command = f'Compile {MasterFile}'

elem_inj_pot = []
elem_flux_pot = []
elem_tensao = []

iniciar_medidores(elem_inj_pot, elem_flux_pot, elem_tensao)

DSSText.Command = 'Solve'

baseva =  33.3 * 10**6

barras, num_medidas = medidas(baseva)
nodes = organizar_nodes()

Ybus = sp.sparse.csc_matrix(DSSObj.YMatrix.GetCompressedYMatrix())
Ybus = Conserta_Ybus(Ybus)

#Inicializar o vetor de estados com perfil de tensão neutro
vet_estados = np.zeros(len(barras)*6)
for i in range(len(barras)*3, len(barras)*6):
    vet_estados[i] = 1

ang = np.array([])
tensoes = np.array([])
for barra in DSSCircuit.AllBusNames:
    DSSCircuit.SetActiveBus(barra)
    ang = np.concatenate([ang, DSSCircuit.Buses.puVmagAngle[1::2]*2*np.pi / 360])
    tensoes = np.concatenate([tensoes, DSSCircuit.Buses.puVmagAngle[::2]])

ang = np.concatenate([ang[3:], ang[0:3]])
tensoes = np.concatenate([tensoes[3:], tensoes[0:3]])

gabarito = np.concatenate([ang, tensoes])
teste = gabarito.copy()

vet_estados = teste

matriz_pesos = Calcula_pesos(barras, num_medidas)

vet_estados = EE(barras, vet_estados, matriz_pesos, baseva, 10**-3, 10)

print(gabarito)
print(vet_estados)