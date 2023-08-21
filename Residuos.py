import numpy as np
import pandas as pd
from Derivadas import Shuntmatrix

def Residuo_inj_pot_at(vetor_residuos: np.array, vet_estados: np.array, residuo_atual: int, index_barra: int,
                       num_buses: int, barras: pd.DataFrame, baseva, nodes: dict, Ybus) -> int:
    fases = barras['Fases'][index_barra]
    barra1 = barras['nome_barra'][index_barra]
    basekv = barras['Bases'][index_barra]
    baseY = baseva / ((basekv*1000)**2)
    print(baseY, basekv)
    #baseY = 1
    for fase in range(3):
        inj_pot_est = 0
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]
        if fase in fases:
            no1 = nodes[barra1+f'.{fase+1}']
            for index_barra2 in range(len(barras['nome_barra'])):
                barra2 = barras['nome_barra'][index_barra2]
                fases2 = barras['Fases'][index_barra2]
                for m in range(3):
                    if m in fases2:
                        no2 = nodes[barra2+f'.{m+1}']
                        Yij = Ybus[no1, no2] / baseY
                        if Yij != 0:
                            Gs = np.real(Yij)
                            Bs = np.imag(Yij)
                            tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                            ang_estimado2 = vet_estados[(index_barra2)*3+m]
                            inj_pot_est += tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2)+Bs*np.sin(ang_estimado-ang_estimado2))
            
        inj_pot_med = barras['Inj_pot_at'][index_barra][fase]

        inj_pot_est = tensao_estimada*inj_pot_est

        vetor_residuos.append(np.abs(inj_pot_med) - np.abs(inj_pot_est))
    
    residuo_atual += 1
        
    return residuo_atual
    
def Residuo_inj_pot_rat(vetor_residuos: np.array, vet_estados: np.array, residuo_atual: int, index_barra: int,
                        num_buses: int, barras: pd.DataFrame, baseva, nodes: dict, Ybus) -> int:
    fases = barras['Fases'][index_barra]
    basekv = barras['Bases'][index_barra]
    baseY = baseva / ((basekv*1000)**2)
    #baseY = 1
    barra1 = barras['nome_barra'][index_barra]
    for fase in range(3):
        inj_pot_est = 0
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]
        if fase in fases:
            no1 = nodes[barra1+f'.{fase+1}']
            for index_barra2 in range(len(barras['nome_barra'])):
                barra2 = barras['nome_barra'][index_barra2]
                fases2 = barras['Fases'][index_barra2]
                for m in range(3):
                    if m in fases2:
                        no2 = nodes[barra2+f'.{m+1}']
                        Yij = Ybus[no1, no2] / baseY
                        if Yij != 0:
                            Gs = np.real(Yij)
                            Bs = np.imag(Yij)
                            tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                            ang_estimado2 = vet_estados[(index_barra2)*3+m]
                            inj_pot_est += tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
        inj_pot_med = barras['Inj_pot_at'][index_barra][fase]
        inj_pot_est = tensao_estimada*inj_pot_est
        vetor_residuos.append(np.abs(inj_pot_med) - np.abs(inj_pot_est))
    
    residuo_atual += 1
        
    return residuo_atual

def Residuo_tensao(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra: int,
                   barras: pd.DataFrame, num_barras: int) -> int:
        tensao_estimada = vet_estados[(num_barras+index_barra)*3:(num_barras+index_barra)*3 + 3]
        
        for fase in fases:
            tensao = barras['Tensao'][index_barra][fase]
            vetor_residuos.append(tensao - tensao_estimada[fase])
            residuo_atual += 1
            
        return residuo_atual

def Residuo_fluxo_pot_at(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra1: int,
                         elemento: str, baseva, barras: pd.DataFrame, DSSCircuit, nodes: dict, Ybus) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = baseva / ((basekv*1000)**2)
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
  
def Residuo_fluxo_pot_rat(vetor_residuos: np.array, vet_estados: np.array, fases: np.array, residuo_atual: int, index_barra1: int,
                         elemento: str, baseva, barras: pd.DataFrame, DSSCircuit, nodes: dict, Ybus) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = baseva / ((basekv*1000)**2)
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
