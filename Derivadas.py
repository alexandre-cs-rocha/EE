import numpy as np
import pandas as pd

def Shuntmatrix(DSSCircuit, fases:np.array) -> np.array:
    Rs = DSSCircuit.Lines.Rmatrix * DSSCircuit.Lines.Length
    Xs = DSSCircuit.Lines.Xmatrix*1j * DSSCircuit.Lines.Length
    Ys_vetor = 1 / (Rs+Xs)
    Bsh_vetor = (DSSCircuit.Lines.Cmatrix * DSSCircuit.Lines.Length * (10**-9) * 2*np.pi * 60)*1j
    Ys_matriz = np.zeros((3, 3), dtype=complex)
    Bsh_matriz = np.zeros((3, 3), dtype=complex)
    i = 0
    for p in fases:
        for m in fases:
            Ys_matriz[p, m] = Ys_vetor[i]
            Bsh_matriz[p, m] = Bsh_vetor[i]
            i += 1
            
    Bsh_matriz = np.imag(Bsh_matriz)
    Gs_matriz = np.real(Ys_matriz)
    Bs_matriz = np.imag(Ys_matriz)
    
    return Gs_matriz, Bs_matriz, Bsh_matriz
        
def Derivadas_inj_pot_at(jacobiana: np.array, medida_atual: int, index_barra: int, num_buses: int,
                         vet_estados: np.array, barras: pd.DataFrame, nodes: dict, medidas: np.array, Ybus, baseva) -> int:
    barra1 = barras['nome_barra'][index_barra]
    fases = barras['Fases'][index_barra]
    basekv = barras['Bases'][index_barra]
    baseY = baseva / ((basekv*1000)**2)
    
    for fase in fases:
        no1 = nodes[barra1+f'.{fase+1}']
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]

        #Derivada da injeção de potência ativa com relação as tensões
        for index_barra2 in range(len(barras['nome_barra'])):
            barra2 = barras['nome_barra'][index_barra2]
            fases2 = barras['Fases'][index_barra2]
            # Derivada considerando a barra terminal
            '''if barras['Geracao'][index_barra]==True and index_barra2==5:
                # implementar funcao
                delta = 0
                for m in fases2:
                    no2 = nodes[barra2+f'.{m+1}']
                    Yij = Ybus[no1, no2] / baseY
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    delta += tensao_estimada*(Gs*np.cos(ang_estimado-ang_estimado2)+Bs*np.sin(ang_estimado-ang_estimado2))

                jacobiana[medida_atual][(num_buses+index_barra2)*3+m] = delta
                continue'''

            for m in fases2:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij) 
                Bs = np.imag(Yij) 
                if no1 == no2:
                    delta = ((tensao_estimada**2)*Gs+medidas[fase]) / tensao_estimada

                else:
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    delta = tensao_estimada*(Gs*np.cos(ang_estimado-ang_estimado2)+Bs*np.sin(ang_estimado-ang_estimado2))

                jacobiana[medida_atual][(num_buses+index_barra2)*3+m] = delta
        
        #Derivadas de injeção de potência ativa com relação aos ângulos
        for index_barra2 in range(len(barras['nome_barra'])):
            
            # Derivada considerando a barra terminal
            '''if barras['Geracao'][index_barra]==True and index_barra2==5:
                # implementar funcao
                print('f')
                barra2 = barras['nome_barra'][index_barra2]
                print(barra2)
                fases2 = barras['Fases'][index_barra2]
                delta=0
                for m in fases2:
                    delta=0
                    no2 = nodes[barra2+f'.{m+1}']
                    Yij = Ybus[no1, no2] / baseY
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    delta += tensao_estimada*tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))

                    print(delta, medida_atual, index_barra, index_barra2, m)
                    jacobiana[medida_atual][(index_barra2)*3 + m] = delta
                
                continue'''

            barra2 = barras['nome_barra'][index_barra2]
            fases2 = barras['Fases'][index_barra2]
            for m in fases2:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij) 
                Bs = np.imag(Yij) 
                if no1 == no2:
                    delta = -Bs*(tensao_estimada**2)
                    delta2 = 0
                    for i in range(len(barras['nome_barra'])):
                        barra3 = barras['nome_barra'][i]
                        fases3 = barras['Fases'][i]
                        for n in fases3:
                            no3 = nodes[barra3+f'.{n+1}']
                            Yij = Ybus[no1, no3] / baseY
                            if Yij != 0:
                                tensao_estimada2 = vet_estados[(num_buses+i)*3 + n]
                                ang_estimado2 = vet_estados[(i)*3 + n]
                                Gs = np.real(Yij)
                                Bs = np.imag(Yij)
                                delta2 += tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
                    delta = delta - tensao_estimada*delta2
                else:
                    tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    delta = tensao_estimada*tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
                    
                jacobiana[medida_atual][(index_barra2)*3 + m] = delta
            
        medida_atual += 1
        
    return medida_atual
        
def Derivadas_inj_pot_rat(jacobiana: np.array, medida_atual: int, index_barra: int, num_buses: int,
                         vet_estados: np.array, barras: pd.DataFrame, nodes: dict, medidas: np.array, Ybus, baseva) -> int:
    barra1 = barras['nome_barra'][index_barra]
    fases = barras['Fases'][index_barra]
    basekv = barras['Bases'][index_barra]
    baseY = baseva / ((basekv*1000)**2)
    
    #Derivada da injeção de potência reativa com relação as tensões
    for fase in fases:
        no1 = nodes[barra1+f'.{fase+1}']
        tensao_estimada = vet_estados[(num_buses+index_barra)*3+fase]
        ang_estimado = vet_estados[(index_barra)*3+fase]
        for index_barra2 in range(len(barras['nome_barra'])):
            barra2 = barras['nome_barra'][index_barra2]
            fases2 = barras['Fases'][index_barra2]
            for m in fases2:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij)
                Bs = np.imag(Yij)
                if no1 == no2:
                    delta = (-(tensao_estimada**2)*(Bs)+medidas[fase]) / tensao_estimada
                else:
                    ang_estimado2 = vet_estados[(index_barra2)*3+m]
                    delta = tensao_estimada*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
                    
                jacobiana[medida_atual][(num_buses+index_barra2)*3+m] = delta

        #Derivadas de injeção de potência reativa com relação aos ângulos
        for index_barra2 in range(len(barras['nome_barra'])):
            if index_barra2 != 0:
                barra2 = barras['nome_barra'][index_barra2]
                fases2 = barras['Fases'][index_barra2]
                for m in fases2:
                    no2 = nodes[barra2+f'.{m+1}']
                    Yij = Ybus[no1, no2] / baseY
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    if no1 == no2:
                        medida_at = barras['Inj_pot_at'][index_barra][fase]
                        delta = -Gs*(tensao_estimada**2) + medida_at

                    else:
                        tensao_estimada2 = vet_estados[(num_buses+index_barra2)*3+m]
                        ang_estimado2 = vet_estados[(index_barra2)*3+m]
                        delta = -tensao_estimada*tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2)+Bs*np.sin(ang_estimado-ang_estimado2))
                        
                    jacobiana[medida_atual][(index_barra2)*3+m] = delta
                
        medida_atual += 1
                
    return medida_atual

def Derivadas_tensao(jacobiana: np.array, barras: pd.DataFrame, medida_atual: int, index_barra: int, num_buses: int) -> int:   
    fases = barras['Fases'][index_barra]
        
    for fase in fases:
        jacobiana[medida_atual][(num_buses*3) + (index_barra*3) + fase ] = 1
        medida_atual += 1
    
    return medida_atual

def Derivadas_fluxo_pot_at(jacobiana: np.array, fases: np.array, medida_atual: int, index_barra1: int, elemento: str,
                           barras: pd.DataFrame, nodes: dict, vet_estados: np.array, DSSCircuit, Ybus, baseva) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = baseva / ((basekv*1000)**2)
    num_buses = DSSCircuit.NumBuses
    DSSCircuit.SetActiveElement(elemento)
    Bshmatrix = np.zeros((3, 3))
    if 'line' in elemento:
        Gsmatrix, Bsmatrix, Bshmatrix = Shuntmatrix(DSSCircuit, fases)
    barra2 = DSSCircuit.ActiveCktElement.BusNames[1]
    index_barra2 = barras[barras['nome_barra'] == barra2].index.values[0]
    
    for fase in fases:
        no1 = nodes[barra1+f'.{fase+1}']

        tensao_estimada = vet_estados[(num_buses*3) + (index_barra1*3) + fase]
        ang_estimado = vet_estados[(index_barra1*3) + fase]

        #Derivada do fluxo de Potência ativa com relação a tensão na barra inicial
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2] / baseY
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m] / baseY
            
            if m == fase:
                delta = tensao_estimada*Gs
                for n in fases:
                    no2 = nodes[barra2+f'.{n+1}']
                    Yij = Ybus[no1, no2] / baseY
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    Bsh = Bshmatrix[fase, n] / baseY
                    tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + n]
                    tensao_estimada3 = vet_estados[(num_buses*3) + (index_barra2*3) + n]
                    ang_estimado2 = vet_estados[(index_barra1*3) + n]
                    ang_estimado3 = vet_estados[(index_barra2*3) + n]
                    delta += tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2)+(Bs+Bsh)*np.sin(ang_estimado-ang_estimado2))
                    delta -= tensao_estimada3*(Gs*np.cos(ang_estimado-ang_estimado3)+Bs*np.sin(ang_estimado-ang_estimado3))

            else:
                ang_estimado2 = vet_estados[(index_barra1*3) + m]
                delta = tensao_estimada*(Gs*np.cos(ang_estimado-ang_estimado2) + (Bs+Bsh)*np.sin(ang_estimado-ang_estimado2))
                
            jacobiana[medida_atual][(num_buses+index_barra1)*3 + m - 3] = delta
            
        if index_barra1 != 0:
            #Derivada do fluxo de Potência ativa com relação ao ângulo na barra inicial
            for m in fases:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij)
                Bs = np.imag(Yij)
                Bsh = Bshmatrix[fase, m] / baseY
                if m == fase:
                    delta = -(tensao_estimada**2)*(Bs+Bsh)
                    for n in fases:
                        no2 = nodes[barra2+f'.{n+1}']
                        Yij = Ybus[no1, no2] / baseY
                        Gs = np.real(Yij)
                        Bs = np.imag(Yij)
                        Bsh = Bshmatrix[fase, n] / baseY
                        tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + n]
                        tensao_estimada3 = vet_estados[(num_buses*3) + (index_barra2*3) + n]
                        ang_estimado2 = vet_estados[(index_barra1*3) + n]
                        ang_estimado3 = vet_estados[(index_barra2*3) + n]
                        delta -= tensao_estimada*tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-(Bs+Bsh)*np.cos(ang_estimado-ang_estimado2))
                        delta += tensao_estimada*tensao_estimada3*(Gs*np.sin(ang_estimado-ang_estimado3)-Bs*np.cos(ang_estimado-ang_estimado3))

                else:
                    tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + m]
                    ang_estimado2 = vet_estados[(index_barra1*3) + m]
                    delta = tensao_estimada*tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2) - (Bs+Bsh)*np.cos(ang_estimado-ang_estimado2))
                    
                jacobiana[medida_atual][(index_barra1*3) + m - 3] = delta
            
        #Derivada do fluxo de Potência ativa com relação a tensão na barra final
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2] / baseY
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m] / baseY
            ang_estimado2 = vet_estados[(index_barra2*3) + m]
            delta = -tensao_estimada*(Gs*np.cos(ang_estimado-ang_estimado2) + Bs*np.sin(ang_estimado-ang_estimado2))
            jacobiana[medida_atual][num_buses*3 + (index_barra2*3) + m - 3] = delta
            
        if index_barra2 != 0:
            #Derivada do fluxo de Potência ativa com relação ao ângulo na barra final
            for m in fases:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij)
                Bs = np.imag(Yij)
                Bsh = Bshmatrix[fase, m] /baseY
                tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra2*3) + m]
                ang_estimado2 = vet_estados[(index_barra2*3) + m]
                delta = tensao_estimada*tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2) - Bs*np.cos(ang_estimado-ang_estimado2))
                jacobiana[medida_atual][(index_barra2*3) + m - 3] = delta
        
        medida_atual += 1
        
    return medida_atual
      
def Derivadas_fluxo_pot_rat(jacobiana: np.array, fases: np.array, medida_atual: int, index_barra1: int, elemento: str,
                           barras: pd.DataFrame, nodes: dict, vet_estados: np.array, DSSCircuit, Ybus, baseva) -> int:
    barra1 = barras['nome_barra'][index_barra1]
    basekv = barras['Bases'][index_barra1]
    baseY = baseva / ((basekv*1000)**2)
    num_buses = DSSCircuit.NumBuses
    DSSCircuit.SetActiveElement(elemento)
    Bshmatrix = np.zeros((3, 3))
    if 'line' in elemento:
        Gsmatrix, Bsmatrix, Bshmatrix = Shuntmatrix(DSSCircuit, fases)
    barra2 = DSSCircuit.ActiveCktElement.BusNames[1]
    index_barra2 = barras[barras['nome_barra'] == barra2].index.values[0]
    
    for fase in fases:  
        no1 = nodes[barra1+f'.{fase+1}']

        tensao_estimada = vet_estados[(num_buses*3) + (index_barra1*3) + fase]
        ang_estimado = vet_estados[(index_barra1*3) + fase]
        
        #Derivada do fluxo de Potência reativa com relação a tensão na barra inicial
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2] / baseY
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m] / baseY
            if m == fase:
                delta = -tensao_estimada*(Bs+Bsh)
                for n in fases:
                    no2 = nodes[barra2+f'.{n+1}']
                    Yij = Ybus[no1, no2] / baseY
                    Gs = np.real(Yij)
                    Bs = np.imag(Yij)
                    Bsh = Bshmatrix[fase, n] / baseY
                    tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + n]
                    tensao_estimada3 = vet_estados[(num_buses*3) + (index_barra2*3) + n]
                    ang_estimado2 = vet_estados[(index_barra1*3) + n]
                    ang_estimado3 = vet_estados[(index_barra2*3) + n]
                    delta += tensao_estimada2*(Gs*np.sin(ang_estimado-ang_estimado2)-(Bs+Bsh)*np.cos(ang_estimado-ang_estimado2))
                    delta -= tensao_estimada3*(Gs*np.sin(ang_estimado-ang_estimado3)-Bs*np.cos(ang_estimado-ang_estimado3))

            else:
                ang_estimado2 = vet_estados[(index_barra1*3) + m]
                delta = tensao_estimada*(Gs*np.sin(ang_estimado-ang_estimado2)-(Bs+Bsh)*np.cos(ang_estimado-ang_estimado2))
                
            jacobiana[medida_atual][num_buses*3 + (index_barra1*3) + m - 3] = delta
            
        if index_barra1 != 0:
            #Derivada do fluxo de Potência reativa com relação ao ângulo na barra inicial
            for m in fases:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij)
                Bs = np.imag(Yij)
                Bsh = Bshmatrix[fase, m] / baseY
                if m == fase:
                    delta = -(tensao_estimada**2)*Gs
                    for n in fases:
                        no2 = nodes[barra2+f'.{n+1}']
                        Yij = Ybus[no1, no2] / baseY
                        Gs = np.real(Yij)
                        Bs = np.imag(Yij)
                        Bsh = Bshmatrix[fase, n] / baseY
                        tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + n]
                        tensao_estimada3 = vet_estados[(num_buses*3) + (index_barra2*3) + n]
                        ang_estimado2 = vet_estados[(index_barra1*3) + n]
                        ang_estimado3 = vet_estados[(index_barra2*3) + n]
                        delta += tensao_estimada*tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2)+(Bs+Bsh)*np.sin(ang_estimado-ang_estimado2))
                        delta -= tensao_estimada*tensao_estimada3*(Gs*np.cos(ang_estimado-ang_estimado3)+Bs*np.sin(ang_estimado-ang_estimado3))
                    
                else:
                    tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra1*3) + m]
                    ang_estimado2 = vet_estados[(index_barra1*3) + m]
                    delta = -tensao_estimada*tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2) + (Bs+Bsh)*np.sin(ang_estimado-ang_estimado2))
                
                jacobiana[medida_atual][(index_barra1*3) + m - 3] = delta
            
        #Derivada do fluxo de Potência reativa com relação a tensão na barra final
        for m in fases:
            no2 = nodes[barra2+f'.{m+1}']
            Yij = Ybus[no1, no2] / baseY
            Gs = np.real(Yij)
            Bs = np.imag(Yij)
            Bsh = Bshmatrix[fase, m] / baseY
            ang_estimado2 = vet_estados[(index_barra2*3) + m]
            delta = -tensao_estimada*(Gs*np.sin(ang_estimado-ang_estimado2)-Bs*np.cos(ang_estimado-ang_estimado2))
            jacobiana[medida_atual][num_buses*3 + (index_barra2*3) + m - 3] = delta
            
        if index_barra2 != 0:
            #Derivada do fluxo de Potência reativa com relação ao ângulo na barra final
            for m in fases:
                no2 = nodes[barra2+f'.{m+1}']
                Yij = Ybus[no1, no2] / baseY
                Gs = np.real(Yij)
                Bs = np.imag(Yij)
                Bsh = Bshmatrix[fase, m] / baseY
                tensao_estimada2 = vet_estados[(num_buses*3) + (index_barra2*3) + m]
                ang_estimado2 = vet_estados[(index_barra2*3) + m]
                delta = tensao_estimada*tensao_estimada2*(Gs*np.cos(ang_estimado-ang_estimado2) + Bs*np.sin(ang_estimado-ang_estimado2))
                jacobiana[medida_atual][(index_barra2*3) + m - 3] = delta
            
        medida_atual += 1
    
    return medida_atual
