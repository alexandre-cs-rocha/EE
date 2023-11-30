import numpy as np
import pandas as pd
import scipy as sp
import timeit
import csv

from pathlib import Path
from dss import DSS as dss_engine
from Jacobiana import Jacobiana
from Residuos import Residuo
from Filtragem import Pos_filtragem, Pre_filtragem

class EE():
    def __init__(self, master_path, Y_path, baseva: float = 10**6) -> None:
        self.DSSCircuit, self.DSSText, self.DSSObj, self.DSSMonitors = self.InitializeDSS()
        self.Y_path = Y_path
        self.baseva = baseva
        self.MasterFile = master_path

    def resolve_fluxo_carga(self):
        self.DSSText.Command = 'Clear'
        self.DSSText.Command = f'Compile {self.MasterFile}'

        self.iniciar_medidores()

        self.DSSText.Command = 'Solve'

    def InitializeDSS(self) -> tuple:
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

    def iniciar_medidores(self) -> None:
        for i, barra in enumerate(self.DSSCircuit.AllBusNames):
            self.DSSCircuit.SetActiveBus(barra)
            for j, elem in enumerate(self.DSSCircuit.Buses.AllPCEatBus):
                if 'Load' in elem or 'Generator' in elem or 'Vsource' in elem:
                    self.DSSText.Command = f'New Monitor.pqi{i}{j} element={elem}, terminal=1, mode=1, ppolar=no'
                    
            elem = self.DSSCircuit.Buses.AllPDEatBus[0]
            if elem != 'None':
                self.DSSCircuit.SetActiveElement(elem)
                if self.DSSCircuit.ActiveCktElement.BusNames[0].split('.')[0] == barra:
                    self.DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=1, mode=32'
                    
                elif self.DSSCircuit.ActiveCktElement.BusNames[1].split('.')[0] == barra:
                    self.DSSText.Command = f'New Monitor.v{i} element={elem}, terminal=2, mode=32'
                    
                else:
                    print('Deu errado')

    def organizar_nodes(self) -> dict:
        nodes = {}
        for i, node in enumerate(self.DSSCircuit.YNodeOrder):
            nodes[node.lower()] = i
        
        return nodes

    def indexar_barras(self) -> pd.DataFrame:
        #Designa indíces às barras
        nomes = []
        bases = []
        geracao = []
        for barra in self.DSSCircuit.AllBusNames:
            #if barra.isdigit(): è possível que o sourcebus e o reg não entrem para a EE
            self.DSSCircuit.SetActiveBus(barra)
            #Base é em fase-neutro
            base = self.DSSCircuit.Buses.kVBase
            nomes.append(barra)
            bases.append(base)
            geracao.append(self.DSSCircuit.Buses.AllPCEatBus[0] == 'Vsource.source')

        nomes = np.concatenate([nomes[1:], [nomes[0]]])
        bases = np.concatenate([bases[1:], [bases[0]]])
        geracao = np.concatenate([geracao[1:], [geracao[0]]])

        idx = [i for i in range(len(nomes))]
        inicial1 = [[0, 0, 0] for _ in range(len(nomes))]
        inicial2 = [[0, 0, 0] for _ in range(len(nomes))]
        
        barras = pd.DataFrame(columns=['nome_barra', 'Bases', 'Fases', 'Inj_pot_at', 'Inj_pot_rat', 'Flux_pot_at', 'Flux_pot_rat', 'Tensao', 'Inj_pot_at_est', 'Inj_pot_rat_est', 'Geracao'],
                            index=idx)
        
        barras['nome_barra'] = nomes
        barras.loc[idx, 'Bases'] = bases
        barras.loc[idx, 'Geracao'] = geracao
        
        for i in idx:
            barras['Inj_pot_at_est'][i] = inicial1[i]
            barras['Inj_pot_rat_est'][i] = inicial2[i]

        return barras

    def gera_medida_imperfeita(self, media: float) -> None:
        # Gerar fatores aleatórios com base na distribuição normal
        fatores = np.random.normal(media, self.dp, self.num_medidas)
        
        for i, medidas in enumerate(self.barras['Inj_pot_at']):
            self.barras['Inj_pot_at'][i] = medidas + medidas * fatores[i*3:(i+1)*3]

    def iniciar_vet_estados(self) -> np.array:
        vet_estados = np.zeros((len(self.barras)-1)*6)
        for i in range((len(self.barras)-1)*3, (len(self.barras)-1)*6):
            vet_estados[i] = 1

        for i in range((len(self.barras)-1)*3):
            if (i+1) % 3 == 0:
                vet_estados[i-1] = -120 * 2 * np.pi / 360
                vet_estados[i] = 120 * 2 * np.pi / 360
                
        return vet_estados
    
    def achar_index_barra(self, barras: pd.DataFrame, barra: int) -> int:
        #Retorna o index da barra do monitor ativo
        self.DSSCircuit.SetActiveElement(self.DSSMonitors.Element)
        
        self.DSSCircuit.SetActiveBus(self.DSSCircuit.ActiveCktElement.BusNames[barra])
        nome = self.DSSCircuit.Buses.Name
        
        return barras.index[barras['nome_barra'] == nome].to_list()[0]

    def pegar_fases(self) -> np.array:
        fases = self.DSSCircuit.ActiveCktElement.NodeOrder - 1
        fases = set(fases)
        fases.discard(-1)
        
        return fases

    def medidas(self, baseva: int) -> pd.DataFrame: 
        barras = self.indexar_barras()
        
        num_medidas = 0
        for idx in range(len(self.DSSCircuit.AllBusNames)):
            if not barras['Geracao'][idx]:
                barras['Inj_pot_at'][idx] = np.array([0, 0, 0])
                barras['Inj_pot_rat'][idx] = np.array([0, 0, 0])
                num_medidas += 6
        
        #Amostra e salva os valores dos medidores no sistema
        self.DSSMonitors.SampleAll()
        self.DSSMonitors.SaveAll()

        self.DSSMonitors.First
        for _ in range(self.DSSMonitors.Count):
            barra = self.DSSMonitors.Terminal - 1
            index_barra = self.achar_index_barra(barras, barra)
            fases = self.pegar_fases()
            barras['Fases'][index_barra] = fases
            matriz_medidas = self.DSSMonitors.AsMatrix()[0][2:]
            
            if 'pqij' in self.DSSMonitors.Name:
                if type(barras['Flux_pot_at'][index_barra]) != list and type(barras['Flux_pot_rat'][index_barra]) != list:
                    barras['Flux_pot_at'][index_barra] = []
                    barras['Flux_pot_rat'][index_barra] = []
                    
                elemento = self.DSSMonitors.Element
                self.DSSCircuit.ActiveCktElement.BusNames[1]
                medidas_at = np.full([3], np.NaN)
                medidas_rat = np.full([3], np.NaN)
                
                for i, fase in enumerate(fases):
                    medidas_at[fase] = matriz_medidas[i*2]*1000 / baseva
                    medidas_rat[fase] = matriz_medidas[i*2+1]*1000 / baseva
                    num_medidas += 2
                    
                barras['Flux_pot_at'][index_barra].append((elemento, medidas_at))
                barras['Flux_pot_rat'][index_barra].append((elemento, medidas_rat))
            
            elif 'pqi' in self.DSSMonitors.Name:
                medidas_at = np.zeros(3)
                medidas_rat = np.zeros(3)
                
                for i, fase in enumerate(fases):
                    medidas_at[fase] = matriz_medidas[i*2]
                    medidas_rat[fase] = matriz_medidas[i*2+1]
                    
                barras['Inj_pot_at'][index_barra] = -medidas_at*1000 / baseva
                barras['Inj_pot_rat'][index_barra] = -medidas_rat*1000 / baseva
                
            elif 'v' in self.DSSMonitors.Name:
                if type(barras['Tensao'][index_barra]) != np.ndarray:
                    medidas = np.zeros(3)

                    for i, fase in enumerate(fases):
                        medidas[fase] = matriz_medidas[i]
    
                    basekv = self.DSSCircuit.Buses.kVBase
                    barras['Tensao'][index_barra] = medidas / (basekv*1000)
                    if not barras['Geracao'][index_barra]:
                        num_medidas += 3
            
            self.DSSMonitors.Next
        return barras, num_medidas

    def Conserta_Ybus(self, Y_path):
        # Garantir que as matrizes primitivas sejam 3x3:
        with open(Y_path, newline='') as csvfile: # le o arquivo csv da ymatrix
            Y_array = list(csv.reader(csvfile))
        
        #Constroi a matrix Y
        Y_array = Y_array[1:]
        Y_matrix=[]
        barras_fases_existentes = []
        barras_fases_inexistentes = []
        for line in Y_array:
            aux_array = []
            line.remove('') # eliminar a ultima coluna
            aux_array.append(line[0])
            barras_fases_existentes.append(line[0])
            for elem in range(len(line[1:])//2):
                j, imag = (line[elem*2+2]).rsplit('+j ')
                aux_array.append(complex(line[elem*2+1])+complex(imag)*1j)
            Y_matrix.append(aux_array)
        
        barras_existentes=[]
        for i in barras_fases_existentes:
            x = i.rsplit('.')
            if x[0] in barras_existentes:
                continue
            else:
                barras_existentes.append(x[0])

        for barra in barras_existentes:
            for n in [1,2,3]:
                nome_barra = (f'{barra}'+'.'+f'{n}')
                nova_barra=[]
                if nome_barra in barras_fases_existentes: # verifica se existem as fases
                    continue
                else:
                    barras_fases_inexistentes.append(nome_barra)
                    (l,c) = np.shape(np.asarray(Y_matrix))
                    nova_barra.append(nome_barra)
                    nova_barra.extend(np.zeros(c-1))
                    idx=(barras_existentes.index(f'{barra}')*3 + n - 1)
                    Y_matrix.insert((idx), nova_barra) # insere uma linha de zeros
                    for var in Y_matrix:
                        var.insert(idx+1, 0) # insere uma coluna de zeros

        # Eimina o nome das barras da matriz:
        Y_auxiliar=[]
        o_Y = []
        for var in Y_matrix:
            Y_auxiliar.append(var[1:])
            x=(var[0]).lower()
            o_Y.append(x)
        # funcao para substituir nodes
        ordenacao_Y = {}
        for idx in range(len(o_Y)):
            aux=o_Y[idx]
            ordenacao_Y[aux] = idx

        Y_auxiliar=np.asarray(Y_auxiliar)
        #Ybus = sp.sparse.csc_matrix(self.DSSObj.YMatrix.GetCompressedYMatrix())
        Ybus=sp.sparse.csr_matrix(Y_auxiliar)


        # Conserta a parte associada a barra de geracao
        self.DSSCircuit.Transformers.First
        for _ in range(self.DSSCircuit.Transformers.Count):
            trafo = self.DSSCircuit.Transformers.Name
            self.DSSCircuit.SetActiveElement(trafo)
            num_phases = self.DSSCircuit.ActiveCktElement.NumPhases
            barras_conectadas = self.DSSCircuit.ActiveCktElement.BusNames
            self.DSSCircuit.SetActiveBus(barras_conectadas[0])
            basekv1 = self.DSSCircuit.Buses.kVBase
            self.DSSCircuit.SetActiveBus(barras_conectadas[1])
            basekv2 = self.DSSCircuit.Buses.kVBase
            if '.' in barras_conectadas[0] or '.' in barras_conectadas[1]:
                barras_conectadas[0] = barras_conectadas[0].split('.')[0]
                barras_conectadas[1] = barras_conectadas[1].split('.')[0]
                
            no1 = self.nodes[f"{barras_conectadas[0]}.{1}"]
            no2 = self.nodes[f"{barras_conectadas[1]}.{1}"]
            
            if basekv1 > basekv2:
                n = basekv1 / basekv2
                Ybus[no1:no1+num_phases, no2:no2+num_phases] = (Ybus[no1:no1+num_phases, no2:no2+num_phases])/n
                Ybus[no2:no2+num_phases, no1:no1+num_phases] = (Ybus[no2:no2+num_phases, no1:no1+num_phases])*n
            else:
                n = basekv2 / basekv1
                Ybus[no1:no1+num_phases, no2:no2+num_phases] = (Ybus[no1:no1+num_phases, no2:no2+num_phases])*n
                Ybus[no2:no2+num_phases, no1:no1+num_phases] = (Ybus[no2:no2+num_phases, no1:no1+num_phases])/n
                
            self.DSSCircuit.Transformers.Next

        self.DSSCircuit.Loads.First
        if self.DSSCircuit.Loads.IsDelta: #Caso carga em delta
            pass
        else: #Caso carga em estrela
            self.DSSMonitors.First
            for _ in range(self.DSSMonitors.Count):
                elemento = self.DSSMonitors.Element
                if 'load' in elemento:
                    self.DSSCircuit.SetActiveElement(elemento)
                    Yij = (self.DSSCircuit.Loads.kW - self.DSSCircuit.Loads.kvar*1j)*1000 / ((self.DSSCircuit.Loads.kV*1000)**2)
                    barra_correspondente = self.DSSCircuit.ActiveCktElement.BusNames[0]

                    for k in self.DSSCircuit.Buses.Nodes:
                        no1 = self.nodes[f"{barra_correspondente}.{k}"]
                        Ybus[no1, no1] -= Yij

                    self.DSSCircuit.Loads.Next

                if self.DSSMonitors.Next == None: #Critério de parada
                    break

        self.DSSCircuit.SetActiveElement('Vsource.source')
        Yprim = self.DSSCircuit.ActiveCktElement.Yprim
        real = Yprim[::2]
        imag = Yprim[1::2]*1j
        Yprim = real+imag
        Yprim = np.reshape(Yprim, (6, 6))
        Ybus[:3, :3] -= Yprim[:3, :3]

        return Ybus, ordenacao_Y

    def Calcula_pesos(self) -> tuple:
        inj_pot_at = np.vstack(self.barras['Inj_pot_at'].to_numpy()).flatten()
        inj_pot_rat = np.vstack(self.barras['Inj_pot_rat'].to_numpy()).flatten()
        tensao = np.vstack(self.barras['Tensao'].to_numpy()).flatten()
        
        medidas = np.concatenate([inj_pot_at[:-3], inj_pot_rat[:-3], tensao[:-3]])

        dp = (medidas * 0.01) / (3 * 100)
        dp[dp == 0] = 10**-5
        pesos = dp**-2
        pesos[pesos > 10**10] = 10**10
            
        matriz_pesos = np.diag(pesos)
        
        return matriz_pesos, np.abs(dp)
    
    def Calcula_residuo(self) -> np.array:
        residuo = Residuo()
        count = self.barras['Geracao'].value_counts()[1]
        
        angs = self.vet_estados[:(self.DSSCircuit.NumBuses-count)*3]
        tensoes = self.vet_estados[(self.DSSCircuit.NumBuses-count)*3:]
        ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
        tensoes_ref = self.barras['Tensao'][self.DSSCircuit.NumBuses-1]
        angs = np.concatenate((angs, ang_ref))
        tensoes = np.concatenate((tensoes, tensoes_ref))
        vet_estados_aux = np.concatenate((angs, tensoes))
        
        for idx, geracao in enumerate(self.barras['Geracao']):
            if not geracao:
                residuo.Residuo_inj_pot_at(vet_estados_aux, idx, self.barras, self.DSSCircuit.NumBuses, baseva, self.ordenacao_Y, self.Ybus)

                residuo.Residuo_inj_pot_rat(vet_estados_aux, idx, self.barras, self.DSSCircuit.NumBuses, baseva, self.ordenacao_Y, self.Ybus)
                
                residuo.Residuo_tensao(vet_estados_aux, idx, self.barras, self.DSSCircuit.NumBuses)
                
        for idx1, medidas in enumerate(self.barras['Flux_pot_at']):
            if type(medidas) == list:
                for medida in medidas:
                    elemento = medida[0]
                    fases = np.where((np.isnan(medida[1]) == False))[0]
                    residuo.Residuo_fluxo_pot_at(vet_estados, fases, idx1, elemento, 
                                        baseva, self.barras, self.DSSCircuit, self.ordenacao_Y, self.Ybus)
                
        for idx1, medidas in enumerate(self.barras['Flux_pot_rat']):
            if type(medidas) == list:
                for medida in medidas:
                    elemento = medida[0]
                    fases = np.where((np.isnan(medida[1]) == False))[0]
                    residuo.Residuo_fluxo_pot_rat(vet_estados, fases, idx1, elemento, 
                                        baseva, self.barras, self.DSSCircuit, self.ordenacao_Y, self.Ybus)
            
        return np.concatenate([residuo.vet_inj_at, residuo.vet_inj_rat, residuo.vet_tensao])

    def Calcula_Jacobiana(self) -> np.array:
        jacobiana = np.zeros((self.num_medidas, len(self.vet_estados)))
        
        jac = Jacobiana()

        count = self.barras['Geracao'].value_counts()[1]
        
        angs = self.vet_estados[:(self.DSSCircuit.NumBuses-count)*3]
        tensoes = self.vet_estados[(self.DSSCircuit.NumBuses-count)*3:]
        ang_ref = np.array([0, -2*np.pi/3, 2*np.pi/3])
        tensoes_ref = self.barras['Tensao'][self.DSSCircuit.NumBuses-1]
        angs = np.concatenate((angs, ang_ref))
        tensoes = np.concatenate((tensoes, tensoes_ref))
        vet_estados_aux = np.concatenate((angs, tensoes))
        
        medida_atual = 0
        for idx, medida in enumerate(self.barras['Inj_pot_at']):
            if type(medida) == np.ndarray and not self.barras['Geracao'][idx]:
                medida_atual = jac.Derivadas_inj_pot_at(jacobiana, medida_atual, idx, self.DSSCircuit.NumBuses, vet_estados_aux, self.barras,
                                                    self.ordenacao_Y, self.Ybus, baseva, count)
        
        for idx, medida in enumerate(self.barras['Inj_pot_rat']):
            if type(medida) == np.ndarray and not self.barras['Geracao'][idx]:
                medida_atual = jac.Derivadas_inj_pot_rat(jacobiana, medida_atual, idx, self.DSSCircuit.NumBuses, vet_estados_aux, self.barras,
                                                    self.ordenacao_Y, self.Ybus, baseva, count)
                
        for idx1, medida in enumerate(self.barras['Flux_pot_at']):
            if type(medida) == list:
                elemento = medida[0][0]
                fases = np.where((np.isnan(medida[0][1]) == False))[0]
                medida_atual = jac.Derivadas_fluxo_pot_at(jacobiana, fases, medida_atual, idx1, elemento, self.barras, self.ordenacao_Y, vet_estados_aux,
                                                    self.DSSCircuit, self.Ybus, baseva)
                
        for idx1, medida in enumerate(self.barras['Flux_pot_rat']):
            if type(medida) == list:
                elemento = medida[0][0]
                fases = np.where((np.isnan(medida[0][1]) == False))[0]
                medida_atual = jac.Derivadas_fluxo_pot_rat(jacobiana, fases, medida_atual, idx1, elemento, self.barras, self.ordenacao_Y, vet_estados_aux,
                                                    self.DSSCircuit, self.Ybus, baseva)
                
        for idx, medida in enumerate(self.barras['Tensao']):
            if type(medida) == np.ndarray and not self.barras['Geracao'][idx]:
                medida_atual = jac.Derivadas_tensao(jacobiana, self.barras, medida_atual, idx, self.DSSCircuit.NumBuses, count)
            
        return jacobiana

    def run(self, erro_max: float, lim_iter: int) -> np.array:
        self.resolve_fluxo_carga()
        
        self.barras, self.num_medidas = self.medidas(self.baseva)
        self.vet_estados = self.iniciar_vet_estados()
        self.nodes = self.organizar_nodes()
    
        #Ybus = sp.sparse.csc_matrix(self.DSSObj.YMatrix.GetCompressedYMatrix())
        self.Ybus, self.ordenacao_Y = self.Conserta_Ybus(self.Y_path)
        
        k = 0
        delx = 1
        while(np.max(np.abs(delx)) > erro_max and lim_iter > k):

            self.residuo = self.Calcula_residuo()

            self.jacobiana = self.Calcula_Jacobiana()
            
            self.matriz_pesos, self.dp = self.Calcula_pesos()
            
            self.gera_medida_imperfeita(0)
            
            #Calcula a matriz ganho
            matriz_ganho = self.jacobiana.T @ self.matriz_pesos @ self.jacobiana
            matrizganhopd= pd.DataFrame(matriz_ganho)
            pd.DataFrame.to_excel(matrizganhopd, r"C:\Users\souza\OneDrive\Documentos\Unb\pibic\EE-main-Class\matriz_ganho.xlsx")
            
            #Calcula o outro lado da Equação normal
            seinao = self.jacobiana.T @ self.matriz_pesos @ self.residuo

            delx = np.linalg.inv(matriz_ganho) @ seinao
            
            #Atualiza o vetor de estados
            self.vet_estados += delx
            
            k += 1

        return self.vet_estados

#Achar o path do script do OpenDSS
path = Path(__file__)
CurrentFolder = path.parent
Y_path = CurrentFolder / '13bus' / 'IEEE13Nodeckt_EXP_Y.CSV' #path da Ybus
MasterFile = CurrentFolder / '13bus' / 'IEEE13Nodeckt.dss'

baseva =  33.3 * 10**6

eesd = EE(MasterFile, Y_path, baseva)

'''time = timeit.timeit(lambda: eesd.run(10**-5, 10), number=1)
print(time)'''

vet_estados = eesd.run(10**-5, 10)

def Calcula_Gabarito(eesd):
    ang = np.array([])
    tensoes = np.array([])
    for barra in eesd.DSSCircuit.AllBusNames:
        eesd.DSSCircuit.SetActiveBus(barra)
        ang = np.concatenate([ang, eesd.DSSCircuit.Buses.puVmagAngle[1::2]*2*np.pi / 360])
        tensoes = np.concatenate([tensoes, eesd.DSSCircuit.Buses.puVmagAngle[::2]])

    gabarito = np.concatenate([ang[3:], tensoes[3:]])

    print(gabarito)

Calcula_Gabarito(eesd)
print(vet_estados)

pre_fitragem = Pre_filtragem()

medidas_tensao = pre_fitragem.vetoriza_medidas_tensao(eesd.barras)
medidas_tensao = np.insert(medidas_tensao, 10, 1.5)

print(pre_fitragem.mad_based_outlier(medidas_tensao, 3.5))

'''filtragem = Pos_filtragem(eesd.vet_estados, eesd.jacobiana, eesd.residuo, eesd.barras,
                          eesd.num_medidas, eesd.matriz_pesos, eesd.dp)

text = filtragem.teste_maior_residuo()

print(text)'''
