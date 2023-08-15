import numpy as np
import pandas as pd
import scipy as sp
from pathlib import Path
from dss import DSS as dss_engine
from Derivadas import *
from EESD import *

erro_max = 10**-3
lim_iter = 1

def vetor_medidas(barras: pd.DataFrame):
    vet_med = []
    dp = []
    for medidas in barras['Inj_pot_at']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    vet_med.append(0)

                    continue
                dp.append((medida * 0.01) / (3 * 100))
                vet_med.append((medida))
    
    for medidas in barras['Inj_pot_rat']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    vet_med.append(0)

                    continue
                dp.append((medida * 0.01) / (3 * 100))
                vet_med.append((medida))
                
    for medidas in barras['Flux_pot_at']:
        if type(medidas) == list:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    vet_med.append(0)

                    continue
                dp.append((medida * 0.01) / (3 * 100))
                vet_med.append((medida))
                
    for medidas in barras['Flux_pot_rat']:
        if type(medidas) == list:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    vet_med.append(0)

                    continue
                dp.append((medida * 0.01) / (3 * 100))
                vet_med.append((medida))
    
    for medidas in barras['Tensao']:
        if type(medidas) == np.ndarray:
            for medida in medidas:
                if medida == 0:
                    dp.append(1/10**5)
                    vet_med.append(0)

                    continue
                dp.append((medida * 0.002) / (3 * 100))
                basekv = DSSCircuit.Buses.kVBase
                vet_med.append((medida))

    vet_med = np.asarray(vet_med)
    dp = np.asarray(dp)
    print(len(dp))
    print(len(vet_med))

    return vet_med, dp

def teste_maior_residuo(vet_estados, barras, num_medidas, baseva, erro_max, lim_iter):
    # Calculo das matrizes de covariancia
    l=num_medidas
    Covar_medidas = np.zeros((l,l))
    vet_med, des_pad = vetor_medidas(barras)
    H = Calcula_Jacobiana(barras, vet_estados, num_medidas, baseva)
    W = Calcula_pesos(barras, num_medidas)    
    G = np.dot(np.dot(H.T, W), H)
    for i in range(l):
        Covar_medidas[i][i]= des_pad[i]**2
    print(np.diag(Covar_medidas))
    Covar_estados_estimados = np.linalg.inv(np.dot(np.dot(np.transpose(H),np.linalg.inv(Covar_medidas)),H))
    Covar_medidas_estimadas = np.dot(np.dot(H,Covar_estados_estimados),np.transpose(H))
    Covar_residuos = Covar_medidas-Covar_medidas_estimadas

    # Normalização das Covariâncias
    diag = np.diag(abs(Covar_residuos))

    # Matrix de covariancias normalizadas
    Rn = np.zeros((len(diag),len(diag)))
    for i in range(len(diag)):
        Rn[i][i] = float(diag[i])**(-1/2)

    #Covar_residuos_normalizados = np.dot(np.dot(Rn,Covar_residuos),Rn)
    Matriz_Sensibilidade=np.identity(l)-np.dot(np.dot(np.dot(H,np.linalg.inv(G)),H.T),W)

    # Vetor de covarancias normalizadas
    vetor_residuos = np.dot(Matriz_Sensibilidade,vet_med)

    #vetor_residuos_normalizados = np.dot(Rn,vetor_residuos)
    
    # Análise de erro e de b^
    Matriz_erros = vetor_residuos/np.diag(Matriz_Sensibilidade)
    Matriz_b = np.zeros((1,l))

    for i in range(len(des_pad)):
        Matriz_b[0][i] = abs(Matriz_erros[i]/des_pad[i])

    vetor_b = Matriz_b[0]
    maxb=np.max(vetor_b)
    index_maxb = np.argmax(vetor_b)
    med_b = vet_med[index_maxb]

    #Identificar a medida errada na função barras

    if abs(maxb)>3:
        vet_med_novo = vet_med
        vet_med_novo[index_maxb] = med_b-Matriz_erros[index_maxb]
        txt1 = f'A medida {med_b} provavelmente contém um erro grosseiro,'
        txt2 = f'uma estimativa para ela seria: {vet_med_novo[index_maxb]}.'
        #txt3 = f'Com essa nova medida, o resutado do estimador de estados seria: {EE(barras, vet_estados, matriz_pesos, baseva, erro_max, lim_iter)}'
        #txt4 = f'Nesse sentido, um novo teste de maior residuo teria como resultado: {teste_maior_residuo(EE(barras, vet_estados, matriz_pesos, W, 10**-3, 1),barras, num_medidas, baseva, erro_max, lim_iter)}'
        # {txt3}\n {txt4}
        print(f'{txt1}\n {txt2}\n ')
        return 'A estimação provavelmente contém erros grosseiros.'
    
    else:
        return 'A estimação provavelmente não contém erros grosseiros'

teste_maior_residuo(vet_estados, barras, num_medidas, baseva, erro_max, lim_iter)