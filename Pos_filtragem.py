import math
import numpy as np

class Pos_filtragem():
    def __init__(self, vet_estados, jacobiana, residuo, barras, num_medidas, matriz_pesos, des_padrao) -> None:
        self.vet_estados = vet_estados
        self.jacobiana = jacobiana
        self.residuo = residuo
        self.barras = barras
        self.num_medidas = num_medidas
        self.matriz_pesos = matriz_pesos
        self.des_padrao = des_padrao

    def teste_maior_residuo(self):
        # Calculo das matrizes de covariancia
        l=self.num_medidas
        Covar_medidas = np.zeros((l,l))
        vet1, vet2, vet3 = [],[],[]
        for barra in self.barras['Inj_pot_at']:
            for medida in barra:
                vet1.append(medida)
        
        for barra in self.barras['Inj_pot_rat']:
            for medida in barra:
                vet2.append(medida)
        
        for barra in self.barras['Tensao']:
            for medida in barra:
                vet3.append(medida)

        vet_med = np.asarray(vet1+vet2+vet3)
        H = self.jacobiana
        W = self.matriz_pesos   
        #Calcula a matriz ganho
        matriz_ganho = H.T @ W @ H
        G = matriz_ganho

        for i in range(l):
            Covar_medidas[i][i]= self.des_padrao[i]**2
        Covar_estados_estimados = np.linalg.inv(G)
        Covar_medidas_estimadas = H @ Covar_estados_estimados @ H.T
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
        #vetor_residuos = np.dot(Matriz_Sensibilidade,vet_med)

        #vetor_residuos_normalizados = np.dot(Rn,vetor_residuos)
        
        # Análise de erro e de b^
        Matriz_erros = self.residuo/np.diag(Matriz_Sensibilidade)
        Matriz_b = np.zeros((1,l))

        for i in range(len(self.des_padrao)):
            Matriz_b[0][i] = abs(Matriz_erros[i]/self.des_padrao[i])

        vetor_b = Matriz_b[0]
        maxb=np.max(vetor_b)
        index_maxb = np.argmax(vetor_b)
        med_b = vet_med[index_maxb]

        #Identificar a medida errada na função barras

        if abs(maxb)>3:
            vet_med_novo = vet_med
            vet_med_novo[index_maxb] = med_b-Matriz_erros[index_maxb]
            if math.isclose(med_b, vet_med_novo[index_maxb]) is True or med_b == 0:
                return 'A estimação provavelmente não contém erros grosseiros'
            else:
                txt1 = f'A medida  de índice {index_maxb} e valor {med_b} provavelmente contém um erro grosseiro, pois apresenta um b igual a {maxb}'
                txt2 = f'uma estimativa para ela seria: {vet_med_novo[index_maxb]}.'
                return f'A estimação provavelmente contém erros grosseiros.{txt1} \n {txt2}'
        
        else:
            return 'A estimação provavelmente não contém erros grosseiros'