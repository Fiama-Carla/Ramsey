import networkx as nx
import matplotlib.pyplot as plt
import random
import statistics
from itertools import combinations
import numpy as np
import time
import csv

def criar_grafo(n):
    
    #cria grafo aleatorio 
    grafo_aleatorio = nx.erdos_renyi_graph(n, 0.5)

    #Matriz adjacencia do grafo
    adj_matrix = nx.to_numpy_array(grafo_aleatorio, dtype=int)

    #Elementos acima da diagonal principal pois a matriz é simetrica
    vetor_grafo = adj_matrix[np.triu_indices(len(adj_matrix), k=1)]
    return vetor_grafo

def gerar_populacao(tamanho_populacao, k, s, n):
#Função para gerar uma população inicial de grafos aleatórios.
    populacao = []
    for _ in range(tamanho_populacao):
        individuo = criar_grafo(n)
        populacao.append((individuo, fitness(individuo, k, s, n))) 
    return populacao

def indice_vetor(i, j, n):
    #Função para acessar os indices Aij da matriz no vetor binario
    if i > j:
      i, j = j, i
    return (i * (2 * n - i - 1) // 2) + (j - i - 1)

def fitness(vetor_binario, k, s, n):
    
    #Inicializa os indices
    vertices = range(n)
    count_clique = 0
    count_nclique = 0
    
    #Contar quantas cliques de tamanho k tem no grafo
    for comb in combinations(vertices, k):
        clique = True
        for i in range(k):
            for j in range(i + 1, k):
                index = indice_vetor(comb[i], comb[j],n)
                if vetor_binario[index] != 1:
                    clique = False
                    break
            if not clique:
                break
        if clique:
            count_clique += 1
            
#Contar quantas ncliques de tamanho s tem no grafo
    for comb in combinations(vertices, s):
        nclique = True
        for i in range(s):
            for j in range(i + 1, s):
                index = indice_vetor(comb[i], comb[j],n)
                if vetor_binario[index] != 0:
                    nclique = False
                    break
            if not nclique:
                break
        if nclique:
            count_nclique += 1
    
    # Soma dos valores
    fitness_score = count_clique + count_nclique 
    
    return fitness_score

def roleta(populacao):
   #Função para realizar a seleção por roleta.
    
   # Verifica se há um indivíduo com fitness zero
    for individuo in populacao:
        if individuo[1] == 0:
            return [individuo]
        
    # Inverte os valores de fitness para que fitness menores sejam melhores
    inverso_fitness = [1.0 / individuo[1] for individuo in populacao]
    
    # Calcula a soma total dos fitness invertidos
    soma_total_inverso_fitness = sum(inverso_fitness)
    
    # Calcula a probabilidade de seleção para cada indivíduo
    probabilidades = [fit / soma_total_inverso_fitness for fit in inverso_fitness]
    
    # Realiza a seleção por roleta
    selecionados = random.choices(populacao, weights=probabilidades, k=len(populacao))
    return selecionados

def crossover(pai1, pai2):
   
    # Garantir que os pais têm o mesmo comprimento
    assert len(pai1) == len(pai2), "Os pais devem ter o mesmo comprimento"
    
    # Determinar o ponto de corte (metade do vetor)
    ponto = len(pai1) // 2
    
    # Criar os filhos
    filho1 = pai1[:ponto].copy()
    filho1 = np.append(filho1, pai2[ponto:])
    
    filho2 = pai2[:ponto].copy()
    filho2 = np.append(filho2, pai1[ponto:])
    
    return filho1, filho2

def mutacao(cromossomo, prob2):
    
    probabilidade_mutacao = 0.1 if random.random() < prob2 else 0.5
    
    for i in range(len(cromossomo)):
        if random.random() < probabilidade_mutacao:
            # Inverter o bit
            cromossomo[i] = 1 - cromossomo[i]
    return cromossomo

def evoluir(populacao, taxa_mutacao, prob2, k, s, n, percentual_ruido):
    """
    Função para evoluir a população usando seleção, crossover e mutação.
    Também adiciona ruído na população conforme o percentual de ruído especificado.
    """
    filhos = []
    
    # Selecionar os pais
    pais_selecionados = roleta(populacao)

    # Se roleta encontrar um indivíduo com fitness zero, retorna imediatamente
    if len(pais_selecionados) == 1 and pais_selecionados[0][1] == 0:
        return [pais_selecionados[0]]
    
    # Realizar crossover e mutação nos pais selecionados
    for i in range(0, len(pais_selecionados), 2):
        pai1 = pais_selecionados[i][0]
        pai2 = pais_selecionados[i + 1][0]

        # Crossover
        filho1, filho2 = crossover(pai1, pai2)

        # Mutação
        if random.random() < taxa_mutacao:
            filho1 = mutacao(filho1, prob2)
        if random.random() < taxa_mutacao:
            filho2 = mutacao(filho2, prob2)

        # Adiciona os filhos à lista
        filhos.append((filho1, fitness(filho1, k, s, n)))
        filhos.append((filho2, fitness(filho2, k, s, n)))
    
    # Adicionar ruído à população (indivíduos aleatórios)
    tamanho_populacao = len(populacao)
    qtd_ruido = int(percentual_ruido * tamanho_populacao)  # Percentual de ruído da população
    
    # Gerar indivíduos aleatórios
    individuos_aleatorios = []
    for _ in range(qtd_ruido):
        individuo_aleatorio = criar_grafo(n)
        individuos_aleatorios.append((individuo_aleatorio, fitness(individuo_aleatorio, k, s, n)))

    # Reduzir a população evoluída para garantir espaço para o ruído
    filhos = filhos[:tamanho_populacao - qtd_ruido]
    
    # Adicionar indivíduos aleatórios
    filhos.extend(individuos_aleatorios)

    return filhos


def encontrar_solucao_com_ruido(tamanho_populacao, n, k, s, taxa_mutacao, prob2, max_geracoes, tempo_limite, ruido_percentual):
    """
    Função para encontrar a solução evoluindo a população com incremento de ruído de 10% até 100%.
    """
    populacao = gerar_populacao(tamanho_populacao, k, s, n)
    geracao = 0
    inicio = time.time()
    melhor_fitness_por_geracao = []
    pior_fitness_por_geracao = []
   

    while geracao < max_geracoes and time.time() - inicio < tempo_limite:
        print("gerationes", geracao)
        #while geracao < max_geracoes and (time.time() - inicio) < tempo_limite:
        populacao = evoluir(populacao, taxa_mutacao, prob2, k, s, n, ruido_percentual)

        melhor_individuo = min(populacao, key=lambda x: x[1])
        print("melior fitnesses", melhor_individuo[1])
        pior_individuo = max(populacao, key=lambda x: x[1])

        melhor_fitness_por_geracao.append(melhor_individuo[1])
        pior_fitness_por_geracao.append(pior_individuo[1])
        

        # Verifica se há um indivíduo com fitness zero
        if melhor_individuo[1] == 0:
            total_tempo = time.time() - inicio
            
            return melhor_individuo[0], melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao
        
        geracao += 1
        """
        if time.time() - inicio > tempo_limite:
            total_time = time.time() - inicio
            break
        """
    total_tempo = time.time() - inicio
    return melhor_individuo[0], melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao


def gerar_populacao(tamanho_populacao, k, s, n):
    populacao = [(criar_grafo(n), fitness(criar_grafo(n), k, s, n)) for _ in range(tamanho_populacao)]
    return populacao


def main():
    # ---------------------------
    # Parâmetros do problema
    # ---------------------------
    
        
    
    parametros = [
        #(10,3, 5) #35 resposta é 14 - testei ate 9
        #(10, 4, 4) #resposta é 18 testeiu ate 9
        #(17, 4,5),
        #(range (16,17),4,4)
        #(range (5,9),3,5), #resposta é 9
        (13,4,4) #resaposta é 14
        #(range (5,9),3,4) # resposta é 14
    ]
    
   
    # ---------------------------
    # Parâmetros do teste
    # ---------------------------
    taxa_mutacao = 0.20  # taxa para mutação do vetor binário
    prob2 = 0.7  # taxa para mutação da aresta
    tamanho_populacao = 100 # Tamanho da população
    max_geracoes = 100  # Critério de parada do algoritmo genético
    tempo_limite_global = 3600  #Tempo limite de 1 dia por cada R(k,s) (em segundos)
    tempo_limite_individual = 1800  # Tempo limite de cada execução 5 horas
    x = 1  # Número de execuções por teste
    percentual_ruido = 0.5
    
    melhor_fitness = []  # Lista para armazenar os melhores fitness por execução
    pior_fitness = []  # Lista para armazenar os piores fitness por execução
    geracoes_por_execucao = []
    tempos_execucao = []
    # ---------------------------
    # Salvar os dados
    # ---------------------------
    with open('resultados_geneticoteste.csv', mode='a', newline='') as arquivo:
        escritor_csv = csv.writer(arquivo,delimiter=';')
        escritor_csv.writerow(['n', 'k', 's', 'Media_Tempo', 'Media_Geracoes', 'Melhor_Fitness', 'Pior_Fitness', 'ruido',
                               'desvio_padrao_geracoes', 'variancia_geracoes', 'desvio_tempo', 'variancia_tempo'])

        for parametro in parametros:
            if isinstance(parametro[0], range):
                n_values = parametro[0]
            else:
                n_values = [parametro[0]]

            k = parametro[1]
            s = parametro[2]
            
            for n in n_values:
                
                tempo_inicio_teste = time.time()  # Inicia contagem do tempo do teste
                
                #for percentual_ruido in np.arange(0.1, 1.1, 0.1):  # Aumentando o ruído de 10% até 100%
                for execucao in range(x):
                    print(f"Execução {execucao} {n} {k} {s}.")
                    if (time.time() - tempo_inicio_teste) >= tempo_limite_global:
                        print("Parou devido ao tempo limite global.")
                        break

                    inicio_execucao = time.time()

                    melhor_individuo, melhor_fitness_por_geracao, pior_fitness_por_geracao, tempo_execucao, geracao = encontrar_solucao_com_ruido(
                        tamanho_populacao, n, k, s, taxa_mutacao, prob2, max_geracoes, tempo_limite_individual, percentual_ruido)

                    #fim_execucao = time.time()
                    #tempo_execucao = fim_execucao - inicio_execucao

                    tempos_execucao.append(tempo_execucao)
                    geracoes_por_execucao.append(geracao)
                    melhor_fitness.append(melhor_fitness_por_geracao[-1])
                    pior_fitness.append(pior_fitness_por_geracao[-1])   

                    if tempo_execucao >= tempo_limite_individual:
                        print(f"Iteração {execucao + 1} parou devido ao tempo limite individual.")
                        break

                # ---------------------------
                # Calculando as médias
                # ---------------------------

                media_geracoes = np.mean(geracoes_por_execucao)
                media_tempo = np.mean(tempos_execucao)
                melhor_fitness_global = min(melhor_fitness)
                pior_fitness_global = max(pior_fitness)
                
                desvio_padrao_geracoes = np.std(geracoes_por_execucao)
                variancia_geracoes = np.var(geracoes_por_execucao) 
                
                desvio_tempo = np.std(tempos_execucao) 
                variancia_tempo = np.var(tempos_execucao) 

                # Salvando resultados no CSV
                escritor_csv.writerow([n, k, s, f"{media_tempo:.4f}", int(media_geracoes), melhor_fitness_global, pior_fitness_global, percentual_ruido,
                                       int(desvio_padrao_geracoes), int(variancia_geracoes), f"{desvio_tempo:.4f}", f"{variancia_tempo:.4f}"])
                
               # Exibe os resultados
                print(f"\n--- Resultados {geracao} ---")
                print(f"Vertices {n}, k {k} e s {s}")
                print(f"Média de geracoes: {media_geracoes}")
                print(f"Desvio padrão das geracoes: {desvio_padrao_geracoes}")
                print(f"Variância das geracoes: {variancia_geracoes}")
                print(f"Tempo médio de execução: {media_tempo:.4f} segundos")
                print(f"Desvio padrão do tempo de execução: {desvio_tempo:.4f} segundos")
                print(f"Variância do tempo de execução: {variancia_tempo:.4f} segundos")
                print(f"Grafo melhor  encontrado {melhor_individuo} ")
                print(f"Melhor fitnes  encontrado {melhor_fitness_global} ")
                print(f"Pior fitness encontrado  {pior_fitness_global}")

if __name__ == "__main__":
    main()
                
