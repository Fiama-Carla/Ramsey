import networkx as nx
import matplotlib.pyplot as plt
import random
from itertools import combinations, product
import numpy as np
import pickle
import os
import time
import csv


# Variável global para armazenar a biblioteca e matriz de alelos
biblioteca = None
alelos_dict = None

def carregar_dados(n,r):
    nome_arquivo= f'biblioteca_n={n}_r={r}.pkl'
    arquivo_alelo= f'alelo={n}_r={r}.pkl'
    """Carrega biblioteca e alelos_dict de arquivos, se existirem."""
    global biblioteca, alelos_dict

    # Verifica se os arquivos já existem para carregar dados
    if os.path.exists(nome_arquivo):
        with open(nome_arquivo, "rb") as file:
            biblioteca = pickle.load(file)
    
    if os.path.exists(arquivo_alelo):
        with open(arquivo_alelo, "rb") as file:
            alelos_dict = pickle.load(file)
    
    return biblioteca, alelos_dict

def salvar_dados(n,r):
    nome_arquivo= f'biblioteca_n={n}_r={r}.pkl'
    arquivo_alelo= f'alelo={n}_r={r}.pkl'
    with open(nome_arquivo, "wb") as file:
        pickle.dump(biblioteca, file)

    with open(arquivo_alelo, "wb") as file:
        pickle.dump(alelos_dict, file)

def inicializar_dados(n, r, k, s):
    """
    Inicializa a biblioteca e a matriz de alelos, carregando de arquivos 
    ou gerando novas estruturas, se necessário.
    """
    global biblioteca, alelos_dict
    biblioteca, alelos_dict = carregar_dados(n,r)
    
    if biblioteca is None or alelos_dict is None:
        print("Criando nova biblioteca e matriz de alelos...")
        biblioteca = gerar_biblioteca(n, r)
        alelos_dict = gerar_matriz_alelos(biblioteca, r,k,s)
        #aqui ja vai ser os alelos com os fitness
        salvar_dados(n,r)
    else:
        print("Dados carregados da memória!")

def gerar_biblioteca(n, r):
    global biblioteca
    if biblioteca is not None:
        # Retorna a biblioteca já carregada
        return biblioteca

    # Dicionário para armazenar todas as combinações de ligações nas células e na caspa
    biblioteca = {}
    indice = 0

    # Número de células completas e número de vértices na "caspa"
    num_celulas = n // r
    num_vertices_caspa = n % r

    # Configurações para as células de tamanho r
    combinacoes_celula = list(combinations(range(r), 2))
    for config in product([0, 1], repeat=len(combinacoes_celula)):
        biblioteca[indice] = {
            'tipo': 'celula',
            'arestas': [(combinacoes_celula[i], config[i]) for i in range(len(config))]
        }
        indice += 1

    # Configurações para a "caspa" com num_vertices_caspa vértices (se houver "caspa")
    if num_vertices_caspa > 1:  # Caspa com pelo menos 2 vértices para haver ligações
        combinacoes_caspa = list(combinations(range(num_vertices_caspa), 2))
        for config in product([0, 1], repeat=len(combinacoes_caspa)):
            biblioteca[indice] = {
                'tipo': 'caspa',
                'arestas': [(combinacoes_caspa[i], config[i]) for i in range(len(config))]
            }
            indice += 1
    elif num_vertices_caspa == 1:
        biblioteca[indice] = {
            'tipo': 'caspa',
            'arestas': []
        }

    return biblioteca

def gerar_alelos(celula_i, celula_j):
    #saber quantos vertices tem cada celula
    vertices_i = round((1+np.sqrt(1+8*len(celula_i['arestas'])))) // 2
    vertices_j = round((1+np.sqrt(1+8*len(celula_j['arestas'])))) // 2
    
    combinacoes_ligacoes = list(product(range(vertices_i), range(vertices_j)))
    alelos = list(product([0, 1], repeat=len(combinacoes_ligacoes)))
    return alelos

def indice_vetor(i, j, n):
    #Função para acessar os indices Aij da matriz no vetor binario
    if i > j:
      i, j = j, i
    return (i * (2 * n - i - 1) // 2) + (j - i - 1)

def criar_grafo_alelo(alelo, celula_i, celula_j):
    
    if len(celula_j['arestas']) < len(celula_i['arestas']):
      celula_i, celula_j = celula_j, celula_i
      
    
    #saber quantos vertices tem cada celula
    vertices_i = round((1+np.sqrt(1+8*len(celula_i['arestas'])))) // 2
    vertices_j = round((1+np.sqrt(1+8*len(celula_j['arestas'])))) // 2
                   
    n = vertices_i + vertices_j # total vertices
    arestas = len(celula_i['arestas']) + len(celula_j['arestas']) # Total de arestas no grafo
    

    # Inicializa o vetor binário (matriz de adjacência linearizada)
    grafo_disjunto = [0 for i in range ((n * (n - 1) // 2)) ] #Cada posição do vetor representa uma aresta

    # Adicionar arestas internas de celula_i
    for u, v in celula_i['arestas']:
        index = indice_vetor(u[0],u[1], n) #u é o par ordenado 
        grafo_disjunto[index] = v 

    # Adicionar arestas internas de celula_j
    for u, v in celula_j['arestas']:
        u_offset = u[0] + vertices_i  # Offset para os índices de celula_j
        v_offset = u[1] + vertices_i
        index = indice_vetor(u_offset, v_offset, n)
        grafo_disjunto[index] = v
  
    
    for i in range(vertices_i):  # Itera sobre todos os vértices de celula_i
      for j in range(vertices_j):  # Itera sobre todos os vértices de celula_j
        indexa = vertices_j * i + j  # Lineariza a combinação (i, j) no vetor alelo
        index = indice_vetor(i, j + vertices_i, n)  # Ajusta para índice global no grafo
        if indexa < len(alelo):  # Garantia de que indexa está no intervalo válid
            grafo_disjunto[index] = alelo[indexa]
                    
    return grafo_disjunto, n

def criar_grafo(n):
    """
    Essa função cria um grafo aleatorio {grafo_aleatorio} e passa ele pra um matriz adjacencia
    A matriz adjacencia {adj_matrix} é transformada em um vetor binario cujos elementos
    estão acima da diagonal principal pois a matriz é simetrica.
    Retorna o vetor binario {vetor_grafo} que representa o grafo e é o individuo da popilação inicial
    """
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

def roleta(populacao):
    
   #Função para realizar a seleção por roleta.
   # Verifica se há um indivíduo com fitness zero
    for individuo in populacao:
        if individuo[1] == 0: #posição referente ao fitness
            return [individuo]
        
    # Inverte os valores de fitness para que fitness menores sejam melhores
    #Extrai somente o valor da fitness
    inverso_fitness = [1.0 / individuo[1] for individuo in populacao]
    
    # Calcula a soma total dos fitness invertidos
    soma_total_inverso_fitness = sum(inverso_fitness)
    
    # Calcula a probabilidade de seleção para cada indivíduo
    probabilidades = [fit / soma_total_inverso_fitness for fit in inverso_fitness]
    #print("POPULATION", populacao)
    # Realiza a seleção por roleta
    selecionados = random.choices(populacao, weights=probabilidades, k=len(populacao))
    #print("seleciondoos", selecionados)
    return selecionados

def fitness(grafo, k, s, n):
    # Inicializa os índices dos vértices
    vertices = range(n)
    count_clique = 0
    count_nclique = 0
    
    # Contar quantos cliques de tamanho k existem no grafo
    for comb in combinations(vertices, k):
        clique = True
        for i in range(k):
            for j in range(i + 1, k):
                index = indice_vetor(comb[i], comb[j], n)
                if grafo[index] != 1:
                    clique = False
                    break
            if not clique:
                break
        if clique:
            count_clique += 1
            
    # Contar quantos cliques de tamanho s existem no grafo
    for comb in combinations(vertices, s):
        nclique = True
        for i in range(s):
            for j in range(i + 1, s):
                index = indice_vetor(comb[i], comb[j], n)
                if grafo[index] != 0:
                    nclique = False
                    break
            if not nclique:
                break
        if nclique:
            count_nclique += 1
    
    # A soma de cliques e ncliques (quanto menor, melhor)
    fitness_score = count_clique + count_nclique
    
    return fitness_score

# Gerar a matriz de alelos usando um dicionário
def gerar_matriz_alelos(biblioteca, r, k, s):
   
    global alelos_dict
    if alelos_dict is None:  # Inicializar o dicionário se estiver vazio
        alelos_dict = {}
        
    num_celulas = len(biblioteca)
    matriz_alelos = [[None for _ in range(num_celulas)] for _ in range(num_celulas)]
    
    alelo_index = 0
    
    # Preencher a matriz com todas as combinações de alelos entre as células, mas apenas para i != j
    for i in range(num_celulas):
        for j in range(num_celulas):
            # Acessar as células correspondentes aos índices i e j da biblioteca
            celula_i = biblioteca[i]
            celula_j = biblioteca[j]
                
            # Gerar os alelos entre as células i e j
            alelos = gerar_alelos(celula_i, celula_j)

            # Lista para armazenar os valores de fitness dos alelos
            fitness_list = []
            alelo_indices = []


            for alelo in alelos:
                # Criar um grafo com o alelo específico
                grafo_disjunto, num_vertices = criar_grafo_alelo(alelo, celula_i, celula_j)
                
                # Avaliar fitness do grafo
                fit = fitness(grafo_disjunto, k, s, num_vertices)
                fitness_list.append(fit)  # Adicionar o fitness à lista

                # Atribuir índice ao alelo
                alelo_indices.append(alelo_index)
                alelo_index += 1  # Incrementar o índice do alelo


            # Salvar os alelos e seus fitness no dicionário
            alelos_dict[(i, j)] = {
                'indices': alelo_indices,
                'alelos': alelos,
                'fitness': fitness_list
            }
            
    return alelos_dict                   
           


def extrair_alelo(grafo_disjunto, n, r):
   
    #Extrai as células e os alelos de um grafo representado por um vetor binário linearizado.
    
    # Inicializar as células
    num_celulas = n // r  # Número de células no grafo
    tamanho_caspa = n % r #colocar como var global
    celulas = []
    
  # Extrair arestas internas de cada célula
    for c in range(num_celulas):
        celula = {'id': [], 'vertices': list(range(c * r, (c + 1) * r)), 'arestas': []}
        for i in range(r):
            for j in range(i + 1, r):
                u = c * r + i
                v = c * r + j
                index = indice_vetor(u, v, n)
                peso = grafo_disjunto[index]  # 1 se houver aresta, 0 se não houver
                celula['arestas'].append(((u - c * r, v - c * r), peso))
        celulas.append(celula)

   # Extrair a caspa, se existir
    caspa = None
    if tamanho_caspa > 0:
        caspa = {'id': [], 'vertices': list(range(num_celulas * r, n)), 'arestas': []}
        for i in range(tamanho_caspa):
            for j in range(i + 1, tamanho_caspa):
                u = num_celulas * r + i
                v = num_celulas * r + j
                index = indice_vetor(u, v, n)
                peso = grafo_disjunto[index]
                caspa['arestas'].append(((i, j), peso))
        celulas.append(caspa)
        
    # Atribuir IDs às células com base na biblioteca
    for celula in celulas:
        celula['arestas'] = sorted(celula['arestas'])  # Ordenar as arestas para comparação
        for idx, data in biblioteca.items():
            # Ordenar as arestas da biblioteca antes de comparar
            if sorted(data['arestas']) == celula['arestas']:
                celula['id'].append(idx)  # Adicionar o índice correspondente
                break  # Encontrou o ID, não precisa continuar

    # Extrair os alelos (conexões entre células)
      
    alelos = []
    for i in range(num_celulas):
        for j in range(i + 1, num_celulas):
            alelo = []
            for u in range(r):
                for v in range(r):
                    index = indice_vetor(i * r + u, j * r + v, n)
                    alelo.append(grafo_disjunto[index])
            alelos.append(alelo)
            
    # Conexões entre células e a caspa
    if caspa:
        for i in range(num_celulas):
            alelo = []
            for u in range(r):
                for v in range(tamanho_caspa):
                    index = indice_vetor(i * r + u, num_celulas * r + v, n)
                    alelo.append(grafo_disjunto[index])
            alelos.append(alelo)
    #print("celulas extrair", celulas)
    #print("alelos extrair", alelos)
    
    return celulas, alelos

def encontrar_fitness_crossover(chave, alelo):
    alelo = tuple(alelo)
    chave = tuple(item[0] if isinstance(item, list) else item for item in chave)
    # Verifica se a chave existe no dicionário
    if chave in alelos_dict:
        valores = alelos_dict[chave]
        if alelo in valores['alelos']:
            # Obtém o índice do alelo
            indice = valores['alelos'].index(alelo)
            # Retorna o fitness correspondente
            return valores['fitness'][indice]
    return None  # Caso o alelo ou a chave não sejam encontrados


def criar_grafo_completo(celulas, alelos,n):
    """
    Cria um grafo combinando n células e n alelos.
    """
    # Determinar o número total de vértices

    total_vertices = n
    
    # Inicializar o vetor do grafo completo
    novo_grafo = [0 for _ in range((total_vertices * (total_vertices - 1)) // 2)]
    
    # Adicionar arestas internas de cada célula
    for celula in celulas:
        vertices = celula['vertices']
        for (u, v), valor in celula['arestas']:
            index = indice_vetor(vertices[u], vertices[v], total_vertices)
            novo_grafo[index] = valor
    
    # Adicionar arestas entre células usando os alelos
    alelo_index = 0
    for i in range(len(celulas) - 1):
        for j in range(i + 1, len(celulas)):
            if alelo_index >= len(alelos):
                raise ValueError("Número de alelos insuficiente para conectar todas as células.")
            
            # Obter os vértices e o alelo correspondente
            vertices_i = celulas[i]['vertices']
            vertices_j = celulas[j]['vertices']
            alelo = alelos[alelo_index]
            alelo_index += 1
            
            # Adicionar as arestas entre as células
            for local_i, u in enumerate(vertices_i):  # Índices locais e globais da célula i
                for local_j, v in enumerate(vertices_j):  # Índices locais e globais da célula j
                    alelo_pos = local_i * len(vertices_j) + local_j  # Índice no vetor alelo
                    if alelo_pos < len(alelo):  # Garantir que o índice está dentro do vetor alelo
                        valor_aresta = alelo[alelo_pos]
                        index = indice_vetor(u, v, total_vertices)
                        novo_grafo[index] = valor_aresta
    
    return novo_grafo

def combinar_grafos(grafo_pai, grafo_mama, n, r, p):
    """
    Combina dois grafos para gerar um novo grafo com base em uma probabilidade de escolha de células.
    """
    # Extrair as células e alelos dos dois grafos
    celulas_pai, alelo_pai = extrair_alelo(grafo_pai, n, r)
    celulas_mama, alelo_mama = extrair_alelo(grafo_mama, n, r)
    #print("celula papis", celulas_pai)
    #print("celula mamis", celulas_mama)
    #print("alelo papis", alelo_pai)
    #print("alelo mamis", alelo_mama)
  
    # Inicializar o novo grafo (vetor binário) com o mesmo tamanho
    novo_grafo = [0] * len(grafo_pai)
    
    # Passo 1: Escolher as células do novo grafo
    #TESTAR DE NOVO PRA DUAS CELULAS
    celulas_escolhidas = []
    for i in range(len(celulas_pai)):
        if random.random() < p:
            celulas_escolhidas.append(celulas_pai[i])
        else:
            celulas_escolhidas.append(celulas_mama[i])
        
    #print("Celulas Escolhida", celulas_escolhidas) 
    
    """
    #Passo 2: Escolher os alelos do novo grafo

    alelo_escolhido = []
    for i, (id1, id2) in enumerate(combinations([c['id'] for c in celulas_escolhidas], 2)):
        tipo1 = biblioteca[id1[0]]['tipo']
        tipo2 = biblioteca[id2[0]]['tipo']
        #print("tipooo", tipo1)
        #print("tipooo2", tipo2)
        
    """ 
    #Passo 2: Escolher os alelos do novo grafo
  
    alelo_escolhido = []
    #celula celula
    indice = 0
    for i in range(len(celulas_escolhidas) - 1):
        for j in range(i + 1, len(celulas_escolhidas)-1):
            chave_busca = [celulas_escolhidas[i]['id'], celulas_escolhidas[j]['id']]
            # Alelos do pai e da mãe
            alelo_pai_atual = alelo_pai[indice]
            alelo_mama_atual = alelo_mama[indice]
            #print("chave busca alelo combinar grafo1", chave_busca)
            # Fitness do pai e da mãe
            fit_pai = encontrar_fitness_crossover(chave_busca, alelo_pai_atual)
            #print("fit pai terminado")
            fit_mama = encontrar_fitness_crossover(chave_busca, alelo_mama_atual)
            #print("fit mae terminado")
            #print("Chave:", chave_busca)
            #print("Alelo Pai:", alelo_pai_atual, "Fitness Pai:", fit_pai)
            #print("Alelo Mãe:", alelo_mama_atual, "Fitness Mãe:", fit_mama)

            # Condicional para escolher o alelo com base no fitness
            if fit_pai < fit_mama:
                alelo_escolhido.append(alelo_pai_atual)
            else:
                alelo_escolhido.append(alelo_mama_atual)
            
            indice += 1

        #print("Alelo Escolhido:", alelo_escolhido)
    for i in range(len(celulas_escolhidas) - 1):
        j = len(celulas_escolhidas)-1
        chave_busca = [celulas_escolhidas[i]['id'], celulas_escolhidas[j]['id']]
        # Alelos do paie
        alelo_pai_atual = alelo_pai[indice]
        alelo_mama_atual = alelo_mama[indice]
        #print("chave busca alelo ombinar grafo2", chave_busca)
        # Fitness do pai e da mãe
        fit_pai = encontrar_fitness_crossover(chave_busca, alelo_pai_atual)
        #print("fit pai terminado")
        fit_mama = encontrar_fitness_crossover(chave_busca, alelo_mama_atual)
        #print("fit mae terminado")
        #print("Chave:", chave_busca)
        #print("Alelo Pai:", alelo_pai_atual, "Fitness Pai:", fit_pai)
        #print("Alelo Mãe:", alelo_mama_atual, "Fitness Mãe:", fit_mama)

        # Condicional para escolher o alelo com base no fitness
        if fit_pai < fit_mama:
            alelo_escolhido.append(alelo_pai_atual)
        else:
            alelo_escolhido.append(alelo_mama_atual)
        
        indice += 1

        #print("Alelo Escolhido:", alelo_escolhido)

    
    # Passo 3: Construir o novo grafo combinando células e alelos escolhidos
    
    novo_grafo = criar_grafo_completo(celulas_escolhidas,alelo_escolhido,n)
    #print("Grafito da mamae:",novo_grafo )
    
    return novo_grafo
    
    


def encontrar_solucao(tamanho_populacao, k, s, n, r, p, max_geracoes, tempo_limite):
    """
    Função principal para encontrar o contraexemplo usando algoritmo genético.
    Parâmetros:
    """
    
    populacao = gerar_populacao(tamanho_populacao,k, s, n)
    geracao = 0
    inicio = time.time()
    melhor_fitness_por_geracao = []
    pior_fitness_por_geracao = []
    
    while geracao < max_geracoes:
    #while geracao < max_geracoes and (time.time() - inicio) < tempo_limite:
        print("gerationes", geracao)    
        melhor_individuo = min(populacao, key=lambda x: x[1]) #veraqui
        #print("populatiopnes", populacao)
        #print("fitness", melhor_individuo[0])
        print(melhor_individuo[1])
        #print("population", populacao)
        pior_individuo = max(populacao, key=lambda x: x[1])#veraqui
        #print("pior", pior_individuo)
        melhor_fitness_por_geracao.append(melhor_individuo[1])#veraqui
        pior_fitness_por_geracao.append(pior_individuo[1]) #veraqui
        
        if melhor_individuo[1] == 0:
            total_tempo = time.time() - inicio
            print("melhor_fitness", melhor_individuo[1])
            #print("teste1")
            return melhor_individuo[0], melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao
        else:
            #print("teste two")
            populacao = evoluir(populacao, k, s, n,r,p)
            
        geracao += 1
            
    total_tempo = time.time() - inicio
    melhor_dos_piores = min(populacao, key=lambda x: x[1])
    # Se atingir o máximo de gerações sem solução, retorna o melhor encontrado
    #print("melhor_fitness por geração while", melhor_fitness)
   
    return melhor_dos_piores, melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao

    
def evoluir(populacao, k, s, n,r,p):
    """
    Função para evoluir a população usando seleção, crossover e mutação.
    """
    filhos = []
    
    #Selecionar os pais
    pais_selecionados = roleta(populacao)
    #print("pais selecionado",pais_selecionados)
    # Realizar crossover e mutação nos pais selecionados
    for i in range(0, len(pais_selecionados), 2):
        #print("i range", i)
        #pais_selecionados é uma lista do tipo pai[(vetor binario),(fitness)]
        #Para o crossover me interessa o indice [0] - referente ao vetor binario
        papis = pais_selecionados[i][0]
        #print("pai", papis)
        
        mama = pais_selecionados[i + 1][0]
        
       # print("mae", mama)
        #Crossover
        filho = combinar_grafos(papis, mama,n,r,p)

        #Mutacao
        filhom = mutacao(filho)
        #filho2 = mutacao(filho2, taxa_mutacao)

        #Adiciona na lista de filhos, os filhos gerados com suas respectivas fitness
        filhos.append((filhom, fitness(filhom, k, s, n)))

    # Substituir os menos aptos da população pelos filhos
    # Ordenar a população por fitness (ascendente: os piores no início)
    # Ordenar a população por fitness (decrescente: os piores no início)
    populacao.sort(key=lambda x: x[1], reverse=True)
    
    # Substituir os menos aptos pelos filhos
    for i in range(len(filhos)):
        populacao[i] = filhos[i]
    
    # Retornar a nova população
    #print("popualtiooooon", populacao)
    return populacao

def mutacao(filho):
    probabilidade = 0.7
    prob_mut = 0.2 if random.random() < probabilidade else 0.5
    
    for i in range(len(filho)):
        if random.random() < prob_mut:
            # Inverter o bit
            filho[i] = 1 - filho[i]
    return filho


def main():
    # ---------------------------
    # Parâmetros do problema
    # ---------------------------
    random.seed(time.time())
    parametros = [
        #(range (5,9),3,5), #resposta é 9
        #(range(12,18),4,4)#resaposta é 14
        (range (10,18),4,4) # resposta é 14
    ]
    
    p = 0.5
    r = 3
    tamanho_populacao = 100
    max_geracoes = 100
    tempo_limite_global = 3600
    tempo_limite_individual = 3600
    x = 1

    # Abrindo o CSV uma vez no início
    with open('resultados_geneticotesteteste2.csv', mode='a', newline='') as arquivo:
        escritor_csv = csv.writer(arquivo, delimiter=';')
        escritor_csv.writerow(['n', 'k', 's', 'Media_Tempo', 'Media_Geracoes', 
                               'Melhor_Fitness', 'Pior_Fitness', 'Desvio_Padrao_Geracoes', 
                               'Variancia_Geracoes', 'Desvio_Tempo', 'Variancia_Tempo'])

        for parametro in parametros:
            if isinstance(parametro[0], range):
                n_values = parametro[0]
            else:
                n_values = [parametro[0]]

            k = parametro[1]
            s = parametro[2]

            for n in n_values:
                inicializar_dados(n, r, k, s)

                # Variáveis para armazenar métricas por conjunto de execuções
                tempos_execucao = []
                geracoes_por_execucao = []
                melhor_fitness = []
                pior_fitness = []

                tempo_inicio_teste = time.time()
                print(f"\n--- Resultados para n = {n}, k = {k}, s = {s} ---")
                for execucao in range(x):
                    if (time.time() - tempo_inicio_teste) >= tempo_limite_global:
                        print("Parou devido ao tempo limite global.")
                        break

                    inicio_execucao = time.time()

                    melhor_individuo, melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao = encontrar_solucao(
                        tamanho_populacao, k, s, n, r, p, max_geracoes, tempo_limite_individual
                    )

                    tempos_execucao.append(total_tempo)
                    geracoes_por_execucao.append(geracao)
                    melhor_fitness.append(melhor_fitness_por_geracao[-1])
                    pior_fitness.append(pior_fitness_por_geracao[-1])

                    if total_tempo >= tempo_limite_individual:
                        print(f"Execução {execucao + 1} parou devido ao tempo limite individual.")
                        break

                # ---------------------------
                # Calculando as médias e estatísticas
                # ---------------------------
                if geracoes_por_execucao:
                    media_geracoes = np.mean(geracoes_por_execucao)
                    media_tempo = np.mean(tempos_execucao)
                    melhor_fitness_global = min(melhor_fitness)
                    pior_fitness_global = max(pior_fitness)

                    desvio_padrao_geracoes = np.std(geracoes_por_execucao)
                    variancia_geracoes = np.var(geracoes_por_execucao)

                    desvio_tempo = np.std(tempos_execucao)
                    variancia_tempo = np.var(tempos_execucao)

                    # Escrevendo os resultados imediatamente no CSV
                    escritor_csv.writerow([
                        n, k, s,
                        f"{media_tempo:.4f}",
                        int(media_geracoes),
                        melhor_fitness_global,
                        pior_fitness_global,
                        int(desvio_padrao_geracoes),
                        int(variancia_geracoes),
                        f"{desvio_tempo:.4f}",
                        f"{variancia_tempo:.4f}"
                    ])

                    # Exibindo resultados no terminal
                    print(f"\n--- Resultados para n = {n}, k = {k}, s = {s} ---")
                    print(f"Média de gerações: {media_geracoes}")
                    print(f"Desvio padrão das gerações: {desvio_padrao_geracoes}")
                    print(f"Variância das gerações: {variancia_geracoes}")
                    print(f"Tempo médio de execução: {media_tempo:.4f} segundos")
                    print(f"Desvio padrão do tempo de execução: {desvio_tempo:.4f} segundos")
                    print(f"Variância do tempo de execução: {variancia_tempo:.4f} segundos")
                    print(f"Melhor fitness encontrado: {melhor_fitness_global}")
                    print(f"Pior fitness encontrado: {pior_fitness_global}")

if __name__ == "__main__":
    main()
"""
def main():
    # ---------------------------
    # Parâmetros do problema
    # ---------------------------
    
    parametros = [
        #(range (9,15),4,4)
        (11,4,4)
        #(range (5,14),4,4), #resposta é 9
        #(range (5,15),4,4), #resaposta é 18
        #(range (5,13),3,5) # resposta é 14  
    ]
   
    # ---------------------------
    # Parâmetros do teste
    # ---------------------------
    
    
    p = 0.5 #Probabilidade crossover
    r = 3   # Tamanho da célula
    tamanho_populacao = 100 # Tamanho da população
    max_geracoes = 100  # Critério de parada do algoritmo genético
    tempo_limite_global = 25*3600  #Tempo limite de 1 dia por cada R(k,s) (em segundos)
    tempo_limite_individual = 5*3600  # Tempo limite de cada execução 5 horas
    x = 10  # Número de execuções por teste
    melhor_fitness = []  # Lista para armazenar os melhores fitness por execução
    pior_fitness = []  # Lista para armazenar os piores fitness por execução
    geracoes_por_execucao = []
    tempos_execucao = []
    # ---------------------------
    # Salvar os dados
    # ---------------------------
    with open('resultados_genetico2.csv', mode='w', newline='') as arquivo:
        escritor_csv = csv.writer(arquivo,delimiter=';')
        escritor_csv.writerow(['n', 'k', 's', 'Media_Tempo', 'Media_Geracoes', 'Melhor_Fitness', 'Pior_Fitness','desvio_padrao_geracoes', 'variancia_geracoes', 'desvio_tempo', 'variancia_tempo'])

        for parametro in parametros:
            if isinstance(parametro[0], range):
                n_values = parametro[0]
            else:
                n_values = [parametro[0]]

            k = parametro[1]
            s = parametro[2]
            

            for n in n_values:
                # ---------------------------
                # Criar biblioteca e alelos
                # ---------------------------

                inicializar_dados(n, r, k, s)
               
                #Exibir configurações de biblioteca e alelos para verificação
              
                print("Biblioteca gerada:")
                for b, i in biblioteca.items():
                    print(b, i)
               
                #exibir_alelos()
                
                # ---------------------------
                # Variaveis para salvar as metricas
                # ---------------------------
                
                tempo_inicio_teste = time.time()  # Inicia contagem do tempo do teste
                
                for execucao in range(x):
                    print(f"execução {x}, n = {n}, {k} {s}")
                    if (time.time() - tempo_inicio_teste) >= tempo_limite_global:
                        print("Parou devido ao tempo limite global.")
                        break

                    inicio_execucao = time.time()

                    melhor_individuo, melhor_fitness_por_geracao, pior_fitness_por_geracao, total_tempo, geracao = encontrar_solucao(tamanho_populacao, k, s, n, r, p, max_geracoes, tempo_limite_individual)
                    
                    tempos_execucao.append(total_tempo)
                    geracoes_por_execucao.append(geracao)
                    melhor_fitness.append(melhor_fitness_por_geracao[-1]) #veraqui
                    pior_fitness.append(pior_fitness_por_geracao[-1])

                    if total_tempo >= tempo_limite_individual:
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
                escritor_csv.writerow([n, k, s, f"{media_tempo:.4f}", int(media_geracoes), melhor_fitness_global, pior_fitness_global,
                                       int(desvio_padrao_geracoes), int(variancia_geracoes), f"{desvio_tempo:.4f}", f"{variancia_tempo:.4f}"])
                
               # Exibe os resultados
                
                print(f"\n--- Resultados ---")
                print(f"Vertices {n}, k {k} e s {s}")
                print(f"Média de geracoes: {media_geracoes}")
                print(f"Desvio padrão das geracoes: {desvio_padrao_geracoes}")
                print(f"Variância das geracoes: {variancia_geracoes}")
                print(f"Tempo médio de execução: {media_tempo:.4f} segundos")
                print(f"Desvio padrão do tempo de execução: {desvio_tempo:.4f} segundos")
                print(f"Variância do tempo de execução: {variancia_tempo:.4f} segundos")
                print(f"Melhor fitnes  encontrado {melhor_fitness_global} ")
                print(f" Pior fitness encontrado  {pior_fitness_global}")
                print(f"Melhor grafo encontrado  {melhor_individuo}")
                

if __name__ == "__main__":
    main()

"""



    
    
