import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import csv
import random

# Dicionario dos números de Ramsey conhecidos ou aproximados
# Para valores aproximados, usaremos o limite maximo
ramsey = {
    (1, 1): 1,
    (1, 3): 1,
    (3, 1): 1,
    (2, 1): 1,
    (1, 2): 1,
    (2, 2): 2,
    (2, 3): 3,
    (3, 2): 3,
    (2, 4): 4,
    (4, 2): 4,
    (2, 5): 5,
    (5, 2): 5,
    (2, 6): 6,
    (6, 2): 6,
    (2, 7): 7,
    (7, 2): 7,
    (2, 8): 8,
    (8, 2): 8,
    (3, 3): 6,
    (4, 4): 18,
    (3, 4): 9,
    (4, 3): 9,
    (5, 5): 48,
    (6, 6): 165,
    (7, 7): 540,
    (8, 8): 1870,
    (3, 5): 14,
    (5, 3): 14,
    (4, 5): 25,
    (5, 4): 25,
    (6, 3): 18,
    (3, 6): 18,
    (6, 4): 41,
    (4, 6): 41,
    (6, 5): 87,
    (5, 6): 87,
    (7, 3): 23,
    (3, 7): 23,
    (7, 4): 61,
    (4, 7): 61,
    (7, 5): 143,
    (5, 7): 143,
    (7, 6): 298,
    (6, 7): 298,
    (8, 3): 28,
    (3, 8): 28,
    (8, 4): 84,
    (4, 8): 84,
    (8, 5): 216,
    (5, 8): 216,
    (8, 6): 495,
    (6, 8): 495,
    (8, 7): 1031,
    (7, 8): 1031,   
}

def turan(k, n):

    """
    Calcula o número máximo de arestas em um grafo com n vértices
    que não contém um subgrafo completo K_k usando a fórmula de Turán.
    :return: Número máximo de arestas no grafo sem K_k
    """   

    if k >= 2:
        s= n%k
        turan = ((1-1/k)*(n**2 - s**2))//2 + (s*(s-1))//2
        
    else:
        turan = 0
    return turan


def verifica(grafo, k, s):

    #Variaveis

    n = len(grafo.nodes)

    num_arestas_G = grafo.number_of_edges()

    grafo_complementar = nx.complement(grafo)

    num_arestas_G_complementar = grafo_complementar.number_of_edges()

 

    #Verifica condição de Turan
    if num_arestas_G > turan(k-1, n):
        return False

    if num_arestas_G_complementar > turan(s-1, n):
        return False

    # Verifica se cada vértice participa de um clique

    for v in list(grafo.nodes):

        if grafo.degree[v] == n - 1:

            novo_grafo = grafo.subgraph(set(grafo.nodes) - {v})

            return verifica(novo_grafo, k - 1, s)

    # Verifica se existe algum vértice com grau 0

    for v in list(grafo.nodes):

        if grafo.degree[v] == 0:

            novo_grafo = grafo.subgraph(set(grafo.nodes) - {v})

            return verifica(novo_grafo, k, s - 1)


    # Verifica se cada vértice satisfaz as condições de Ramsey

    for v in list(grafo.nodes):

        #neighbors_v = set(grafo.neighbors(v))
        #non_neighbors_v = grafo.subgraph(set(grafo.nodes) - neighbors_v - {v})
        #print(v,neighbors_v)

        if grafo.degree[v] >= ramsey[k-1, s]:

            return False

        if n - grafo.degree[v]-1 >= ramsey[k, s-1]:

            return False

        neighbors = set(grafo.neighbors(v))

        induced_subgraph = grafo.subgraph(neighbors)

 

        if not verifica(induced_subgraph, k-1, s):

           return False 

        todos_vertices = set(grafo.nodes)

        nao_vizinhos_v = todos_vertices - neighbors - {v}

        subgrafo_nao_vizinhos = grafo.subgraph(nao_vizinhos_v)

        if not verifica(subgrafo_nao_vizinhos, k, s-1):

            return False
 

    # Se todas as condições forem satisfeitas para todos os vértices, retorna True

    return True


def cria_grafo(n, k, s, tempo_limite):
    
    contador_tentativas = 0
    inicio = time.time()

    while time.time() - inicio < tempo_limite:

        # Cria um grafo aleatório com n vértices, prob 1/2
        contador_tentativas += 1
        grafo = nx.erdos_renyi_graph(n, 0.5)

        # Verifica se o grafo satisfaz a propriedade de Ramsey para R(k,s)
        if verifica(grafo, k, s):
            '''
            # Exibindo o grafo
            nx.draw(grafo, with_labels=True, node_color='skyblue', font_color='black', font_weight='bold')
            plt.show()
            '''
            
            print(verifica(grafo, k, s))
            print(f"Número de tentativas: {contador_tentativas}")
            return grafo, contador_tentativas
        
        
    print(f"Falha em encontrar um grafo após {contador_tentativas} tentativas em {tempo_limite} segundos.")
    return None, contador_tentativas
        
def salvar_media_csv(n, k, s, media_tempo, media_tentativas, nome_arquivo='Resultados_Aleatorio.csv'):
    """Salva os parâmetros n, k, s, a média do tempo de execução e a média das tentativas em um arquivo CSV."""
    with open(nome_arquivo, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')  # Usando ponto e vírgula como delimitador
        # Escreve o cabeçalho
        writer.writerow(['n', 'k', 's', 'Média do Tempo de Execução (s)', 'Média de Tentativas'])
        # Escreve os valores
        writer.writerow([n, k, s, media_tempo, media_tentativas])
    
    print(f"Parâmetros e médias salvos em {nome_arquivo} com sucesso!")

    
def main():
    random.seed(time.time())
   
    
    parametros = [
        #(9,4,4)
        #(14,4,4) #resposta é 9
        (10,3,5) # resposta é 14
        #fiz ate n = 12
        #(range (5,13),4,4) #resaposta é 18
    ]
    """
    n = 5 #numero de vértices
    k = 3 # Parâmetro k da propriedade de Ramsey
    s = 3 # Parâmetro s da propriedade de Ramsey
    """
    x = 1 # Número de execuções
    tempo_limite = 3600*3 #30 minutos tempo limite para execução em segundos
    
    tempos_execucao = []  # Armazena o tempo de cada execução
    tentativas_por_execucao = []  # Armazena o número de tentativas em cada execução
    grafos_gerados = []  # Lista para armazenar os grafos geradoso

    melhor_grafo = None
    pior_grafo = None
    min_tentativas = float('inf')
    max_tentativas = 0
    
    nome_arquivo = 'aleatorio44.csv'
    with open(nome_arquivo, mode='a', newline='') as arquivo:
        escritor_csv = csv.writer(arquivo,delimiter=';')
        escritor_csv.writerow(['n', 'k', 's', 'Média Tempo (s)', 'Desvio Padrão Tempo (s)', 'Variância Tempo (s²)',
                               'Média Tentativas', 'Desvio Padrão Tentativas', 'Variância Tentativas', 'Melhor Tentativa', 'Pior Tentativa'])
        for parametro in parametros:
            if isinstance(parametro[0], range):
                n_values = parametro[0]
            else:
                n_values = [parametro[0]]

            k = parametro[1]
            s = parametro[2]
            

            for n in n_values:
                
                    for i in range(x):
                        print(f"Execução {i + 1}:")

                        inicio = time.time()
                        grafo_resultado, tentativas = cria_grafo(n, k, s, tempo_limite)
                        fim = time.time()
                        tempo_de_execucao = fim - inicio

                        if grafo_resultado is not None:
                            grafos_gerados.append((grafo_resultado, tentativas))
                            tempos_execucao.append(tempo_de_execucao)
                            tentativas_por_execucao.append(tentativas)

                            # Atualiza o melhor e o pior grafo
                            if tentativas < min_tentativas:
                                min_tentativas = tentativas
                                melhor_grafo = grafo_resultado

                            if tentativas > max_tentativas:
                                max_tentativas = tentativas
                                pior_grafo = grafo_resultado

                    # Calcula as estatísticas
                    media_tempo = np.mean(tempos_execucao) if tempos_execucao else 0
                    desvio_tempo = np.std(tempos_execucao) if tempos_execucao else 0
                    variancia_tempo = np.var(tempos_execucao) if tempos_execucao else 0

                    media_tentativas = np.mean(tentativas_por_execucao) if tentativas_por_execucao else 0
                    desvio_tentativas = np.std(tentativas_por_execucao) if tentativas_por_execucao else 0
                    variancia_tentativas = np.var(tentativas_por_execucao) if tentativas_por_execucao else 0
                    
                    escritor_csv.writerow([n, k, s,
                                           f"{media_tempo:.4f}", desvio_tempo, variancia_tempo,
                                 int(media_tentativas), int(desvio_tentativas), int(variancia_tentativas),
                                 min_tentativas, max_tentativas])
                    
                    # Exibe os resultados
                    print(f"\n--- Resultados {n}, {k}, {s} ---")
                    print(f"Média de tentativas: {media_tentativas}")
                    print(f"Desvio padrão das tentativas: {desvio_tentativas}")
                    print(f"Variância das tentativas: {variancia_tentativas}")
                    print(f"Tempo médio de execução: {media_tempo:.4f} segundos")
                    print(f"Desvio padrão do tempo de execução: {desvio_tempo:.4f} segundos")
                    print(f"Variância do tempo de execução: {variancia_tempo:.4f} segundos")
                    print(f"Grafo encontrado em {min_tentativas} tentativas melhor")
                    print(f"Grafo encontrado em {max_tentativas} tentativas pior")
               
            
if __name__ == "__main__":

    main()
