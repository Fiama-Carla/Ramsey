# Algoritmos para Busca de Contraexemplos no Problema de Ramsey

Este repositório contém o código-fonte desenvolvido como parte do **Trabalho de Conclusão de Curso (TCC)**, cujo objetivo é a aplicação de **heurísticas na busca por contraexemplos no problema de Ramsey**. O estudo envolve algoritmos **aleatórios e genéticos** para a coloração de arestas e a identificação de grafos monocromáticos.

## 📌 Descrição
O objetivo deste projeto é explorar diferentes heurísticas para encontrar contraexemplos dos números de Ramsey. A implementação considera algoritmos aleatórios e genéticos, aplicando estratégias de otimização baseadas em reprodução e seleção de alelos. Foram implementadas duas abordagens principais:

- **🔹 Algoritmo Aleatório:** Geração e teste de grafos aleatórios para verificar se satisfazem as condições do problema de Ramsey.
- **🔹 Algoritmo Genético:** Evolução de populações de grafos com operadores de seleção, crossover e mutação, explorando a estrutura de células e alelos para melhorar a busca por contraexemplos.

## 🚀 Como Executar

1. **Clone o repositório**:
   ```bash
   git clone https://github.com/Fiama-Carla/Ramsey.git

2. **Acesse a pasta do projeto**:
   ```bash
   cd repositorio
   ```
3. **Instale as dependências **:
   ```bash
   pip install -r requirements.txt
   ```
4. **Execute o código**:
   ```bash
   python src/main.py
   ```

## 📊 Resultados e Análises

Os experimentos realizados testaram a eficácia dos algoritmos para encontrar contraexemplos nos seguintes casos:
- **R(3,4)  (limite conhecido: 9)**
- **R(4,4)  (limite conhecido: 18)**
- **R(3,5)  (exploração de novos limites)**

As métricas analisadas incluem:
- Tempo médio de execução
- Número médio de tentativas
- Melhor e pior desempenho dos algoritmos
- Média móvel do fitness ao log das gerações


## 🛠️ Tecnologias Utilizadas
- **Python** 🐍
- **NetworkX** para manipulação de grafos
- **Bibliotecas matemáticas** (NumPy, itertools, etc.)

## 🐝 Licença
Este projeto é de código aberto e pode ser utilizado para fins acadêmicos e de pesquisa.

## 📌 Contato
Para dúvidas ou sugestões, entre em contato através do e-mail: **fiamacmsousa@gmail.com**

