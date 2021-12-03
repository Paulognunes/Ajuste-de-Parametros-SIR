# Lista de Imports
import csv
import time
import matplotlib
import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc
import networkx as nx
import matplotlib.pyplot as plt
from random import random, sample
from sklearn.metrics import mean_squared_error


# from sklearn.metrics import mean_absolute_error


def converter_cromossomo_para_decimal(cromossomo):
    # Obtém o valor decimal referente ao cromossomo. Basicamente transforma o cromossomo "binário" em decimal.
    cromossomo.reverse()
    potencia = 0
    valor = 0
    for i in cromossomo:
        if i == '1':
            valor = valor + (2 ** potencia)
        potencia += 1
    cromossomo.reverse()
    return valor


def interpolar_valor(n, start1, stop1, start2, stop2):
    # Interpolação de valores. Converte um valor contido em um intervalo para outro intervalo.
    return ((n - start1) / (stop1 - start1)) * (stop2 - start2) + start2


def leitura_dados(nome_arquivo):
    # Leitura dos dados
    sucetiveis_sjdr = []
    infectados_sjdr = []
    recuperados_sjdr = []

    arq = open(nome_arquivo)
    dados_sjdr = csv.DictReader(arq, fieldnames=["S", "I", "R"])

    for row in dados_sjdr:
        sucetiveis_sjdr.append(int(row['S']))
        infectados_sjdr.append(int(row['I']))
        recuperados_sjdr.append(int(row['R']))

    arq.close()
    return sucetiveis_sjdr, infectados_sjdr, recuperados_sjdr


def simulacao(infectados_dados_reais, recuperados_dados_reais, grafo, dias, beta, gamma, nos_infectados,
              nos_recuperados):
    # Simula a difusão da doença e retorna o erro da simulação em relação aos dados.
    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Selecionando o modelo
    modelo = ep.SIRModel(grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nos infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    cfg.add_model_initial_configuration("Removed", nos_recuperados)

    modelo.set_initial_status(cfg)

    # Simulando execução
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(i['node_count'][1])
        recuperados_simulacao.append(i['node_count'][2])

    mse_i = round(mean_squared_error(infectados_dados_reais, infectados_simulacao), 5)
    mse_r = round(mean_squared_error(recuperados_dados_reais, recuperados_simulacao), 5)

    # mse_i = round(mean_absolute_error(infectados_dados_reais, infectados_simulacao), 5)
    # mse_r = round(mean_absolute_error(recuperados_dados_reais, recuperados_simulacao), 5)

    # Condicional para lidar com estouro de memória provocado pelo mean square error
    if mse_i < 0:
        mse_i = 10000000

    if mse_r < 0:
        mse_r = 10000000

    return mse_i, mse_r, sucetiveis_simulacao, infectados_simulacao, recuperados_simulacao


class Individuo:
    # Cada indivíduo é uma possível solução para o problema.

    def __init__(self, nro_bits, geracao=0):
        # Construtor da classe, responsável por inicializar aleatoriamente os parâmetros beta e gamma
        self.beta = ["0", "0"]
        self.gamma = ["0", "0"]
        self.nota_avaliacao_beta = 0
        self.nota_avaliacao_gamma = 0
        self.nota_avaliacao = 0
        self.sucetiveis_simulacao = []
        self.infectados_simulacao = []
        self.recuperados_simulacao = []
        self.geracao = geracao

        # Inicialização aleatória dos cromossomos beta e gamma
        for i in range(nro_bits - 2):
            if random() < 0.5:
                self.beta.append("0")
            else:
                self.beta.append("1")

        for i in range(nro_bits - 2):
            if random() < 0.5:
                self.gamma.append("0")
            else:
                self.gamma.append("1")

    def fitness(self, infectados_sjdr, recuperados_sjdr, nro_bits, grafo, dias, nos_infectados, nos_recuperados):
        """ Realização da simulação de difusão da doença na rede. Em seguida, os resultados da simulação são comparados
         com os dados reais. A nota do indivíduo será a média dos erros em relação aos infectados e recuperados"""

        # Obtém-se o valor decimal referente aos cromossomos beta e gamma.
        beta = converter_cromossomo_para_decimal(self.beta)
        gamma = converter_cromossomo_para_decimal(self.gamma)

        # Interpolamos esses valores para o intervalo entre 0 e 1.
        beta = interpolar_valor(beta, 0, ((2 ** nro_bits) - 1), 0, 1)
        gamma = interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1)

        # Simulação
        mse_beta, mse_gamma, suce, inf, rec = simulacao(infectados_sjdr, recuperados_sjdr, grafo, dias, beta, gamma,
                                                        nos_infectados, nos_recuperados)

        # A nota final é dada pela média dos erros de beta e gamma. Quanto menor esse erro, melhor é a solução!
        media = (mse_beta + mse_gamma) / 2

        self.nota_avaliacao_beta = mse_beta
        self.nota_avaliacao_gamma = mse_gamma

        self.sucetiveis_simulacao = suce.copy()
        self.infectados_simulacao = inf.copy()
        self.recuperados_simulacao = rec.copy()
        self.nota_avaliacao = round((1000000 / media), 5)

    def crossover(self, outro_individuo):
        # Realiza o crossover "de um único ponto" nos cromossomos beta e gamma.
        nro_bits = len(self.beta)
        corte_beta = round(random() * len(self.beta))

        # Crossover Beta.
        filho1_beta = outro_individuo.beta[0:corte_beta] + self.beta[corte_beta::]
        filho2_beta = self.beta[0:corte_beta] + outro_individuo.beta[corte_beta::]

        corte_gamma = round(random() * len(self.gamma))

        # Crossover Gamma.
        filho1_gamma = outro_individuo.gamma[0:corte_gamma] + self.gamma[corte_gamma::]
        filho2_gamma = self.gamma[0:corte_gamma] + outro_individuo.gamma[corte_gamma::]

        # Criação de dois novos indivíduos que herdaram os genes dos "pais".
        filhos = [Individuo(nro_bits, (self.geracao + 1)),
                  Individuo(nro_bits, (self.geracao + 1))]

        filhos[0].beta = filho1_beta
        filhos[1].beta = filho2_beta
        filhos[0].gamma = filho1_gamma
        filhos[1].gamma = filho2_gamma
        return filhos

    def mutacao(self, taxa_mutacao):
        # Realiza a mutação nos cromossomos.
        for i in range(len(self.beta)):
            if random() < taxa_mutacao:
                if self.beta[i] == '1':
                    self.beta[i] = '0'
                else:
                    self.beta[i] = '1'

            if random() < taxa_mutacao:
                if self.gamma[i] == '1':
                    self.gamma[i] = '0'
                else:
                    self.gamma[i] = '1'
        return self


class AlgoritmoGenetico:
    # Definição da classe AG, responsável por "resolver" o problema.

    def __init__(self, tamanho_populacao):
        # Construtor da classe.
        self.tamanho_populacao = tamanho_populacao
        self.populacao = []
        self.geracao = 0
        self.melhor_solucao = 0
        self.melhores_solucoes = []

    def inicializar_populacao(self, nro_bits):
        # Inicializa a população inicial. Basicamente gera soluções aleatórias, que serão "evoluidas" posteriormente.
        for i in range(self.tamanho_populacao):
            self.populacao.append(Individuo(nro_bits))

        self.melhor_solucao = self.populacao[0]

    def seleciona_pai(self, soma_avaliacao):
        # Simulação da roleta viciada
        pai = -1
        valor_sorteado = random() * soma_avaliacao
        soma = 0
        i = 0

        while i < len(self.populacao) and soma < valor_sorteado:
            soma += self.populacao[i].nota_avaliacao
            pai += 1
            i += 1

        return pai

    def soma_avaliacoes(self):
        soma = 0
        for individuo in self.populacao:
            soma += individuo.nota_avaliacao

        return soma

    def ordenar_populacao(self):
        """Ordena a população conforme a nota de avaliação. Como o objetivo é minimizar o erro,
        o "melhor" indivíduo indivíduo estará na posição zero do vetor de populacao."""
        self.populacao = sorted(self.populacao, key=lambda populacao: populacao.nota_avaliacao, reverse=True)

    def melhor_individuo(self, individuo):
        """ Compara se o indivíduo passado como parâmetro é melhor que o melhor indivíduo encontrado até o momento."
        Se sim, o melhor indivíduo passa a ser o indivíduo do parâmetro."""

        # Compara a nota desses indivíduos. Quanto menor a nota, melhor é a solução. Ela representa o erro da simulação
        if self.melhor_solucao.nota_avaliacao < individuo.nota_avaliacao:
            self.melhor_solucao = individuo
            self.melhores_solucoes.append(individuo)

    def visualizar_geracao(self, nro_bits):
        # Imprime o melhor indivíduo da geração. Ele sempre estará na posição zero.
        melhor = self.populacao[0]

        # Obtém-se o valor decimal referente aos cromossomos beta e gamma.
        beta = converter_cromossomo_para_decimal(melhor.beta)
        gamma = converter_cromossomo_para_decimal(melhor.gamma)

        # Interpolamos esses valores para o intervalo entre 0 e 1.
        beta = interpolar_valor(beta, 0, ((2 ** nro_bits) - 1), 0, 1)
        gamma = interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1)

        print(f'\nGeração: {self.populacao[0].geracao}'
              f'\nBeta : {beta}'
              f'\nGamma: {gamma}'
              f'\nAvaliação Beta: {melhor.nota_avaliacao_beta}'
              f'\nAvaliação Gamma: {melhor.nota_avaliacao_gamma}'
              f'\nNOTA: {melhor.nota_avaliacao}')

    def evoluir_solucao(self, taxa_mutacao, nro_bits, nro_geracoes, inf, rec, grafo, dias, nos_infectados,
                        nos_recuperados):
        # Gera a população inicial
        self.inicializar_populacao(nro_bits)

        # Avalia os individuos na população
        for individuo in self.populacao:
            individuo.fitness(inf, rec, nro_bits, grafo, dias, nos_infectados, nos_recuperados)

        self.ordenar_populacao()
        self.visualizar_geracao(nro_bits)

        # Execução do loop responsável pela convergência dos parâmetros
        for i in range(nro_geracoes):
            soma_avaliacao = self.soma_avaliacoes()
            nova_populacao = []

            # Os dois melhores indivíduos da população são selecionados como "pais" da nova população.
            for individuo_gerados in range(0, self.tamanho_populacao, 2):
                # Seleção dos pais
                pai1 = self.seleciona_pai(soma_avaliacao)
                pai2 = self.seleciona_pai(soma_avaliacao)

                # Aplicação do crossover, gerando novos filhos com os genes dos pais
                filhos = self.populacao[pai1].crossover(self.populacao[pai2])

                # Aplicação do processo de mutação nos filhos.
                nova_populacao.append(filhos[0].mutacao(taxa_mutacao))
                nova_populacao.append(filhos[1].mutacao(taxa_mutacao))

            # A população antiga é substituida pela nova população, que foi criada anteriormente.
            self.populacao = list(nova_populacao)

            for individuo in self.populacao:
                # Avaliação dos invidíduos (Fitness).
                individuo.fitness(inf, rec, nro_bits, grafo, dias, nos_infectados, nos_recuperados)

            self.ordenar_populacao()
            self.visualizar_geracao(nro_bits)

            melhor = self.populacao[0]
            self.melhor_individuo(melhor)


def main(nome_arq, probabilidade_criacao_aresta, nro_simul):
    inicio = time.time()

    # Lendo os dados e gerando o grafo
    suce, inf, rec = leitura_dados(nome_arq)

    # Variáveis de controle
    tamanho_populacao = 20
    nro_bits = 10
    nro_geracoes = 300
    taxa_mutacao = 0.05
    nos = 90497
    dias = len(suce)
    quantidade_rec = rec[0]
    quantidade_infec = inf[0]

    # Leitura do grafo do vinícius após a remoção das arestas
    grafo = nx.erdos_renyi_graph(nos, probabilidade_criacao_aresta)
    path = str(nro_simul) + "grafo.graphml"
    nx.write_graphml(grafo, path)

    # Sorteia um conjunto de nós que será os nós infectados e recuperados da simulação
    aux = sample(grafo.nodes, (quantidade_infec + quantidade_rec))

    # Seleciona aleatoriamente quais nós serão os nós infectados na simulação
    nos_infectados = sample(aux, quantidade_infec)
    nos_recuperados = []

    for i in aux:
        if i not in nos_infectados:
            nos_recuperados.append(i)

    deg = nx.degree(grafo)
    media_grau = 0

    for _, j in deg:
        media_grau += j

    media_grau = media_grau / nos
    ag = AlgoritmoGenetico(tamanho_populacao)
    ag.evoluir_solucao(taxa_mutacao, nro_bits, nro_geracoes, inf, rec, grafo, dias, nos_infectados, nos_recuperados)

    nome_arq_result = str(nro_simul) + '_resultados.txt'
    dados = open(nome_arq_result, 'w')
    dados.write('########################################\n')
    dados.write('            Configuração AG\n')
    dados.write('########################################\n')
    dados.write('\nTamanho População: ' + str(tamanho_populacao))
    dados.write('\nNúmero de Gerações: ' + str(nro_geracoes))
    dados.write('\nBits Cromossomo: ' + str(nro_bits))
    dados.write('\nDias: ' + str(dias))
    dados.write('\nTaxa de Mutação: ' + str(taxa_mutacao))
    dados.write('\nQuantidade de nós: ' + str(nos))
    dados.write('\nNós infectados: ' + str(nos_infectados))
    dados.write('\nProbabilidae de criação de aresta: ' + str(probabilidade_criacao_aresta))
    dados.write('\nNome do arquivo de amostra: ' + str(nome_arq))
    dados.write('\nGrau médio dos nós: ' + str(media_grau))
    dados.write('\n\n')

    dados.write('########################################\n')
    dados.write('     Evolução da solução encontrada\n')
    dados.write('########################################\n')

    for i in range(len(ag.melhores_solucoes)):
        beta = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].beta)
        gamma = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].gamma)
        beta = round(interpolar_valor(beta, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
        gamma = round(interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
        r0 = 0
        if gamma != 0:
            r0 = beta / gamma

        dados.write('\nGeração: ' + str(ag.melhores_solucoes[i].geracao))
        dados.write('\nBeta : ' + str(beta))
        dados.write('\nGamma: ' + str(gamma))
        dados.write('\nAvaliação Beta: ' + str(ag.melhores_solucoes[i].nota_avaliacao_beta))
        dados.write('\nAvaliação Gamma: ' + str(ag.melhores_solucoes[i].nota_avaliacao_gamma))
        dados.write('\nNOTA: ' + str(ag.melhores_solucoes[i].nota_avaliacao))
        dados.write('\nNúmero Básico de Reprodução (R0): ' + str(r0))
        dados.write('\nSucetiveis Simulação: ' + str(ag.melhores_solucoes[i].sucetiveis_simulacao))
        dados.write('\nInfectados Simulação: ' + str(ag.melhores_solucoes[i].infectados_simulacao))
        dados.write('\nRecuperados Simulação: ' + str(ag.melhores_solucoes[i].recuperados_simulacao))
        dados.write('\n')

    beta = converter_cromossomo_para_decimal(ag.melhor_solucao.beta)
    gamma = converter_cromossomo_para_decimal(ag.melhor_solucao.gamma)
    beta = round(interpolar_valor(beta, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
    gamma = round(interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1), 8)

    dados.write('\n\nMelhor Solução Encontrada:')
    dados.write('\nGeração: ' + str(ag.melhor_solucao.geracao))
    dados.write('\nBeta : ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nAvaliação Beta: ' + str(ag.melhor_solucao.nota_avaliacao_beta))
    dados.write('\nAvaliação Gamma: ' + str(ag.melhor_solucao.nota_avaliacao_gamma))
    dados.write('\nNOTA: ' + str(ag.melhor_solucao.nota_avaliacao))
    dados.write('\nSucetiveis Simulação: ' + str(ag.melhor_solucao.sucetiveis_simulacao))
    dados.write('\nInfectados Simulação: ' + str(ag.melhor_solucao.infectados_simulacao))
    dados.write('\nRecuperados Simulação: ' + str(ag.melhor_solucao.recuperados_simulacao))
    dados.write('\n\n')

    resultados = open("resultados_simulacao.txt", 'w')
    for i in range(dias):
        resultados.write(str(ag.melhor_solucao.sucetiveis_simulacao[i]) + ',' +
                         str(ag.melhor_solucao.infectados_simulacao[i]) + ',' +
                         str(ag.melhor_solucao.recuperados_simulacao[i]) + '\n')
    resultados.close()

    lista_de_dias = []
    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    # Plotando o gráfico de Infectados
    plt.plot(lista_de_dias, inf, label='Infectados Reais', linestyle='-')
    plt.plot(lista_de_dias, ag.melhor_solucao.infectados_simulacao, label='Infectados Simulados', linestyle='--')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Comparação Infectados")
    plt.legend()
    plt.close()
    # plt.show()

    # Salvando o gráfico
    nome_figura_inf = str(nro_simul) + '_Figura_Infectados.pdf'
    fig1.savefig(nome_figura_inf, format='pdf')

    # Plotando o gráfico de Recuperados
    plt.plot(lista_de_dias, rec, label='Recuperados Reais', linestyle='-')
    plt.plot(lista_de_dias, ag.melhor_solucao.recuperados_simulacao, label='Recuperados Simulados', linestyle='--')

    fig2 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Comparação Recuperados")
    plt.legend()
    plt.close()
    # plt.show()

    # Salvando o gráfico
    nome_figura_rec = str(nro_simul) + '_Figura_Recuperados.pdf'
    fig2.savefig(nome_figura_rec, format='pdf')

    fim = time.time()
    dados.write('Tempo de Execução: ' + str(fim - inicio) + ' segundos')
    dados.close()


if __name__ == '__main__':

    main("dados_sjdr_12_03_2021__22_04_2021.csv", 0.00007, 1)
    main("dados_sjdr_12_03_2021__22_04_2021.csv", 0.00007, 2)
    main("dados_sjdr_12_03_2021__22_04_2021.csv", 0.00007, 3)
    main("dados_sjdr_12_03_2021__22_04_2021.csv", 0.00007, 4)
    main("dados_sjdr_12_03_2021__22_04_2021.csv", 0.00007, 5)
