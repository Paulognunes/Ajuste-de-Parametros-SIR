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

    if dias == 0:
        return 0, 0, [0], [0], [0]

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
        self.beta1 = ["0", "0", "0"]
        self.beta2 = ["0", "0", "0"]
        self.gamma = ["0", "0", "0"]
        self.delta_t = []
        self.nota_avaliacao_beta1 = 0
        self.nota_avaliacao_gamma1 = 0
        self.nota_avaliacao_beta2 = 0
        self.nota_avaliacao_gamma2 = 0
        self.nota_avaliacao = 0
        self.sucetiveis_simulacao1 = []
        self.infectados_simulacao1 = []
        self.recuperados_simulacao1 = []
        self.sucetiveis_simulacao2 = []
        self.infectados_simulacao2 = []
        self.recuperados_simulacao2 = []
        self.geracao = geracao

        # Inicialização aleatória dos cromossomos beta e gamma
        for i in range(nro_bits - 3):
            # Beta1
            if random() < 0.5:
                self.beta1.append("0")
            else:
                self.beta1.append("1")

            # Beta2
            if random() < 0.5:
                self.beta2.append("0")
            else:
                self.beta2.append("1")

            # Gamma
            if random() < 0.5:
                self.gamma.append("0")
            else:
                self.gamma.append("1")

        for i in range(nro_bits):
            # Delta_t
            if random() < 0.5:
                self.delta_t.append("0")
            else:
                self.delta_t.append("1")

    def fitness(self, infectados_sjdr, recuperados_sjdr, nro_bits, grafo, dias, nos_infectados, nos_recuperados):
        """ Realização da simulação de difusão da doença na rede. Em seguida, os resultados da simulação são comparados
         com os dados reais. A nota do indivíduo será a média dos erros em relação aos infectados e recuperados"""

        # Obtém-se o valor decimal referente aos cromossomos beta, gamma e delta_t.
        beta1 = converter_cromossomo_para_decimal(self.beta1)
        beta2 = converter_cromossomo_para_decimal(self.beta2)
        gamma = converter_cromossomo_para_decimal(self.gamma)
        delta_t = converter_cromossomo_para_decimal(self.delta_t)

        # Interpolamos esses valores para o intervalo entre 0 e 1.
        beta1 = interpolar_valor(beta1, 0, ((2 ** nro_bits) - 1), 0, 1)
        beta2 = interpolar_valor(beta2, 0, ((2 ** nro_bits) - 1), 0, 1)
        gamma = interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1)
        delta_t = round(interpolar_valor(delta_t, 0, ((2 ** nro_bits) - 1), 1, dias - 1))

        # Simulação 1
        mse_beta1, mse_gamma1, suce1, inf1, rec1 = simulacao(infectados_sjdr[0:delta_t], recuperados_sjdr[0:delta_t],
                                                             grafo, delta_t, beta1, gamma, nos_infectados,
                                                             nos_recuperados)

        if delta_t == 0:
            aux = sample(grafo.nodes, (infectados_sjdr[0] + infectados_sjdr[0]))
            nos_infectados2 = sample(aux, infectados_sjdr[delta_t])
            nos_recuperados2 = []

            for i in aux:
                if i not in nos_infectados2:
                    nos_recuperados2.append(i)

        else:
            # Selecionar nós a serem infectados/recuperados na segunda simulação
            aux = sample(grafo.nodes, (inf1[-1] + rec1[-1]))
            nos_infectados2 = sample(aux, inf1[-1])
            nos_recuperados2 = []

            for i in aux:
                if i not in nos_infectados2:
                    nos_recuperados2.append(i)

        if delta_t != dias:
            # Simulação 2
            mse_beta2, mse_gamma2, suce2, inf2, rec2 = simulacao(infectados_sjdr[delta_t::],
                                                                 recuperados_sjdr[delta_t::],
                                                                 grafo, (dias - delta_t), beta2, gamma, nos_infectados2,
                                                                 nos_recuperados2)
        else:
            mse_beta2 = 0
            mse_gamma2 = 0
            suce2 = []
            inf2 = []
            rec2 = []

        # A nota final é dada pela média dos erros de beta e gamma. Quanto menor esse erro, melhor é a solução!
        media = (mse_beta1 + mse_gamma1 + mse_beta2 + mse_gamma2) / 4

        self.nota_avaliacao_beta1 = mse_beta1
        self.nota_avaliacao_gamma1 = mse_gamma1

        self.nota_avaliacao_beta2 = mse_beta2
        self.nota_avaliacao_gamma2 = mse_gamma2

        self.sucetiveis_simulacao1 = suce1.copy()
        self.infectados_simulacao1 = inf1.copy()
        self.recuperados_simulacao1 = rec1.copy()

        self.sucetiveis_simulacao2 = suce2.copy()
        self.infectados_simulacao2 = inf2.copy()
        self.recuperados_simulacao2 = rec2.copy()

        self.nota_avaliacao = round((1000000 / media), 8)

    def crossover(self, outro_individuo):
        # Realiza o crossover "de um único ponto" nos cromossomos beta e gamma.
        nro_bits = len(self.beta1)

        # Crossover Beta1
        corte_beta = round(random() * len(self.beta1))
        filho1_beta1 = outro_individuo.beta1[0:corte_beta] + self.beta1[corte_beta::]
        filho2_beta1 = self.beta1[0:corte_beta] + outro_individuo.beta1[corte_beta::]

        # Crossover Beta2
        corte_beta = round(random() * len(self.beta1))
        filho1_beta2 = outro_individuo.beta2[0:corte_beta] + self.beta2[corte_beta::]
        filho2_beta2 = self.beta2[0:corte_beta] + outro_individuo.beta2[corte_beta::]

        # Crossover Gamma
        corte_gamma = round(random() * len(self.gamma))
        filho1_gamma = outro_individuo.gamma[0:corte_gamma] + self.gamma[corte_gamma::]
        filho2_gamma = self.gamma[0:corte_gamma] + outro_individuo.gamma[corte_gamma::]

        # Crossover delta_t
        corte_delta_t = round(random() * len(self.delta_t))
        filho1_delta_t = outro_individuo.delta_t[0:corte_delta_t] + self.delta_t[corte_delta_t::]
        filho2_delta_t = self.delta_t[0:corte_delta_t] + outro_individuo.delta_t[corte_delta_t::]

        # Criação de dois novos indivíduos que herdaram os genes dos "pais".
        filhos = [Individuo(nro_bits, (self.geracao + 1)),
                  Individuo(nro_bits, (self.geracao + 1))]

        filhos[0].beta1 = filho1_beta1
        filhos[0].beta2 = filho1_beta2
        filhos[0].gamma = filho1_gamma
        filhos[0].delta_t = filho1_delta_t

        filhos[1].beta1 = filho2_beta1
        filhos[1].beta2 = filho2_beta2
        filhos[1].gamma = filho2_gamma
        filhos[1].delta_t = filho2_delta_t

        return filhos

    def mutacao(self, taxa_mutacao):
        # Realiza a mutação nos cromossomos.
        for i in range(len(self.beta1)):
            if random() < taxa_mutacao:
                if self.beta1[i] == '1':
                    self.beta1[i] = '0'
                else:
                    self.beta1[i] = '1'

            if random() < taxa_mutacao:
                if self.beta2[i] == '1':
                    self.beta2[i] = '0'
                else:
                    self.beta2[i] = '1'

            if random() < taxa_mutacao:
                if self.gamma[i] == '1':
                    self.gamma[i] = '0'
                else:
                    self.gamma[i] = '1'

            if random() < taxa_mutacao:
                if self.delta_t[i] == '1':
                    self.delta_t[i] = '0'
                else:
                    self.delta_t[i] = '1'
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

    def visualizar_geracao(self, nro_bits, dias):
        # Imprime o melhor indivíduo da geração. Ele sempre estará na posição zero.
        melhor = self.populacao[0]

        # Obtém-se o valor decimal referente aos cromossomos beta e gamma.
        beta1 = converter_cromossomo_para_decimal(melhor.beta1)
        beta2 = converter_cromossomo_para_decimal(melhor.beta2)
        gamma = converter_cromossomo_para_decimal(melhor.gamma)
        delta_t = converter_cromossomo_para_decimal(melhor.delta_t)

        # Interpolamos esses valores para o intervalo entre 0 e 1.
        beta1 = interpolar_valor(beta1, 0, ((2 ** nro_bits) - 1), 0, 1)
        beta2 = interpolar_valor(beta2, 0, ((2 ** nro_bits) - 1), 0, 1)
        gamma = interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1)
        delta_t = round(interpolar_valor(delta_t, 0, ((2 ** nro_bits) - 1), 1, dias))

        print(f'\nGeração: {self.populacao[0].geracao}'
              f'\nBeta1 : {beta1}'
              f'\nBeta2 : {beta2}'
              f'\nGamma: {gamma}'
              f'\nDelta T : {delta_t}'
              f'\nAvaliação Beta1: {melhor.nota_avaliacao_beta1}'
              f'\nAvaliação Gamma1: {melhor.nota_avaliacao_gamma1}'
              f'\nAvaliação Beta2: {melhor.nota_avaliacao_beta2}'
              f'\nAvaliação Gamma21: {melhor.nota_avaliacao_gamma2}'
              f'\nNOTA: {melhor.nota_avaliacao}')

    def evoluir_solucao(self, taxa_mutacao, nro_bits, nro_geracoes, inf, rec, grafo, dias, nos_infectados,
                        nos_recuperados):
        # Gera a população inicial
        self.inicializar_populacao(nro_bits)

        # Avalia os individuos na população
        for individuo in self.populacao:
            individuo.fitness(inf, rec, nro_bits, grafo, dias, nos_infectados, nos_recuperados)

        self.ordenar_populacao()
        self.visualizar_geracao(nro_bits, dias)

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
            self.visualizar_geracao(nro_bits, dias)

            melhor = self.populacao[0]
            self.melhor_individuo(melhor)


def main(nome_arq, probabilidade_criacao_aresta, nro_simul):
    inicio = time.time()

    # Lendo os dados e gerando o grafo
    # nome_arq = 'dados_sjdr_18_05_2021__15_09_2021.csv'
    suce, inf, rec = leitura_dados(nome_arq)

    # Variáveis de controle
    tamanho_populacao = 20
    nro_bits = 10
    nro_geracoes = 150
    taxa_mutacao = 0.1
    nos = 90497
    dias = len(suce)
    quantidade_rec = rec[0]
    quantidade_infec = inf[0]

    # grafo = nx.erdos_renyi_graph(nos, probabilidade_criacao_aresta)
    # path = "grafo" + "_" + str(nro_simul) + ".graphml"
    # nx.write_graphml(grafo, path)
    grafo = nx.read_graphml("grafo_1.graphml")

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
    # print(media_grau)
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
        beta1 = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].beta1)
        beta2 = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].beta2)
        gamma = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].gamma)
        delta_t = converter_cromossomo_para_decimal(ag.melhores_solucoes[i].delta_t)

        beta1 = round(interpolar_valor(beta1, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
        beta2 = round(interpolar_valor(beta2, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
        gamma = round(interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
        delta_t = round(interpolar_valor(delta_t, 0, ((2 ** nro_bits) - 1), 1, dias))

        r0 = 0
        if gamma != 0:
            r0 = (beta1 + beta2) / 2 * gamma

        dados.write('\nGeração: ' + str(ag.melhores_solucoes[i].geracao))
        dados.write('\nBeta1 : ' + str(beta1))
        dados.write('\nBeta2 : ' + str(beta2))
        dados.write('\nGamma: ' + str(gamma))
        dados.write('\nDelta t: ' + str(delta_t))
        dados.write('\nAvaliação Beta1: ' + str(ag.melhores_solucoes[i].nota_avaliacao_beta1))
        dados.write('\nAvaliação Gamma1: ' + str(ag.melhores_solucoes[i].nota_avaliacao_gamma1))
        dados.write('\nAvaliação Beta2: ' + str(ag.melhores_solucoes[i].nota_avaliacao_beta2))
        dados.write('\nAvaliação Gamma2: ' + str(ag.melhores_solucoes[i].nota_avaliacao_gamma2))
        dados.write('\nNOTA: ' + str(ag.melhores_solucoes[i].nota_avaliacao))
        dados.write('\nNúmero Básico de Reprodução (R0): ' + str(r0))
        dados.write('\nSucetiveis Simulação1: ' + str(ag.melhores_solucoes[i].sucetiveis_simulacao1))
        dados.write('\nInfectados Simulação1: ' + str(ag.melhores_solucoes[i].infectados_simulacao1))
        dados.write('\nRecuperados Simulação1: ' + str(ag.melhores_solucoes[i].recuperados_simulacao1))
        dados.write('\nSucetiveis Simulação2: ' + str(ag.melhores_solucoes[i].sucetiveis_simulacao2))
        dados.write('\nInfectados Simulação2: ' + str(ag.melhores_solucoes[i].infectados_simulacao2))
        dados.write('\nRecuperados Simulação2: ' + str(ag.melhores_solucoes[i].recuperados_simulacao2))

        dados.write('\n')

    beta1 = converter_cromossomo_para_decimal(ag.melhor_solucao.beta1)
    beta2 = converter_cromossomo_para_decimal(ag.melhor_solucao.beta2)
    gamma = converter_cromossomo_para_decimal(ag.melhor_solucao.gamma)
    delta_t = converter_cromossomo_para_decimal(ag.melhor_solucao.delta_t)

    beta1 = round(interpolar_valor(beta1, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
    beta2 = round(interpolar_valor(beta2, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
    gamma = round(interpolar_valor(gamma, 0, ((2 ** nro_bits) - 1), 0, 1), 8)
    delta_t = round(interpolar_valor(delta_t, 0, ((2 ** nro_bits) - 1), 1, dias))

    dados.write('\n\nMelhor Solução Encontrada:')
    dados.write('\nGeração: ' + str(ag.melhor_solucao.geracao))
    dados.write('\nBeta1 : ' + str(beta1))
    dados.write('\nBeta2 : ' + str(beta2))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nDelta t : ' + str(delta_t))
    dados.write('\nAvaliação Beta1: ' + str(ag.melhor_solucao.nota_avaliacao_beta1))
    dados.write('\nAvaliação Gamma1: ' + str(ag.melhor_solucao.nota_avaliacao_gamma1))
    dados.write('\nAvaliação Beta2: ' + str(ag.melhor_solucao.nota_avaliacao_beta1))
    dados.write('\nAvaliação Gamma2: ' + str(ag.melhor_solucao.nota_avaliacao_gamma1))
    dados.write('\nNOTA: ' + str(ag.melhor_solucao.nota_avaliacao))
    dados.write('\nSucetiveis Simulação1: ' + str(ag.melhor_solucao.sucetiveis_simulacao1))
    dados.write('\nInfectados Simulação1: ' + str(ag.melhor_solucao.infectados_simulacao1))
    dados.write('\nRecuperados Simulação1: ' + str(ag.melhor_solucao.recuperados_simulacao1))
    dados.write('\nSucetiveis Simulação2: ' + str(ag.melhor_solucao.sucetiveis_simulacao2))
    dados.write('\nInfectados Simulação2: ' + str(ag.melhor_solucao.infectados_simulacao2))
    dados.write('\nRecuperados Simulação2: ' + str(ag.melhor_solucao.recuperados_simulacao2))
    dados.write('\n\n')

    lista_de_dias = []
    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    list_infec1 = []
    list_infec1 += (ag.melhor_solucao.infectados_simulacao1 + ag.melhor_solucao.infectados_simulacao2)

    # Plotando o gráfico de Infectados
    plt.plot(lista_de_dias, inf, label='Infectados Reais', linestyle='-')

    plt.plot(lista_de_dias, list_infec1, label='Infectados Simulados', linestyle='--')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Comparação Infectados")
    plt.legend()
    plt.close()
    # plt.show()

    # Salvando o gráfico
    nome_figura_inf = str(nro_simul) + '_Figura_Infectados.png'
    fig1.savefig(nome_figura_inf, format='png')

    list_infec2 = []
    list_infec2 += (ag.melhor_solucao.recuperados_simulacao1 + ag.melhor_solucao.recuperados_simulacao2)

    # Plotando o gráfico de Recuperados
    plt.plot(lista_de_dias, rec, label='Recuperados Reais', linestyle='-')
    plt.plot(lista_de_dias, list_infec2, label='Recuperados Simulados', linestyle='--')

    fig2 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Comparação Recuperados")
    plt.legend()
    plt.close()
    # plt.show()

    # Salvando o gráfico
    nome_figura_rec = str(nro_simul) + '_Figura_Recuperados.png'
    fig2.savefig(nome_figura_rec, format='png')

    fim = time.time()
    dados.write('Tempo de Execução: ' + str(fim - inicio) + ' segundos')
    dados.close()


if __name__ == '__main__':
    # Configuração 1
    main("dados_sjdr_15_07_2021__15_09_2021.csv", 0.00007, 1)
