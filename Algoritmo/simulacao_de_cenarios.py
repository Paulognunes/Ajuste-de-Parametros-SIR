import matplotlib
import ndlib.models.epidemics as ep
import ndlib.models.ModelConfig as mc
import networkx as nx
import matplotlib.pyplot as plt
from random import sample


def simulacao_simples(nome_grafo, beta, gamma, dias, qtd_nos_infectados, nome_fig):
    # Simulação do processo de propagação do vírus na rede

    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Leitura do grafo
    grafo = nx.read_graphml(nome_grafo)

    # Seleção dos nós iniciais da simulação
    nos_infectados = sample(grafo.nodes, qtd_nos_infectados)

    # Selecionando o modelo
    modelo = ep.SIRModel(grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nós infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    # cfg.add_model_initial_configuration("Removed", nos_recuperados)

    # Simulando execução
    modelo.set_initial_status(cfg)
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(int(i['node_count'][1]))
        recuperados_simulacao.append(i['node_count'][2])

    max_infect = max(infectados_simulacao)
    dia_max_infect = infectados_simulacao.index(max_infect)

    # Visualização
    lista_de_dias = []

    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    plt.plot(lista_de_dias, sucetiveis_simulacao, label='Suscetíveis', color='y', linestyle='-')
    plt.plot(lista_de_dias, infectados_simulacao, label='Infectados', color='r', linestyle='-')
    plt.plot(lista_de_dias, recuperados_simulacao, label='Recuperados', color='g', linestyle='-')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Simulação Convencional")
    plt.legend()
    plt.close()
    # plt.show()

    fig1.savefig(nome_fig, format='pdf')

    dados = open("dados_simulacao.txt", 'w')
    dados.write('########################################\n')
    dados.write('            Simulação\n')
    dados.write('########################################\n')
    dados.write('\nBeta: ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nDias: ' + str(dias))
    dados.write('\nNúmero Máximo de Infectados: ' + str(max_infect))
    dados.write('\nDia Ocorrido: ' + str(dia_max_infect + 1))
    dados.write('\n\nSuscetíveis simulação: ' + str(sucetiveis_simulacao))
    dados.write('\n\nInfectados simulação: ' + str(infectados_simulacao))
    dados.write('\n\nRecuperados simulação: ' + str(recuperados_simulacao))
    dados.write('\n\n')
    dados.close()

    resultados = open("resultados_simulacao_normal.txt", 'w')
    for i in range(dias):
        resultados.write(str(sucetiveis_simulacao[i]) + ',' +
                         str(infectados_simulacao[i]) + ',' +
                         str(recuperados_simulacao[i]) + '\n')
    resultados.close()
    return nos_infectados


def simulacao_isolamento_horizontal(nome_grafo, beta, gamma, dias, nos_infectados, porcent_rem_arestas,
                                    nome_fig):
    # Simulação do processo de propagação do vírus em um cenário de isolamento horizontal

    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Modificação da rede
    path = "grafo_" + str(100 * porcent_rem_arestas) + "_removed.graphml"
    new_grafo = remover_arestas(nome_grafo, porcent_rem_arestas, path)

    # Selecionando o modelo
    modelo = ep.SIRModel(new_grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nós infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    # cfg.add_model_initial_configuration("Removed", nos_recuperados)

    # Simulando execução
    modelo.set_initial_status(cfg)
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(int(i['node_count'][1]))
        recuperados_simulacao.append(i['node_count'][2])

    max_infect = max(infectados_simulacao)
    dia_max_infect = infectados_simulacao.index(max_infect)

    # Visualização
    lista_de_dias = []

    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    plt.plot(lista_de_dias, sucetiveis_simulacao, label='Suscetíveis', color='y', linestyle='-')
    plt.plot(lista_de_dias, infectados_simulacao, label='Infectados', color='r', linestyle='-')
    plt.plot(lista_de_dias, recuperados_simulacao, label='Recuperados', color='g', linestyle='-')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Isolamento Horizontal")
    plt.legend()
    plt.close()
    # plt.show()

    fig1.savefig(nome_fig, format='pdf')

    dados = open("dados_simulacao_isolamento_horizontal.txt", 'w')
    dados.write('########################################\n')
    dados.write('   Simulação de Isolamento Horizontal   \n')
    dados.write('########################################\n')
    dados.write('\nBeta: ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nDias: ' + str(dias))
    dados.write('\nNúmero Máximo de Infectados: ' + str(max_infect))
    dados.write('\nDia Ocorrido: ' + str(dia_max_infect + 1))
    dados.write('\n\nSuscetíveis simulação: ' + str(sucetiveis_simulacao))
    dados.write('\n\nInfectados simulação: ' + str(infectados_simulacao))
    dados.write('\n\nRecuperados simulação: ' + str(recuperados_simulacao))
    dados.write('\n\n')
    dados.close()

    resultados = open("resultados_simulacao_isolamento_vertical.txt", 'w')
    for i in range(dias):
        resultados.write(str(sucetiveis_simulacao[i]) + ',' +
                         str(infectados_simulacao[i]) + ',' +
                         str(recuperados_simulacao[i]) + '\n')
    resultados.close()


def simulacao_aumento_de_interacoes_sociais(nome_grafo, beta, gamma, dias, nos_infectados, porcent_add_arestas,
                                            nome_fig):
    # Simulação do processo de propagação do vírus em um cenário de aumento das interacoes sociais

    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Modificação da rede
    path = "grafo_" + str(100 * porcent_add_arestas) + "_add.graphml"
    new_grafo = add_arestas(nome_grafo, porcent_add_arestas, path)

    # Selecionando o modelo
    modelo = ep.SIRModel(new_grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nós infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    # cfg.add_model_initial_configuration("Removed", nos_recuperados)

    # Simulando execução
    modelo.set_initial_status(cfg)
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(int(i['node_count'][1]))
        recuperados_simulacao.append(i['node_count'][2])

    max_infect = max(infectados_simulacao)
    dia_max_infect = infectados_simulacao.index(max_infect)

    # Visualização
    lista_de_dias = []

    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    plt.plot(lista_de_dias, sucetiveis_simulacao, label='Suscetíveis', color='y', linestyle='-')
    plt.plot(lista_de_dias, infectados_simulacao, label='Infectados', color='r', linestyle='-')
    plt.plot(lista_de_dias, recuperados_simulacao, label='Recuperados', color='g', linestyle='-')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Relaxamento do Distanciamento Social")
    plt.legend()
    plt.close()
    # plt.show()

    fig1.savefig(nome_fig, format='pdf')

    dados = open("dados_simulacao_relaxamento.txt", 'w')
    dados.write('########################################\n')
    dados.write('        Simulação de Relaxamento        \n')
    dados.write('########################################\n')
    dados.write('\nBeta: ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nDias: ' + str(dias))
    dados.write('\nNúmero Máximo de Infectados: ' + str(max_infect))
    dados.write('\nDia Ocorrido: ' + str(dia_max_infect + 1))
    dados.write('\n\nSuscetíveis simulação: ' + str(sucetiveis_simulacao))
    dados.write('\n\nInfectados simulação: ' + str(infectados_simulacao))
    dados.write('\n\nRecuperados simulação: ' + str(recuperados_simulacao))
    dados.write('\n\n')
    dados.close()

    resultados = open("resultados_simulacao_relaxamento.txt", 'w')
    for i in range(dias):
        resultados.write(str(sucetiveis_simulacao[i]) + ',' +
                         str(infectados_simulacao[i]) + ',' +
                         str(recuperados_simulacao[i]) + '\n')
    resultados.close()


def simulacao_vacinacao_aleatoria(nome_grafo, beta, gamma, dias, nos_infectados, porcentagem_vacinados, nome_fig):
    # Simulação do processo de propagação do vírus em uma população parcialmente vacinada.
    # A estratégia de vacinação desse cenário é "Vacinação Aleatória"

    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Leitura do grafo
    grafo = nx.read_graphml(nome_grafo)

    # Seleção dos nós iniciais da simulação
    aux = list(set(grafo.nodes) - set(nos_infectados))
    nos_recuperados = sample(aux, round(porcentagem_vacinados * nx.number_of_nodes(grafo)))

    # Selecionando o modelo
    modelo = ep.SIRModel(grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nós infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    cfg.add_model_initial_configuration("Removed", nos_recuperados)

    # Simulando execução
    modelo.set_initial_status(cfg)
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(int(i['node_count'][1]))
        recuperados_simulacao.append(i['node_count'][2])

    max_infect = max(infectados_simulacao)
    dia_max_infect = infectados_simulacao.index(max_infect)

    # Visualização
    lista_de_dias = []

    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    plt.plot(lista_de_dias, sucetiveis_simulacao, label='Suscetíveis', color='y', linestyle='-')
    plt.plot(lista_de_dias, infectados_simulacao, label='Infectados', color='r', linestyle='-')
    plt.plot(lista_de_dias, recuperados_simulacao, label='Recuperados', color='g', linestyle='-')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Simulação de Vacinação Aleatória")
    plt.legend()
    plt.close()
    # plt.show()

    fig1.savefig(nome_fig, format='pdf')

    dados = open("dados_simulacao_vacinacao_aleatoria.txt", 'w')
    dados.write('########################################\n')
    dados.write('    Simulação de Vacinação Aleatória    \n')
    dados.write('########################################\n')
    dados.write('\nBeta: ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nPorcentagem de vacinados: ' + str(100 * porcentagem_vacinados) + '%')
    dados.write('\nDias: ' + str(dias))
    dados.write('\nNúmero Máximo de Infectados: ' + str(max_infect))
    dados.write('\nDia Ocorrido: ' + str(dia_max_infect + 1))
    dados.write('\n\nSuscetíveis simulação: ' + str(sucetiveis_simulacao))
    dados.write('\n\nInfectados simulação: ' + str(infectados_simulacao))
    dados.write('\n\nRecuperados simulação: ' + str(recuperados_simulacao))
    dados.write('\n\n')
    dados.close()

    resultados = open("resultados_simulacao_vacinacao_aleatoria.txt", 'w')
    for i in range(dias):
        resultados.write(str(sucetiveis_simulacao[i]) + ',' +
                         str(infectados_simulacao[i]) + ',' +
                         str(recuperados_simulacao[i]) + '\n')
    resultados.close()


def simulacao_vacinacao_hubs(nome_grafo, beta, gamma, dias, qtd_nos_infectados, porcentagem_vacinados, nome_fig):
    # Simulação do processo de propagação do vírus em uma população parcialmente vacinada.
    # A estratégia de vacinação desse cenário é "Vacinação dos Hubs"

    sucetiveis_simulacao = []
    infectados_simulacao = []
    recuperados_simulacao = []

    # Leitura do grafo
    grafo = nx.read_graphml(nome_grafo)

    # Seleção dos nós iniciais da simulação
    nos_recuperados = hubs_para_vacinacao(nome_grafo, porcentagem_vacinados)
    aux = list(set(grafo.nodes) - set(nos_recuperados))
    nos_infectados = sample(aux, qtd_nos_infectados)

    # Selecionando o modelo
    modelo = ep.SIRModel(grafo)

    # Configurando o modelo
    cfg = mc.Configuration()
    cfg.add_model_parameter('beta', beta)
    cfg.add_model_parameter('gamma', gamma)

    # Inicializando os nós infectados
    cfg.add_model_initial_configuration("Infected", nos_infectados)
    cfg.add_model_initial_configuration("Removed", nos_recuperados)

    # Simulando execução
    modelo.set_initial_status(cfg)
    iteracoes = modelo.iteration_bunch(dias)

    # Obtendo o número de suscetiveis, infectados e recuperados a cada iteração
    for i in iteracoes:
        sucetiveis_simulacao.append(i['node_count'][0])
        infectados_simulacao.append(int(i['node_count'][1]))
        recuperados_simulacao.append(i['node_count'][2])

    max_infect = max(infectados_simulacao)
    dia_max_infect = infectados_simulacao.index(max_infect)

    # Visualização
    lista_de_dias = []

    for i in range(1, dias + 1):
        lista_de_dias.append(i)

    plt.plot(lista_de_dias, sucetiveis_simulacao, label='Suscetíveis', color='y', linestyle='-')
    plt.plot(lista_de_dias, infectados_simulacao, label='Infectados', color='r', linestyle='-')
    plt.plot(lista_de_dias, recuperados_simulacao, label='Recuperados', color='g', linestyle='-')

    fig1 = plt.gcf()
    matplotlib.pyplot.xlabel('Dias')
    matplotlib.pyplot.ylabel('Número de Pessoas')
    plt.title("Simulação de Vacinação dos Hubs")
    plt.legend()
    plt.close()
    # plt.show()

    fig1.savefig(nome_fig, format='pdf')

    dados = open("dados_simulacao_vacinacao_hubs.txt", 'w')
    dados.write('########################################\n')
    dados.write('    Simulação de Vacinação ods Hubs     \n')
    dados.write('########################################\n')
    dados.write('\nBeta: ' + str(beta))
    dados.write('\nGamma: ' + str(gamma))
    dados.write('\nPorcentagem de vacinados: ' + str(100 * porcentagem_vacinados) + '%')
    dados.write('\nDias: ' + str(dias))
    dados.write('\nNúmero Máximo de Infectados: ' + str(max_infect))
    dados.write('\nDia Ocorrido: ' + str(dia_max_infect + 1))
    dados.write('\n\nSuscetíveis simulação: ' + str(sucetiveis_simulacao))
    dados.write('\n\nInfectados simulação: ' + str(infectados_simulacao))
    dados.write('\n\nRecuperados simulação: ' + str(recuperados_simulacao))
    dados.write('\n\n')
    dados.close()

    resultados = open("resultados_simulacao_vacinacao_hubs.txt", 'w')
    for i in range(dias):
        resultados.write(str(sucetiveis_simulacao[i]) + ',' +
                         str(infectados_simulacao[i]) + ',' +
                         str(recuperados_simulacao[i]) + '\n')
    resultados.close()
    return nos_infectados


def remover_arestas(nome_grafo, porcentagem, path):
    grafo = nx.read_graphml(nome_grafo)

    arestas_removidas = sample(grafo.edges(), round(porcentagem * nx.number_of_edges(grafo)))

    for i in arestas_removidas:
        grafo.remove_edge(*i)

    nx.write_graphml(grafo, path)
    return grafo


def add_arestas(nome_grafo, porcentagem, path):
    grafo = nx.read_graphml(nome_grafo)
    qnt_arestas_add = round(porcentagem * nx.number_of_edges(grafo))
    i = 0

    while i < qnt_arestas_add:
        nos = sample(grafo.nodes(), 2)

        if not (grafo.has_edge(nos[0], nos[1]) or grafo.has_edge(nos[1], nos[0])):
            grafo.add_edge(nos[0], nos[1])
            i = i + 1

    nx.write_graphml(grafo, path)
    return grafo


def hubs_para_vacinacao(nome_grafo, porcentagem_vacinada):
    grafo = nx.read_graphml(nome_grafo)
    deg = nx.degree(grafo)
    lista_nos = []

    for i, j in deg:
        lista_nos.append((i, j))

    lista_nos.sort(key=lambda x: x[1], reverse=True)

    lista_de_hubs = []
    qnt_vacinados = round(porcentagem_vacinada * nx.number_of_nodes(grafo))

    for i in range(qnt_vacinados):
        (node, _) = lista_nos[i]
        lista_de_hubs.append(node)

    return lista_de_hubs


if __name__ == '__main__':
    nos_infectados = []
    nos_recuperados = []

    nos_infectados = simulacao_simples("4grafo.graphml", 0.01173, 0.04007, 365, 685, "simulacao_simples.pdf")

    simulacao_aumento_de_interacoes_sociais("4grafo.graphml", 0.01173, 0.04007, 365, nos_infectados, 0.4,
                                            "simulacao_relaxamento.pdf")

    simulacao_isolamento_horizontal("4grafo.graphml", 0.01173, 0.04007, 365, nos_infectados, 0.4,
                                    "simulacao_isolamento_horizontal.pdf")

    nos_infectados = simulacao_vacinacao_hubs("4grafo.graphml", 0.01173, 0.04007, 365, 685, 0.5,
                                              "simulacao_vacinacao_hubs.pdf")

    simulacao_vacinacao_aleatoria("4grafo.graphml", 0.01173, 0.04007, 365, nos_infectados, 0.5,
                                  "simulacao_vacinacao_aleatoria.pdf")
