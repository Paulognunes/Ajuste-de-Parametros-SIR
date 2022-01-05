def gera_amostra(data_inicio_corte, data_fim_corte, arquivo_confirmados, arquivo_recuperados, arquivo_obitos,
                 arquivo_datas):
    confirmados = open(arquivo_confirmados, 'r')
    recuperados = open(arquivo_recuperados, 'r')
    obitos = open(arquivo_obitos, 'r')
    data = open(arquivo_datas, 'r')

    n = 90497

    c = []
    r = []
    o = []
    d = []
    rt = []
    infec = []
    qnt_infec = 205
    d1 = []
    # Novos Casos Confirmados
    for i in confirmados:
        c.append(i.replace('\n', ''))

    # Número de casos confirmados acumulados
    for i in c:
        qnt_infec += int(i)
        infec.append(qnt_infec)

    # Número de recuperados acumulados
    for i in recuperados:
        r.append(i.replace('\n', ''))

    # Número de óbitos
    for i in obitos:
        o.append(i.replace('\n', ''))

    # Vetor de datas
    for i in data:
        aux = i.replace('/', '_')
        d.append(aux.replace('\n', ''))
        d1.append(i)
    # Desde 01/07/2020
    # Soma os óbitos aos recuperados
    x = 5
    for i in range(len(r)):
        if int(o[i]) == 0:
            aux = x + int(r[i])
            rt.append(aux)
        else:
            x += int(o[i])
            aux = x + int(r[i])
            rt.append(aux)

    data_inicio_corte = data_inicio_corte.replace('/', '_')
    data_fim_corte = data_fim_corte.replace('/', '_')

    corte = 0
    tam_amostra = 0

    for i in range(len(d)):
        if d[i] == data_inicio_corte:
            corte = i
            print(corte)
            print(c[corte])
            print(infec[corte])
            print(rt[corte])

    for i in range(len(d)):
        if d[i] == data_fim_corte:
            tam_amostra = i - corte

    if corte > (442 - tam_amostra):
        print("Data ou tamanho da amostra inválidos")
        return

    a = 'dados_sjdr_' + str(d[corte]) + '__' + str(d[corte + tam_amostra]) + '.csv'
    dados = open(a, 'w')

    for i in range(corte, corte + tam_amostra):
        suc = n - int(infec[i])
        if i != corte + tam_amostra - 1:
            inf_ativos = int(infec[i]) - int(rt[i])
            dados.write(str(suc) + ',' + str(inf_ativos) + ',' + str(rt[i]) + '\n')
        else:
            inf_ativos = int(infec[i]) - int(rt[i])
            dados.write(str(suc) + ',' + str(inf_ativos) + ',' + str(rt[i]))

    for i in range(len(infec)):
        inf_ativos = int(infec[i]) - int(rt[i])
        print(f' {inf_ativos}, {d1[i]}\n')

    confirmados.close()
    recuperados.close()
    obitos.close()
    data.close()
    dados.close()


if __name__ == '__main__':
    gera_amostra('15/07/2021', '15/09/2021', 'confirmados.txt', 'recuperados.txt', 'obitos.txt', 'datas.txt')
