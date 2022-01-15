import matplotlib.pyplot as plt
import sys
import numpy as np

arquivo = sys.argv[1] #nome do arquivo
variavel = sys.argv[2] #variavel a ser plotada como no arquivo .tab
pontos = int(sys.argv[3]) #pontos para plotar

tempo = []
resposta = []
count = 0
valores = []

with open(arquivo, 'r') as fp:
    for linha in fp.readlines():
        if count == 0:
            variaveis = linha.strip().split(' ')
        if count == 2:
            passo = float(linha.strip().split(' ')[0])
        count = count + 1
        if count > 1:
            valores.append(linha.strip().split(' '))
tempoIndice = variaveis.index('t')
variavelIndice = variaveis.index(variavel)
valores = np.array(valores, dtype=float)

print 'Variaveis: ' + str(variaveis)
print 'Plotando variavel %s'%variavel
print 'Passo %f'%passo

plt.figure(figsize=(12,6))
plt.plot(valores[:pontos,tempoIndice], valores[:pontos,variavelIndice], color='b')
axes = plt.gca()
axes.set_xlim([0,np.max(valores[:pontos,tempoIndice])])
axes.set_ylim([np.min(valores[:pontos,variavelIndice]),np.max(valores[:pontos,variavelIndice])])
plt.title('%s x Tempo'%variavel)
plt.ylabel('Amplitude')
plt.xlabel('Tempo (s)')

plt.savefig('variavel(%s)xtempo.png'%variavel)
