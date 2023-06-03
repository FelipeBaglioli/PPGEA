# Funcoes:

def m0t (cc, tt):
  '''
  Calcula o momento temporal de ordem zero de uma serie C(t) atraves do 
  metodo de integracao discreta dos trapezios:

  ENTRADAS:
  - cc (array): serie temporal de concentracoes;
  - tt (array): instantes de tempo relacionados as concentracoes;

  SAIDAS:
  - AA (float64): valor do momento temporal de ordem zero.
  '''
  if len(cc) != len(tt):
    print ("ERRO: Vetor de tempo e concentracao devem ter o mesmo tamanho.")
    return None
  else:
    N = len(cc)
  
  AA = 0

  for i in range(1, N-1):

    fnew = cc[i]
    fold = cc[i-1]
    da = (fnew + fold)*(tt[i] - tt[i-1])/2
    AA += da

  return AA

def t_med (cc, tt):
  '''
  Calcula o tempo medio de uma serie C(t) atraves do metodo de integracao 
  discreta dos trapezios:

  ENTRADAS:
  - cc (array): serie temporal de concentracoes;
  - tt (array): instantes de tempo relacionados as concentracoes;

  SAIDAS:
  - AA (float64): valor do tempo medio.
  '''
  if len(cc) != len(tt):
    print ("ERRO: Vetor de tempo e concentracao devem ter o mesmo tamanho.")
    return None
  else:
    N = len(cc)
  
  AA = 0

  for i in range(1, N-1):

    fnew = cc[i]*tt[i]
    fold = cc[i-1]*tt[i-1]
    da = (fnew + fold)*(tt[i] - tt[i-1])/2
    AA += da

  return AA

def sigmaT (cc, tt, t_med):
  '''
  Calcula a variancia temporal de uma serie C(t) atraves do metodo de integracao
  discreta dos trapezios:

  ENTRADAS:
  - cc (array): serie temporal de concentracoes;
  - tt (array): instantes de tempo relacionados as concentracoes;
  - t_med (float64): tempo medio da serie.
  
  SAIDAS:
  - AA (float64): valor da variancia temporal.
  '''
  if len(cc) != len(tt):
    print ("ERRO: Vetor de tempo e concentracao devem ter o mesmo tamanho.")
    return None
  else:
    N = len(cc)
  
  AA = 0

  for i in range(1, N-1):

    fnew = cc[i]*((tt[i] - t_med)**2)
    fold = cc[i-1]*((tt[i-1] - t_med)**2)
    da = (fnew + fold)*(tt[i] - tt[i-1])/2
    AA += da

  return AA

def propaga_C(cc, tt, t1, t2, tj, u, D):
  '''
  Propaga uma distribuicao inicial de concentracoes de um ponto inicial (1) até
  um ponto a jusante (2), em determinado instante de tempo.

  ENTRADAS:
  - cc (array): serie temporal de concentracoes no ponto inicial;
  - tt (array): instantes de tempo relacionados as concentracoes inicias;
  - t1 (float64): tempo medio de observacao da pluma em (1);
  - t2 (float64): tempo medio de observacao da pluma em (2);
  - tj (float64): tempo em que queremos observar a concentracao propagada;
  - u (float64): velocidade media longitudinal do escoamento;
  - D (float64): coeficiente de dispersao longitudinal do escoamento
  
  SAIDAS:
  - AA (float64): valor concentracao propagada C(x2, tj).
  '''
  if len(cc) != len(tt):
    print ("ERRO: Vetor de tempo e concentracao devem ter o mesmo tamanho.")
    return None
  else:
    N = len(cc)
  
  AA = 0
  ts = t2 - t1 -tj
  deltaT = t2 - t1
  a = u/(np.sqrt(4*D*np.pi*deltaT))


  for i in range(1, N-1):

    bnew = ((u*(ts + tt[i]))**2)/(4*D*deltaT)
    bold = ((u*(ts + tt[i-1]))**2)/(4*D*deltaT)

    fnew = a*cc[i]*np.exp(-bnew)
    fold = a*cc[i-1]*np.exp(-bold)
    da = (fnew + fold)*(tt[i] - tt[i-1])/2
    AA += da

  return AA


def erro(D):
  '''
  Calcula o erro quadratico medio entre uma distribuicao medida
  de concentracoes e uma distribuicao modelada atraves de
  propaga_C(cc, tt, t1, t2, tj, u, D)

  ENTRADAS:
  - Entrada indireta: parametros de propaga_ C(cc, tt, t1, t2, tj, u).
  - D (float64): coeficiente de dispersao longitudinal da pluma
  propagada.
  
  SAIDAS:
  - AA (float64): valor do erro.
  '''

  c_mod = []
  t_mod = []
  err = 0

  for t in tt2:
    ct = propaga_C(cc1*1000000, tt1, t_med1, t_med2, t, v, D)
    c_mod.append(ct)
    t_mod.append(t)

  return mean_squared_error(cc2*1000000, c_mod)


# Importacao de bibliotecas:

import numpy as np
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
mp.rcParams.update(mp.rcParamsDefault) # Resetar os parametros de plotagem

# Leitura do arquivo
file_in = "/content/campo6-23-05-23.xlsx"

data_in = pd.read_excel(file_in, header = 1)

# Cria uma copia do arquivo original para trabalharmos sobre:

df = data_in.copy()

# Separa as informacoes de cada ponto:

x1 = df[['t (s)', 'C (ppb)']]
x2 = df[['t (s).1', 'C (ppb).1']]
x2 = x2.rename(columns = {'t (s).1': 't (s)','C (ppb).1' : 'C (ppb)'})

# Retira linhas com dados faltantes:

x1 = x1.dropna()
x2 = x2.dropna()

# Calculo da concentração de background:

m1 = x1['C (ppb)'].min()
m2 = x2['C (ppb)'].min()

bg = min(m1, m2)
print("Concentracao de background: ", bg, ' ppb.\n')

# Transformando os dataframes em vetores, e retirando o bg das concentracoes:

cc1 = np.array(x1['C (ppb)'])-bg
tt1 = np.array(x1['t (s)'])

cc2 = np.array(x2['C (ppb)'])-bg
tt2 = np.array(x2['t (s)'])

# Plot inicial dos dados medidos (Sem retirar o bg)

fig = plt.figure(figsize = (12, 6))

plt.plot(x1['t (s)']/60, x1['C (ppb)'], color = '#FF6347', label = 'x1')
plt.plot(x2['t (s)']/60, x2['C (ppb)'], color = 'purple', label = 'x2')
plt.axhline(bg, ls = '--', color = 'black', label = 'Background')

plt.ylim(0, 35)
plt.ylabel('C (ppb)')

plt.xlim(0, 200)
plt.xlabel('t (min)')

plt.legend(loc = 'best')
plt.grid(axis = 'both', ls = ':')
plt.savefig('dados.png')  # Salva a figura como dados.png

# Metodo dos Momentos

# Retira os valores de background:

cc1 = np.array(x1['C (ppb)'])-bg
tt1 = np.array(x1['t (s)'])

cc2 = np.array(x2['C (ppb)'])-bg
tt2 = np.array(x2['t (s)'])

# Dados extras de entrada:

M = 50 # g
d1 = 1047 # m
d2 = 2292 # m
 
deltax = d2 - d1

cc1 = cc1/1000000 # Conversao para kg/m³

# x1:

m01 = m0t(cc1, tt1)

Q1 = (M/1000)/m01

t_med1 = t_med(cc1, tt1)/m01
sigma1 = sigmaT(cc1, tt1, t_med1)/m01

# x2:

cc2 = cc2/1000000

m02 = m0t(cc2, tt2)

Q2 = (M/1000)/m02

t_med2 = t_med(cc2, tt2)/m02
sigma2 = sigmaT(cc2, tt2, t_med2)/m02

# Resultados:

Q = (Q1+Q2)/2

print('Vazao media estimada: ', '{:.5f}'.format(Q), 'm³/s.\n')

massa1 = Q*m01
massa2 = Q*m02

print('Em x1:\n ')

print("Estatisticas: ")
print('Massa = ', '{:.5f}'.format(massa1*1000), ' g')
print('Tempo medio = ', '{:.3f}'.format(t_med1), ' s')
print('Variancia temporal = ', '{:.3f}'.format(sigma1), ' s²\n')

print('Em x2:\n ')

print("Estatisticas: ")
print('Massa = ', '{:.5f}'.format(massa2*1000), ' g')
print('Tempo medio = ', '{:.3f}'.format(t_med2), ' s')
print('Variancia temporal = ', '{:.3f}'.format(sigma2), ' s²\n')

v = (deltax)/(t_med2 - t_med1)

print('Velocidade media estimada: ', '{:.5f}'.format(v), 'm/s.\n')

A = Q/v

print('Secao media estimada: ', '{:.5f}'.format(A), 'm².\n')

D = ((v**2)/2)*((sigma2 - sigma1)/(t_med2 - t_med1))

print('Coeficiente de dispersao longitudinal estimado (momentos) : ', '{:.3f}'.format(D), ' m²/s\n')

print('EQM : ', '{:.3f}'.format(erro(D)), ' (kg/m³)².\n')

plt.rcParams.update({'axes.titlesize': 'medium'})

# Calcula um distribuicao propagada usando D0

c_mod = []
t_mod = []

for t in tt2:
  ct = propaga_C(cc1*1000000, tt1, t_med1, t_med2, t, v, D)
  c_mod.append(ct)
  t_mod.append(t)

# Plota uma comparacao entre os dados medidos e modelados

plt.plot(np.array(t_mod)/60, c_mod, color = '#FF6347', label = 'Pluma propagada')
plt.plot(tt2/60, cc2*1000000, color = 'purple', label = 'Pluma medida')
plt.title('D = ' + '{:.3f}'.format(D) + ' m²/s')

plt.grid(axis = 'both', ls = ':')

plt.ylim(0, 16)
plt.ylabel('C (ppb)')

plt.xlim(40, 180)
plt.xlabel('t (min)')

plt.legend(loc = 'best')
plt.savefig('disp_momento.png')

# Define o chute incial como D0, e otimiza usando gradiente conjugado

D0 = D
Opt = minimize(erro, D0, tol = 0.0000001, method = 'CG')

print('Resultado da otimizacao: \n')
print(Opt)
print('\n')

# Salva o valor de D otimizado

Dnew = Opt.x[0]

print('Coeficiente de dispersao longitudinal estimado (propagacao): ', '{:.3f}'.format(Dnew), ' m²/s')

print('EQM : ', '{:.3f}'.format(erro(Dnew)), ' (kg/m³)².\n')

# Calcula uma nova distribuicao propagada, agora com Dnew

c_mod = []
t_mod = []

for t in tt2:
  ct = propaga_C(cc1*1000000, tt1, t_med1, t_med2, t, v, Dnew)
  c_mod.append(ct)
  t_mod.append(t)

plt.plot(np.array(t_mod)/60, c_mod, color = '#FF6347', label = 'Pluma propagada')
plt.plot(tt2/60, cc2*1000000, color = 'purple', label = 'Pluma medida')
plt.title('D = ' + '{:.3f}'.format(Dnew) + ' m²/s')

plt.grid(axis = 'both', ls = ':')

plt.ylim(0, 16)
plt.ylabel('C (ppb)')

plt.xlim(40, 180)
plt.xlabel('t (min)')

plt.legend(loc = 'best')
plt.savefig('disp_propag.png')