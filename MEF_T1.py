import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D as line
from scipy import linalg
from scipy.optimize import fsolve

# 1) Entrada de Dados:

a = 4.0
b = 3.4

coord = np.array([[0.0, 0.0],
                  [a, 0.0],
                  [a, b],
                  [2*a, b],
                  [2*a, b/2],
                  [2*a, 0.0],
                  [2*a+(a/4), b]]) # Coordenadas dos nos (x, y)

nnos = coord.shape[0]

mat = np.array([[2e8, 5e-3, 1.28e-4]]) # Propriedades dos materiais

inci = np.array([[0, 2, 0],
                 [2, 1, 0],
                 [2, 3, 0],
                 [3, 4, 0],
                 [4, 5, 0],
                 [3, 6, 0]]) # Incidencia dos elementos

nelem = inci.shape[0]

fig, ax = plt.subplots()

xx = coord[:, 0]
yy = coord[:, 1]

# 2) Plotagem do portico

for e in inci:
  ni = e[0]
  nf = e[1]
  xi = xx[ni]
  xf = xx[nf]
  yi = yy[ni]
  yf = yy[nf]
  elemento = line((xi, xf), (yi, yf), color = 'black')
  ax.add_line(elemento)

ax.scatter(xx, yy, color = 'blue', zorder = 2)
plt.ylabel('y (m)')
plt.xlabel('x (m)')
plt.grid(axis = 'both', ls = ':', color = 'grey')
plt.axis('scaled')
plt.xlim(-1, 10)
plt.ylim(-1, 4)

plt.savefig('portico.png')

# 3) Matriz de rigidez e vetor de forças

tam = 3*nnos  # Tres GL por no

Kg = np.zeros((tam, tam))
Fg = np.zeros((tam, 1))

q = -14 # kN/m - carga aplicada

for i in range(nelem):

  ni = inci[i, 0]
  nf = inci[i, 1]
  mate = inci[i, 2]

  xi = xx[ni]
  xf = xx[nf]
  yi = yy[ni]
  yf = yy[nf]

  Le = np.sqrt(((xf - xi)**2) + ((yf - yi)**2))
  E = mat[mate, 0]
  A = mat[mate, 1]
  I = mat[mate, 2]

  val1 = E*A/Le
  val2 = 12*E*I/(Le**3)
  val3 = 6*E*I/(Le**2)
  val4 = 4*E*I/Le
  val5 = 2*E*I/Le

  Kel = np.array([[val1, 0, 0, -val1,    0,      0],
                  [0, val2, val3, 0, -val2,   val3],
                  [0, val3, val4, 0, -val3,   val5],
                  [-val1, 0, 0, val1,    0,      0],
                  [0, -val2, -val3, 0, val2, -val3],
                  [0, val3, val5, 0, -val3,   val4]])

  # Transformacao para sistema de coord. global:

  cos = (xf-xi)/Le
  sen = (yf-yi)/Le
  T = np.array([[cos, sen, 0, 0, 0, 0],
                [-sen, cos, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, cos, sen, 0],
                [0, 0, 0, -sen, cos, 0],
                [0, 0, 0, 0, 0, 1]])

  Keg = np.dot(np.transpose(T), Kel)
  Keg = np.dot(Keg, T)

  pos = [3*ni, 3*ni+1, 3*ni+2, 3*nf, 3*nf+1, 3*nf+2]

  # Construcao da matriz de rigidez global:

  for j in range (Keg.shape[0]):
    for k in range (Keg.shape[1]):
      Kg[pos[j], pos[k]] += Keg[j, k]

  # Vetor de forcas elementar

  if i == 1: # Elementos com carga constante
    fel = (Le*q/2)*np.array([[0],
                             [1],
                             [Le/6],
                             [0],
                             [1],
                             [-Le/6]])

  elif i == 2: # Elementos com carga linear
    fel = (Le*q/4)*np.array([[0],
                             [3/5],
                             [2*Le/15],
                             [0],
                             [7/5],
                             [-Le/5]])

  else: # Sem carga
    fel = np.zeros((6, 1))

  # Transformacao para sistema de coord. global:

  feg = np.dot(np.transpose(T), fel)

  # Construcao do vetor de forcas global:

  for j in range (feg.shape[0]):
    Fg[pos[j]] += feg[j]

# 4) Implementacao de forças nodais

p = 24 # kN
h = 26 # kN

Fg[7,0] += -p  # kN
Fg[12,0] += h  # kN
Fg[19,0] += -p # kN

# 5) Condicoes de contorno essenciais e ordenacao do sistema

ncon = 7 # no. de condicoes
ordem = [2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 0, 1, 3, 4, 5, 15, 16] # Ordem de organizacao da matriz

K = np.zeros((tam, tam))
F = np.zeros((tam, 1))

for i in range(K.shape[0]):
  F[i, 0] = Fg[ordem[i], 0]
  for j in range(K.shape[1]):
    K[i, j] = Kg[ordem[i], ordem[j]]

# 6) Calculo dos deslocamentos

dim = tam - ncon

dl = linalg.solve(K[:dim, :dim], F[:dim, 0])
desloc = np.zeros((tam, 1)) # valores prescritos sao nulos
for i in range(dim):
  desloc[ordem[i]] = dl[i]

print('Deslocamentos: \n', desloc, '\n')

dx = []
dy = []
dr = []

for i in range(0, len(desloc), 3):
  dx.append(desloc[i,0])
  dy.append(desloc[i+1,0])
  dr.append(desloc[i+2,0])

print('Deslocamentos horizontais: ', dx, '\n')
print('Deslocamentos verticais: ', dy, '\n')
print('Rotacoes: ', dr, '\n')

# 7) Calculo das reacoes de apoio

fp = np.dot(K[dim:, :dim], dl)

fp = fp - F[dim:, 0]
print('Reacoes de apoio: \n', fp, '\n')

# 8) Funcoes de forma e suas derivadas

#Funcoes de forma:

def phi1(x):
  f = (1/2)*(1-x)
  return f

def phi2(x):
  f = (1/2) - ((3/4)*x) + ((1/4)*(x**3))
  return f

def phi3(x, Le):
  f = (Le/8)*(1 - x - (x**2) + (x**3))
  return f

def phi4(x):
  f = (1/2)*(1+x)
  return f

def phi5(x):
  f = (1/2) + ((3/4)*x) - ((1/4)*(x**3))
  return f

def phi6(x, Le):
  f = (Le/8)*(- 1 - x + (x**2) + (x**3))
  return f

# Derivadas Segundas

def dphi1(x):
  f = (3/2)*x
  return f

def dphi2(x, Le):
  f = (Le/8)*((6*x) - 2)
  return f

def dphi3(x):
  f = (-3/2)*x
  return f

def dphi4(x, Le):
  f = (Le/8)*((6*x) + 2)
  return f

# Derivadas Terceiras

def ddphi1(x):
  f = 3/2
  return f

def ddphi2(x, Le):
  f = (Le/8)*6
  return f

def ddphi3(x):
  f = -3/2
  return f

def ddphi4(x, Le):
  f = (Le/8)*6
  return f

# 9) Esforcos normais, momentos fletores, e esforcos cortantes

disp, ax_disp = plt.subplots(3, 2)

mom, ax_mom = plt.subplots(3, 2)
cort, ax_cort = plt.subplots(3, 2)

N = []
M = []
Q = []

# Fatores de deformacao para os diagramas

fator_def = 1700
f_mom = 1/20
f_cort = 1/20

for i in range(nelem):

  ni = inci[i, 0]
  nf = inci[i, 1]
  mate = inci[i, 2]

  xi = xx[ni]
  xf = xx[nf]
  yi = yy[ni]
  yf = yy[nf]

  Le = np.sqrt(((xf - xi)**2) + ((yf - yi)**2))
  E = mat[mate, 0]
  A = mat[mate, 1]
  I = mat[mate, 2]

  csis = np.linspace(-1, 1, 100)
  deform = np.zeros(len(csis))
  momento = np.zeros(len(csis))
  cortante = np.zeros(len(csis))

  cos = (xf-xi)/Le
  sen = (yf-yi)/Le

  # Re-calculo da matriz de transformacao

  T = np.array([[cos, sen, 0, 0, 0, 0],
                [-sen, cos, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, cos, sen, 0],
                [0, 0, 0, -sen, cos, 0],
                [0, 0, 0, 0, 0, 1]])

  pos = [3*ni, 3*ni+1, 3*ni+2, 3*nf, 3*nf+1, 3*nf+2]

  Ue = np.array([[desloc[pos[0], 0]],
                [desloc[pos[1], 0]],
                [desloc[pos[2], 0]],
                [desloc[pos[3], 0]],
                [desloc[pos[4], 0]],
                [desloc[pos[5], 0]]])

  ue = np.dot(T, Ue) # Vetor de deslocamentos no sist. local

  # Separacao dos deslocamentos em barra e viga:

  barra = np.array([ue[0,0],
                    ue[3,0]])

  viga = np.array([ue[1,0],
                   ue[2,0],
                   ue[4,0],
                   ue[5,0]])

  for j in range(len(csis)):

    # Calculo da posicao deformada ponto a ponto:
    deform[j] = phi1(csis[j])*ue[0,0]*fator_def + phi2(csis[j])*ue[1,0]*fator_def +  \
    phi3(csis[j], Le)*ue[2,0]*fator_def + phi4(csis[j])*ue[3,0]*fator_def + \
    phi5(csis[j])*ue[4,0]*fator_def + phi6(csis[j], Le)*ue[5,0]*fator_def

    # Calculo do momento fletor ponto a ponto:
    momento[j] = (4*E*I/(Le**2))*(dphi1(csis[j])*viga[0] + dphi2(csis[j], Le)*viga[1] \
                                  + dphi3(csis[j])*viga[2] + dphi4(csis[j], Le)*viga[3])

    # Calculo do esforco cortante ponto a ponto:
    cortante[j] = (8*E*I/(Le**3))*(ddphi1(csis[j])*viga[0] + ddphi2(csis[j], Le)*viga[1]\
                                   + ddphi3(csis[j])*viga[2] + ddphi4(csis[j], Le)*viga[3])


  # Calculo de esforcos normais

  normal = (E*A/Le)*(barra[1] - barra[0])

  N.append(normal)
  M.append((momento[0], momento[-1]))
  Q.append((cortante[0], cortante[-1]))

  # 10) Plotagem elemento a elemento

  # Elementos em cada grafico:

  elemento1 = line((-1, 1), (0, 0), color = 'black')
  elemento2 = line((-1, 1), (0, 0), color = 'black')
  elemento3 = line((-1, 1), (0, 0), color = 'black')

  title = 'Elemento ' + str(i)

  # Plotagem dos valores obtidos

  if i <= 2:

    ax_disp[i, 0].add_line(elemento1)
    ax_disp[i, 0].set_title(title)
    ax_disp[i, 0].plot(csis, deform, ls = '--', color = 'red')

    ax_mom[i, 0].add_line(elemento2)
    ax_mom[i, 0].set_title(title)
    ax_mom[i, 0].plot(csis, momento, ls = '--', color = 'red')

    ax_cort[i, 0].add_line(elemento3)
    ax_cort[i, 0].set_title(title)
    ax_cort[i, 0].plot(csis, cortante, ls = '--', color = 'red')

  else:

    ax_disp[i-3, 1].add_line(elemento1)
    ax_disp[i-3, 1].set_title(title)
    ax_disp[i-3, 1].plot(csis, deform, ls = '--', color = 'red')

    ax_mom[i-3, 1].add_line(elemento2)
    ax_mom[i-3, 1].set_title(title)
    ax_mom[i-3, 1].plot(csis, momento, ls = '--', color = 'red')

    ax_cort[i-3, 1].add_line(elemento3)
    ax_cort[i-3, 1].set_title(title)
    ax_cort[i-3, 1].plot(csis, cortante, ls = '--', color = 'red')

print('\n Esforços normais: ', N, '\n')
print('Momentos fletores (nó incial, nó final): ', M, '\n')
print('Esforços cortantes: (nó incial, nó final)', Q, '\n')

# Parametros graficos para os plots elemento a elemento

for (m,n), subplot in np.ndenumerate(ax_disp):
    subplot.set_xlim(-2, 2)
    subplot.set_ylim(-4, 4)
    subplot.grid(axis = 'both', linestyle = ':')

disp.text(-0.01, 0.5, 'Deslocamento', ha='center', va='center', rotation='vertical')
disp.tight_layout()

disp.savefig('deslocamento.png')

for (m,n), subplot in np.ndenumerate(ax_mom):
    subplot.set_xlim(-2, 2)
    subplot.set_ylim(-30, 30)
    subplot.grid(axis = 'both', linestyle = ':')

mom.text(-0.01, 0.5, 'Momento Fletor (kNm)', ha='center', va='center', rotation='vertical')
mom.tight_layout()

mom.savefig('fletor.png')

for (m,n), subplot in np.ndenumerate(ax_cort):
    subplot.set_xlim(-2, 2)
    subplot.set_ylim(-30, 30)
    subplot.grid(axis = 'both', linestyle = ':')

cort.text(-0.01, 0.5, 'Esforço Cortante (kN)', ha='center', va='center', rotation='vertical')
cort.tight_layout()

cort.savefig('cortante.png')

# 11) Plotagem dos diagramas

figDisp, axDisp = plt.subplots()
figMom, axMom = plt.subplots()
figCort, axCort = plt.subplots()

for i in range(nelem):

  ni = inci[i, 0]
  nf = inci[i, 1]
  mate = inci[i, 2]

  xi = xx[ni]
  xf = xx[nf]
  yi = yy[ni]
  yf = yy[nf]

  Le = np.sqrt(((xf - xi)**2) + ((yf - yi)**2))
  E = mat[mate, 0]
  A = mat[mate, 1]
  I = mat[mate, 2]

  csis = np.linspace(-1, 1, 100) # Vetor em xi
  xis = np.linspace(0, Le, len(csis)) # Vetor em x

  deform = np.zeros(len(csis))
  momento = np.zeros(len(csis))
  cortante = np.zeros(len(csis))

  cos = (xf-xi)/Le
  sen = (yf-yi)/Le
  theta = np.arcsin(sen) # Angulo do elemento

  # Re-calculo da matriz de tranformacao

  T = np.array([[cos, sen, 0, 0, 0, 0],
                [-sen, cos, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, cos, sen, 0],
                [0, 0, 0, -sen, cos, 0],
                [0, 0, 0, 0, 0, 1]])

  pos = [3*ni, 3*ni+1, 3*ni+2, 3*nf, 3*nf+1, 3*nf+2]

  Ue = np.array([[desloc[pos[0], 0]],
                [desloc[pos[1], 0]],
                [desloc[pos[2], 0]],
                [desloc[pos[3], 0]],
                [desloc[pos[4], 0]],
                [desloc[pos[5], 0]]])

  ue = np.dot(T, Ue) # Vetor de deslocamentos no sist. local

  # Separacao dos deslocamentos em barra e viga:

  barra = np.array([ue[0,0],
                    ue[3,0]])

  viga = np.array([ue[1,0],
                   ue[2,0],
                   ue[4,0],
                   ue[5,0]])

  for j in range(len(csis)):

    # Calculo da posicao deformada ponto a ponto:
    deform[j] = phi1(csis[j])*ue[0,0]*fator_def + phi2(csis[j])*ue[1,0]*fator_def +  \
    phi3(csis[j], Le)*ue[2,0]*fator_def + phi4(csis[j])*ue[3,0]*fator_def + \
    phi5(csis[j])*ue[4,0]*fator_def + phi6(csis[j], Le)*ue[5,0]*fator_def

    # Calculo do momento fletor ponto a ponto:
    momento[j] = (-4*E*I/(Le**2))*(dphi1(csis[j])*viga[0] + dphi2(csis[j], Le)*viga[1]\
                                   + dphi3(csis[j])*viga[2] + dphi4(csis[j], Le)*viga[3])*f_mom

    # Calculo da esforco cortante ponto a ponto:
    cortante[j] = (-8*E*I/(Le**3))*(ddphi1(csis[j])*viga[0] + ddphi2(csis[j], Le)*viga[1]\
                                    + ddphi3(csis[j])*viga[2] + ddphi4(csis[j], Le)*viga[3])*f_cort

  xdata = np.array([xi, xf]) # Coord. x dos nós
  ydata = np.array([yi, yf]) # Coord. y dos nós

  # PLotgem dos elementos em cada diagrama:

  axDisp.plot (xdata, ydata, color = 'black')
  axMom.plot (xdata, ydata, color = 'black')
  axCort.plot (xdata, ydata, color = 'black')

  dx = xi # Translacao horizontal a ser aplicada para determinado elemento
  dy = yi # Translacao vertical a ser aplicada para determinado elemento

  # Matriz de rotacao para o angulo do elemento:

  rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                              [np.sin(theta), np.cos(theta)]])

  # Rotacao da funcao:
  rotated_disp = np.dot(rotation_matrix, np.vstack((xis, deform)))
  # Translacao da funcao:
  translated_disp = rotated_disp + np.array([[dx], [dy]])

  dispX_trans, dispY_trans = translated_disp # Deslocamentos transformados

  # Plotagem dos dados transformados em cada diagrama

  rotated_mom = np.dot(rotation_matrix, np.vstack((xis, momento)))
  translated_mom = rotated_mom + np.array([[dx], [dy]])

  momX_trans, momY_trans = translated_mom # Momentos transformados

  rotated_cort = np.dot(rotation_matrix, np.vstack((xis, cortante)))
  translated_cort = rotated_cort + np.array([[dx], [dy]])

  cortX_trans, cortY_trans = translated_cort # Cortantes transformados

  axDisp.plot(dispX_trans, dispY_trans, color = 'red', ls = '--')

  axMom.plot(momX_trans, momY_trans, color = 'red', ls = '--')

  axCort.plot(cortX_trans, cortY_trans, color = 'red', ls = '--')

axes = [axDisp, axCort, axMom] # Lista para iterar sobre os diagramas

# Parametros graficos dos diagramas

for ax in axes:

  ax.scatter(xx, yy, color = 'blue', zorder = 2)
  ax.set_ylabel('y (m)')
  ax.set_xlabel('x (m)')
  ax.grid(axis = 'both', ls = ':', color = 'grey')
  ax.axis('scaled')

  ax.set_xlim(-1, 10)
  ax.set_ylim(-1,6)

figDisp.savefig('D.png')
figMom.savefig('M.png')
figCort.savefig('Q.png')
