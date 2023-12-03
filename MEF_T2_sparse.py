# -*- coding: utf-8 -*-
'''
----------------------------------------------------------------------
MEF_T2_sparse.py : resolve um viga bidimensional utilizando MEF e a biblioteca
de scipy para matrizes esparsas.

Autor: Felipe Baglioli
Data de criacao: 04/09/2023
Ultima alteracao em: 16/04/2023 - Obtencao de tensoes horizontais
Contexto: MNUM7122 - 2o trim. de 2023 - Trabalho 2
----------------------------------------------------------------------
ENTRADAS:

- L (float): comprimento em m;
- h (float): altura em m;
- t (float): espessura em m;
- E (float): modulo de elasticidade em kN/m2;
- nu (float): coef. de Poisson;
- qx e qy (float): cargas de superficie em kN;
- px e py (float): cargas de volume em kN;
- nx e ny (int): numero de nos da malha em cada direcao.

----------------------------------------------------------------------
SAIDAS:

- Graficos de cada uma das malhas utilizadas;
- Grraficos da configuracao deformada da viga para cada malha;
- Diagrama de convergencia (deslocamento x no. de nos);
- Diagrama de tensoes na secao transversal central;

----------------------------------------------------------------------
'''
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import scipy
from scipy.spatial import Delaunay
from scipy.sparse.linalg import spsolve

# Parametros gerais para os plots de malha e deformacao:

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

rc('font', **font)

# ---------------------------- FUNCOES ----------------------------

def calc_lenght(x1, y1, x2, y2):
  '''
  --------------------------------------------------------------------
  Calcula o comprimento de um segmento de reta (x1, y1) ate (x2, y2).
  --------------------------------------------------------------------
  ENTRADAS:

  - x1, y1, x2, y2 (int): coordenadas dos pontos que definem a reta.
  --------------------------------------------------------------------
  SAIDAS:

  - side (integer): comprimento do segmento.
  --------------------------------------------------------------------
  '''

  side = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

  return side

# --------------------------------------------------------------------

def delaunay_mesh_mid(L, h, divX, divY, save = False, node_dots = False,
                      node_labels = False, element_labels = False, 
                      target_coords = None):
  '''
  --------------------------------------------------------------------
  Gera uma malha triangular para um dominio retangular, utilizando o 
  algoritmo de triangulacao de Delaunay e garantindo a colocacao
  de um no na posicao central inferior do dominio.
  --------------------------------------------------------------------
  ENTRADAS:

  - L, h (float): comprimento e altura do dominio, respectivamente;
  - divX e divY (int): numero de nos a serem posicionados em cada 
    direcao;
  - target_coords (tupla de 2 val.): contem coordenadas de um no 
    cujo indice deseja ser conhecido.

  FLAGS (boolean):

  - save: padrao False. Caso seja True, a figura da malha sera salva;
  - node_dots: padrao False. Caso seja True, plota os pontos dos nos;
  - node_label: padrao False. Caso seja True, plota os numeros dos nos;
  - element_labels: padrao False. Caso seja True, plota os
    numeros dos elementos;
  --------------------------------------------------------------------
  SAIDAS:

  - nodes (float, array 2D): coordenadas dos nos da malha;
  - inci_sorted (int, array 2D): matriz de incidencia dos elementos;
  - target_index (int): indice do no com as coordenadas informadas.

  - Caso save = True: grafico da malha gerada.
  --------------------------------------------------------------------
  '''

  # Garante que o ponto central inferior seja um no:

  x_points = np.unique(np.concatenate([np.linspace(0, L, divX), [L/2]]))
  y_points = np.unique(np.concatenate([np.linspace(0, h, divY), [0]]))

  # vetores com coordenadas dos pontos:

  X, Y = np.meshgrid(x_points, y_points)

  # conversao em um array 2d:

  vertices = np.column_stack((X.ravel(), Y.ravel()))

  # Criacao da malha:

  tri = Delaunay(vertices)

  # Coleta das coord. dos nos e da matriz de incidencia:

  nodes = tri.points
  inci = tri.simplices

  # Calculos dos centroides de cada elemento:

  centroids = np.mean(nodes[inci], axis=1)

  # Enumeracao dos elementos de acordo com centroides, para melhor condicionamento:

  sorting_order = np.lexsort((centroids[:, 0], centroids[:, 1]))
  inci_sorted = inci[sorting_order]

  fig = plt.figure(figsize = (24, 16))

  # Plotagem da malha:

  plt.triplot(nodes[:, 0], nodes[:, 1], inci_sorted, color = 'green', 
  ls = '-', lw = 0.8)

  # Plotagem dos nos:
  if node_dots == True:
    plt.scatter(nodes[:, 0], nodes[:, 1], color = 'orange', 
    marker = 'o', s = 70)

  offsetx = 0.12
  offsety = 0.2

  # Anotar numeros dos elementos em seus centroides
  for element_number, element in enumerate(inci_sorted):

      centroid_x = np.mean(nodes[element, 0])
      centroid_y = np.mean(nodes[element, 1])

      if element_labels == True:
        plt.text(centroid_x, centroid_y, str(element_number), 
        ha='center', va='center', color='green', fontsize=20)

  # Anotar numeros dos nos:
  if node_labels == True:

    for node_number, node_coordinate in enumerate(nodes):
        plt.text(node_coordinate[0] + offsetx, 
        node_coordinate[1] + offsety, str(node_number), 
        ha='center', va='center', color='orange', fontsize=20)

  N = len(inci_sorted) # Numero de elementos

  # Parametros de plotagem:

  plt.xlabel('X (m)')
  plt.ylabel('Y (m)')
  plt.title('N = ' + str(N))
  plt.axis('scaled')
  plt.xlim(-1, L+1)
  plt.ylim(-1, h+1)
  plt.grid(axis = 'both', ls = ':')
  plt.tight_layout()
  if save == True:
    plt.savefig('grid'+str(N)+'.png')
  else:
    plt.show()

  # Obtencao do indice do no passado como target_index

  if target_coords is not None:
      for i, node in enumerate(nodes):
          if np.allclose(node, target_coords):
              target_index = i
              break

  # Retorno dos resultados de interesse:
  return nodes, inci_sorted, target_index

# -------------------------------------------------------------------

def deform_plot(nodes, desloc, els, save = False, fator_def = 1.0,
                node_dots = False, node_labels = False, 
                element_labels = False):
  '''
  --------------------------------------------------------------------
  Plota a configuracao deformada de uma certa malha de nos.
  --------------------------------------------------------------------
  ENTRADAS:

  - nodes (float, array 2D): array com as coord. inciais dos nos;
  - nodes (float, array 1D): array com os deslocamentos nodais, organizados
    verticalmente por no, com ordem x e depois y;
  - els (int, array 1D): matriz de incidencia dos elementos;
  - fator_def (float): padrao 1. Amplifica ou reduz a deformacao 
    real no plot.

  FLAGS (boolean):

  - save: padrao False. Caso seja True, a figura da malha sera salva;
  - node_dots: padrao False. Caso seja True, plota os pontos dos nos;
  - node_label: padrao False. Caso seja True, plota os numeros dos nos;
  - element_labels: padrao False. Caso seja True, plota os
                    numeros dos elementos;
  --------------------------------------------------------------------
  SAIDAS:

  - Caso save = True: grafico da malha gerada.
  --------------------------------------------------------------------
  '''

  # Reorganizar a matriz de deslocamentos para um array 2D e aplicar fator_def:

  deformed_node_coordinates = nodes + (desloc.reshape((-1, 2))*fator_def)

  # Plotar a grade deformada:
  fig = plt.figure(figsize = (24, 16))

  plt.triplot(deformed_node_coordinates[:, 0], 
              deformed_node_coordinates[:, 1], 
              els, color = 'green', ls = '-', lw = 0.8)

  # Adicionar pontos dos nos:
  if node_dots == True:
    plt.scatter(deformed_node_coordinates[:, 0], 
                deformed_node_coordinates[:, 1], 
                color = 'orange', marker = 'o', s = 70)

  offsetx = 0.1
  offsety = 0.15

  # Adicionar numeros dos elementos:
  for element_number, element in enumerate(els):
      centroid_x = np.mean(deformed_node_coordinates[element, 0])
      centroid_y = np.mean(deformed_node_coordinates[element, 1])
      if element_labels == True:
        plt.text(centroid_x, centroid_y, str(element_number), 
        ha='center', va='center', color='green', fontsize=20)

  # Adicionar numeros dos nos:
  if node_labels == True:
    for node_number, node_coordinate in enumerate(deformed_node_coordinates):
      plt.text(node_coordinate[0] + offsetx, 
      node_coordinate[1] + offsety, str(node_number), 
      ha='center', va='center', color='orange', fontsize=20)

  # Parametros do plot

  N = len(els)

  plt.xlabel('X (m)')
  plt.ylabel('Y (m)')
  plt.title('N = ' + str(N))
  plt.axis('scaled')
  plt.xlim(-1, L+1)
  plt.ylim(-1, h+1)
  plt.grid(axis = 'both', ls = ':')
  plt.tight_layout()
  if save == True:
    plt.savefig('def'+str(N)+'.png')
  else:
    plt.show()

  return None

#--------------------- SOLVER -----------------------------

# Parametros de entrada:

L = 9.0 # m (comprimento)
h = 3.0 # m (altura)
t = 0.1 # m (espessura)

E = 2.0e+7 # kN/m2 (modulo de elasticidade)
nu = 0.3   # adim. (coef. de Poisson)

p = -25.0    # kN/m3 (Peso proprio)
q = - 200.0  # kN/m (Carga distribuida)

# Array com numero de divisoes em x e y para cada malha:

mesh_divs = np.array([[7, 3],
                      [7, 7],
                      [11, 11],
                      [21, 11],
                      [21, 21],
                      [31, 31]])

N_elmnt = []
defmid = []

for div in mesh_divs: # loop para varias malhas

  # 1) Malha, matriz D e nos de topo:

  divX = div[0]
  divY = div[1]

  nodes, els, mid_idx = delaunay_mesh_mid(L, h, divX, divY, 
  target_coords=(L/2, 0), save = True)

  nnos = len(nodes)
  nelem = len(els)

  N_elmnt.append(nelem)

  #print('No central inferior: ', mid_idx)

  D = E/(1-(nu**2))*np.array([[1, nu, 0],
                              [nu, 1, 0],
                              [0, 0, (1-nu)/2]])

  tam = nnos*2 # 2 DoF by node

  Kg = np.zeros(shape = (tam, tam))
  Fg = np.zeros(shape = (tam, 1))

  top_nodes = []

  for i in range(nnos):

    if nodes[i][1] == h:
      top_nodes.append(i)

  #print(top_nodes)

  for i, e in enumerate(els): # Loop em cada elemento:

    # 2) Calculo da matriz de rigidez elementar:

    #print('Elemento ', i)
    #print(e)

    x1, y1 = nodes[e[0]]
    x2, y2 = nodes[e[1]]
    x3, y3 = nodes[e[2]]

    #print('No 1: ', x1, y1, ' (No geral ', e[0], ')')
    #print('No 2: ', x2, y2, ' (No geral ', e[1], ')')
    #print('No 3: ', x3, y3, ' (No geral ', e[2], ')')

    # Area:

    Ae = 0.5*((x2*y3) + (x1*y2) + (x3*y1) - (x2*y1) - (x1*y3) - (x3*y2))

    #print('Area: ', Ae)

    b1 = y2 - y3
    b2 = -y1 + y3
    b3 = y1 - y2

    c1 = -x2 + x3
    c2 = x1 - x3
    c3 = -x1 + x2

    # Matriz deslocamento-deformacao:

    B = (1/(2*Ae))*np.array([[b1, 0, b2, 0, b3, 0],
                             [0, c1, 0, c2, 0, c3],
                             [c1, b1, c2, b2, c3, b3]])

    #print('\nMatriz B: ')
    #print(B)

    Ke = np.dot(D, B)
    Ke = np.dot(np.transpose(B), Ke)
    Ke = t*Ae*Ke

    #print('\nMatriz de rigidez: ')
    #print(Ke)

    # 3) Vetor de forcas de volume elementar

    px = 0
    py = p

    fe_vol = ((t*Ae)/3)*np.array([[px],
                                  [py],
                                  [px],
                                  [py],
                                  [px],
                                  [py]])

    #print('\nForcas de volume: ')
    #print(fe_vol)

    #print('\n')

    # 4)Vetor de forcas de superficie elementar

    qy = q  #  kN/m
    qx = 0

    fe_sup = np.zeros_like(fe_vol) 
    # Nulo em elementos que nao sao de topo

    count_top_nodes = 0
    el_top_nodes = []
    local_node_idx = []

    for j in range(3):

      if e[j] in (top_nodes):
        count_top_nodes += 1
        el_top_nodes.append(e[j])
        local_node_idx.append(j)

        # Elementos com 2 nos de topo tem aresta no topo:

        if count_top_nodes == 2:
          #print('------- Elemento de topo. ------- \n')
          #print('Nos de topo: ', el_top_nodes)

          xi = nodes[el_top_nodes[0]][0]
          yi = nodes[el_top_nodes[0]][1]

          xf = nodes[el_top_nodes[1]][0]
          yf = nodes[el_top_nodes[1]][1]

          # Calcula o comprimento da aresta de topo
          top_side_length = calc_lenght(xi, yi, xf, yf)

          #print('Comprimento aresta de topo: ', top_side_length, '\n')

          fe_vol[(2*local_node_idx[0])] += (top_side_length/2)*qx
          fe_vol[(2*local_node_idx[0])+1] += (top_side_length/2)*qy

          fe_vol[(2*local_node_idx[1])] += (top_side_length/2)*qx
          fe_vol[(2*local_node_idx[1])+1] += (top_side_length/2)*qy

    # 5) Construcao da matriz de rigidez e vetor de forcas globais:

    fe = fe_sup + fe_vol

    #print('Vetor de forcas: ')
    #print(fe, '\n')

    #G.L.'s relacionados a cada elemento:

    pos = [e[0]*2, (e[0]*2) + 1, e[1]*2, 
          (e[1]*2) + 1, e[2]*2, (e[2]*2) + 1]

    #print('GLs :', pos, '\n')

    local_size = len(fe)

    for m in range(local_size): # Loop nas linhas
      Fg[pos[m]] += fe[m]
      for n in range(local_size): # loop nas colunas
        Kg[pos[m], pos[n]] += Ke[m, n]

  #print(Kg, '\n')
  #print(Fg, '\n')

  left_supX = L/10
  right_supX = L - (L/10)

  #print(left_supX)
  #print(right_supX)

  # 6) Condicoes de contorno essenciais:

  sup_gls = []

  for i, coord in enumerate(nodes):
    if (coord[1] == 0 and coord[0] <= left_supX) or (coord[1] == 0 and coord[0] >= right_supX):

      # Apoios atuam no GL vertical:
      glv = (i*2) + 1

      sup_gls.append(glv)
      #print(i, ' e no apoiado.')

  #print('C.C.s em : ', sup_gls)

  ncon = len(sup_gls)

  ordem = np.arange(0, tam, 1).tolist()

  #print(ordem)

  # 7) Ordenacao e resolucao do sistema

  for gl in sup_gls:

    ordem.remove(gl)
    ordem.append(gl)

  #print(ordem)

  # Reorganizar Kg
  K= Kg[np.ix_(ordem, ordem)]

  # Reorganizar Fg
  F = Fg[ordem]

  #print(Fg, '\n')
  #print(F, '\n')

  dim = tam - ncon

  sparse_K = scipy.sparse.csr_matrix(K[:dim, :dim]) # matriz esparsa

  dl = spsolve(sparse_K, F[:dim, 0])

  desloc = np.zeros((tam, 1)) # valores prescritos sao nulos

  for l in range(dim):
    desloc[ordem[l]] = dl[l]

  dx = []
  dy = []

  for d in range(len(desloc)):

    if d%2 == 0:
      dx.append(desloc[d][0])

    else:
      dy.append(desloc[d][0])

  #print(desloc)

  # 8) Processamento dos deslocamentos:

  # Salvar dados relativos ao no central inferior:

  mid_x = desloc[mid_idx*2]
  mid_y = desloc[mid_idx*2 + 1]

  mid_tot = np.sqrt((mid_x**2) + (mid_y**2))

  defmid.append(mid_tot*100) # Conversao para cm

  #print('\nDeslocamento no central inferior (No ', mid_idx, '): ')
  #print('Horiz. : ', mid_x)
  #print('Vert. : ', mid_y)
  #print('Total : ', mid_tot)
  #print('No. de elementos: ', nelem)

  # Plotagem da configuracao deformada:

  deform_plot(nodes, desloc, els, fator_def = 100, save = True)

# 9) Grafico de convergencia:

conv = plt.figure()

plt.rcdefaults()

print(defmid)

plt.plot(N_elmnt, defmid, color = 'black', marker = 'o')
plt.grid(axis = 'both', ls = ':')
plt.xlim(0, 1900)
plt.ylim(0.1, 0.5)
plt.xlabel('Numero de elementos')
plt.ylabel('Deslocamento total (cm)')
plt.tight_layout()
plt.savefig('conv_plot.png')

# 10) Tensoes horizontais:

mid_nodes = []

for i in range(nnos):

  if nodes[i][0] == L/2:
    mid_nodes.append(i)

mid_els = []

for i, e in enumerate(els):
  count_mid_nodes = 0

  for j in range(3):

    if e[j] in (mid_nodes):
      count_mid_nodes += 1

      if count_mid_nodes == 2:

        mid_els.append(i)

mid_els.sort()

mid_inc = els[mid_els]

sigmaX = []

for i, e in enumerate(mid_inc):

    x1, y1 = nodes[e[0]]
    x2, y2 = nodes[e[1]]
    x3, y3 = nodes[e[2]]

    no1, no2, no3 = e[0], e[1], e[2]

    #print('No 1: ', x1, y1, ' (No geral ', e[0], ')')
    #print('No 2: ', x2, y2, ' (No geral ', e[1], ')')
    #print('No 3: ', x3, y3, ' (No geral ', e[2], ')')

    # Area:

    Ae = 0.5*((x2*y3) + (x1*y2) + (x3*y1) - (x2*y1) - (x1*y3) - (x3*y2))

    #print('Area: ', Ae)

    b1 = y2 - y3
    b2 = -y1 + y3
    b3 = y1 - y2

    c1 = -x2 + x3
    c2 = x1 - x3
    c3 = -x1 + x2

    # Matriz deslocamento-deformacao:

    B = (1/(2*Ae))*np.array([[b1, 0, b2, 0, b3, 0],
                              [0, c1, 0, c2, 0, c3],
                              [c1, b1, c2, b2, c3, b3]])

    loc_d = np.array([[desloc[no1*2, 0]],
                      [desloc[no1*2+1, 0]],
                      [desloc[no2*2, 0]],
                      [desloc[no2*2+1, 0]],
                      [desloc[no3*2, 0]],
                      [desloc[no3*2+1, 0]]])

    S = np.dot(np.dot(D, B), loc_d)/1000
    sigmaX.append(S[0,0])


sxMed = []

for i in range(len(sigmaX)-1):
  sigma_med = (sigmaX[i] + sigmaX[i+1])/2
  sxMed.append(sigma_med)

yy = np.linspace(0, h, len(sxMed))

sig = plt.figure()

plt.plot(sxMed, yy, color = 'black')
plt.ylim(-0.5, 3.5)
plt.xlim(-10, 10)
plt.grid(axis = 'both', ls = ':')
plt.xlabel('Tensao horizontal (MPa)')
plt.ylabel('Distancia vertical (m)')
plt.tight_layout()
plt.savefig('sigma_plot.png')
