'''
copia de version antigua de com55.py. hace el primer análisis de sensibilidad para la variación de altura
'''
##############################
###### LIBRERIAS IMPORTADAS
##############################
import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.ndimage import laplace
import json
import sys





plt.close('all')

'''

def calcular_coeficientes_potencial(coordenada, coordenada_im, K, h, R, Vol):
    # Convertimos las listas de coordenadas a arreglos de numpy para realizar operaciones vectorizadas
    coordenada = np.array(coordenada)
    coordenada_im = np.array(coordenada_im)
    
    largo = len(coordenada)
    
    # Calcular matriz D de distancias entre pares de coordenadas reales
    coord_diff = coordenada[:, np.newaxis, :] - coordenada[np.newaxis, :, :]
    D = np.linalg.norm(coord_diff, axis=2)
    np.fill_diagonal(D, 1)  # Aseguramos que no haya división por cero más adelante
    
    # Calcular matriz D_im de distancias entre pares de coordenadas y sus imágenes
    coord_diff_im = coordenada[:, np.newaxis, :] - coordenada_im[np.newaxis, :, :]
    D_im = np.linalg.norm(coord_diff_im, axis=2)
    np.fill_diagonal(D_im, 1)
    
    # Crear matriz de coeficientes P
    P = np.zeros((largo, largo))
    diagonal_indices = np.arange(largo)
    P[diagonal_indices, diagonal_indices] = K * np.log(2 * h / R)
    
    # Calcular los coeficientes de potencial fuera de la diagonal
    mask_off_diag = ~np.eye(largo, dtype=bool)
    P[mask_off_diag] = K * np.log(D_im[mask_off_diag] / D[mask_off_diag])
    
    # Vector de voltajes
    V = np.array([Vol] * largo) if isinstance(Vol, (int, float)) else np.array(Vol)
    
    # Cálculo de cargas
    Q = np.dot(np.linalg.inv(P), V)
    
    return Q
    
Q = calcular_coeficientes_potencial(coordenada, coordenada_im, K, R, h, Vol)
'''

##################################
########## DEFINICIÓN DE FUNCIONES
##################################
def convierte_a_radio(area, conversion=1, es_mcm=False):
    """
    Calcula el radio de un conductor a partir del área.
    Parámetros:
        area: Área del conductor.
        conversion: Factor de conversión (por ejemplo, de MCM a mm²).
        es_mcm: Indica si el área está en MCM (True) o no (False).
    Retorna:
        Radio en las mismas unidades que el área (por defecto en mm).
    """
    factor = conversion if es_mcm else 1
    return np.sqrt(area * factor / np.pi)

def calcula_lados(numero, separacion):
    """
    Calcula las distancias entre los subconductores según su disposición.
    Parámetros:
        numero: Número de subconductores.
        separacion: Distancia entre subconductores vecinos.
    Retorna:
        Lista con las longitudes de los lados del polígono.
    """
    lados = [separacion] * numero
    if numero == 4:
        lados.extend([separacion * np.sqrt(2)] * 2)
    return lados

def distancia_equivalente(lados):
    """
    Calcula la distancia geométrica media (D_eq) a partir de una lista de lados.
    Parámetros:
        lados: Lista de distancias entre subconductores.
    Retorna:
        Distancia geométrica media (D_eq).
    """
    return (ma.prod(lados))**(1 / len(lados))

def radio_eq(radio, numero, distancia):
    """
    Calcula el radio equivalente de un subconductor fasciculado.
    Parámetros:
        radio: Radio del subconductor individual.
        numero: Número de subconductores.
        distancia: Distancia equivalente entre subconductores.
    Retorna:
        Radio equivalente del conductor.
    """
    r_eq = distancia * (radio / distancia)**(1 / numero)
    return 1.09 * r_eq if numero == 4 else r_eq

def calculo_radio_eq(numero, area_sub, sep, conversion=1, es_mcm=False, es_cm=False):
    """
    Calcula el radio equivalente de un conductor fasciculado.
    Parámetros:
        numero: Número de subconductores.
        area_sub: Área de cada subconductor.
        sep: Separación entre subconductores vecinos.
        conversion: Factor de conversión del área si está en MCM.
        es_mcm: Indica si el área está en MCM.
        es_cm: Indica si la separación está en cm.
    Retorna:
        Tuple con el radio equivalente y el radio individual.
    """
    lados = calcula_lados(numero, sep)
    distancia = distancia_equivalente(lados)
    radio = convierte_a_radio(area_sub, conversion, es_mcm)
    if es_cm:
        radio /= 10  # Convertir de mm a cm si es necesario
    radio_equi = radio_eq(radio, numero, distancia)
    return radio_equi, radio


def mod(z1,z2):
    return np.sqrt((z1[0]-z2[0])**2 + (z1[1]-z2[1])**2)

def windDist(wndx, wndy, hr, y, alpha, uni='uniforme'):
    # Input: wndx: valocidad media en x, wndy: velocidad media en y, hr: altura de referencia
    # y: Matriz de alturas, alpha: coeficiente de rugosidad del piso
    if uni=='no uniforme':
        Wx = wndx*(y/hr)**alpha
        Wy = wndy
    elif uni=='uniforme':
        Wx = wndx
        Wy = wndy
    return Wx, Wy

# Encontrar los índices más cercanos
def encuentra_nodos(x0,y0):
    indice_x0 = (np.abs(x - x0)).argmin()
    indice_y0 = (np.abs(y - y0)).argmin()
    return indice_x0, indice_y0

def espacio(anchox,anchoy, nodox, nodoy):
    '''
    Función que entrega las dimensiones espaciales de la malla en base al ancho de separacion de los nodos y en base a
    la cantidad de nodos que se desean en ambas dimensiones
    '''
    sx = np.abs((anchox/2)*(nodox-1))
    sy = np.abs((anchoy)*(nodoy-1))
    return sx, sy

def malla(ancho, Sx, Sy, nodox = False, nodoy = False):
    '''
    Discretiza la malla y define a los  arreglos de distancia que definen la malla
    Si nodox y nodoy no son dados, asume que se calculan a partir de 'ancho'
    Si son dados, solo discretiza el largo x e y
    '''
    if nodox is False and nodoy is False:
        x = np.arange(-Sx, Sx+ancho, ancho)
        y = np.arange(Sy, 0-ancho, -ancho)
        nodox = len(x)
        nodoy = len(y)
    else:
        x = np.linspace(-Sx, Sx, nodox)
        y = np.linspace(Sy, 0, nodoy)
    return x, y, nodox, nodoy

def discretiza(min,x_coor,y_coor, nodox=False, nodoy=False, sx=False, sy=False):
    if (sx is False and sy is False) and (nodox is not False and nodoy is not False):
        Sx, Sy = espacio(min, min, nodox, nodoy)
        x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=nodox, nodoy=nodoy)
        distx=np.abs(x_coor)-Sx
        nodos_sx = int(distx/min) + 1
        disty=np.abs(y_coor)-Sy
        nodos_sy = int(disty/min) + 1
        assert np.abs(x_coor) < Sx, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sx+nodosx}, o bien ubica a menor distancia x_coor, la dif x-Sx es {distx}'
        assert np.abs(y_coor) < Sy, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sy+nodosy}, o bien ubica a menor distancia y_coor, la dif x-Sx es {disty}'
    elif (sx is not False and sy is not False) and (nodox is False and nodoy is False):
        Sx, Sy = sx, sy
        x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=False, nodoy=False)
    elif (sx is not False and sy is not False) and (nodox is not False and nodoy is not False):
        Sx, Sy = sx, sy
        x, y, nodosx, nodosy = malla(min, sx, sy, nodox=nodox,nodoy=nodoy)
    else:
        print('Los parámetros de discretización no están bien ajustados en la función')
    return x, y, nodosx, nodosy, Sx, Sy

def val_rhoA(val):
    r_s=val*np.cos((np.pi-np.pi)/2)
    r_n=val*np.cos((np.pi-0)/2)
    r_e=val*np.cos((np.pi-np.pi/2)/2)
    r_o=val*np.cos((np.pi-np.pi/2)/2)
    return r_s, r_n, r_e, r_o

#def Dist_ValA(ValA, malla, fixedpx, fixedpy, in_condct='no'):
def Dist_ValA(malla, rho_b, ValA, px, py, in_condct='no', copia='si'):
    '''
    Distribuye la densidad de carga rhoA alrededor del conductor que se asume abarca 4 nodos
    sur,norte, este y oeste donde el máximo está en sur y minimo en norte
    Fija la condición de borde
    Si in_condct es 'no' entonces distribuye alrededor del conductor
    Si in_condct es 'si' entonces distribuye en un nodo único
    Cualquier otro caso no es válido
    '''
    if in_condct == 'no':
        if copia=='no':
            rs, rn, re, ro = val_rhoA(ValA)
            malla[py+1, px] = rs
            malla[py-1, px] = rn
            malla[py, px+1] = re
            malla[py, px-1] = ro
        elif copia=='si':
            malla[py+1, px] = rho_b[py+1, px]
            malla[py-1, px] = rho_b[py-1, px]
            malla[py, px+1] = rho_b[py, px+1]
            malla[py, px-1] = rho_b[py, px-1]
    elif in_condct == 'noV':
        malla[py+1, px] = ValA
        malla[py-1, px] = ValA
        malla[py, px+1] = ValA
        malla[py, px-1] = ValA
    elif in_condct == 'si':
        malla[py, px] = ValA
    else:
        print('in_condct escogido es inválido')
    return malla



def calcular_campo_electrico_inicial(nodosx, nodosy, x, y, w, h, Q, K):
    # Inicializar matrices para las componentes y magnitud del campo eléctrico
    E = np.zeros((nodosy, nodosx))
    Exx = np.zeros((nodosy, nodosx))
    Eyy = np.zeros((nodosy, nodosx))
    # Calcular el campo eléctrico en cada punto de la malla
    for i in range(nodosx):
        for j in range(nodosy):
            Ex, Ey = 0, 0  # Inicializar componentes del campo eléctrico en (i, j)
            
            for z in range(len(Q)):
                N1x = x[i] - w[z]
                N1y = y[j] - h[z]
                N2y = y[j] + h[z]
                
                # Contribución de la carga z al campo eléctrico en el punto (i, j)
                Ex += Q[z] * (K) * (N1x / (N1x**2 + N1y**2)**(1) - N1x / (N1x**2 + N2y**2)**(1))
                Ey += Q[z] * (K) * (N1y / (N1x**2 + N1y**2)**(1) - N2y / (N1x**2 + N2y**2)**(1))
            
            # Almacenar los resultados en las matrices
            Exx[j, i] = Ex
            Eyy[j, i] = Ey
            E[j, i] = np.sqrt(Ex**2 + Ey**2)
    
    return Exx, Eyy, E

## Resolución para densidad de carga en el espacio
def rhoA(Sx, Jp, b, a,  m_ind, d):
    # Cálculo condición de borde densidad de carga conductor energizado
    # Depende del valor de Jp
    # Asume que existe efecto corona
    rho_i = (Sx*Jp)*10**(-3)/(np.pi*b*(100*a)*m_ind*(30*d + 9*np.sqrt(d/(100*a)))) # radio está en cm y rho_i en (c/m^3)
    return rho_i

def potencial_electrostático(f_point, f_value, X, Y, radio, ind, carga=None, st = 'noV', copiado='no'):
    '''
    Potencial inicial libre de iones obtenido  de CSM
    '''
    py, px = f_point
    Vm = np.zeros((nodosy,nodosx))
    Vm2 = Vm.copy()
    ## Condiciones de borde potencial inicial
    Vm[f_point] = f_value
    Vm[-1,:] = 0 # borde inferior
    Vm2[f_point] = f_value
    Vm2[-1,:] = 0 # borde inferior
    #Vm = Dist_ValA(f_value, Vm, px,py, in_condct=st)
    Vm = Dist_ValA(Vm, Vm, f_value, px, py, in_condct=st, copia=copiado)
    if carga is not None:
        # Calcular la distancia entre cada punto de la malla y el punto (w, h)
        Mod_coor = np.sqrt((X - w)**2 + (Y - h)**2)
        Mod_coor2 = np.sqrt((X - w)**2 + (Y + h)**2)
        # Calcular el potencial eléctrico en cada nodo de la malla
        #Vm = carga * 0.5 * K*(1/Mod_coor - 1/Mod_coor2)
        Vm = carga*K*(np.log(20/Mod_coor) - np.log(20/Mod_coor2))
        #Vm = calcular_potencial_con_borde_vectorizado(Exx, Eyy,dx, dy)
        #Vm[fixed_point] = Vol
        #Vm = Dist_ValA(f_value, Vm, px,py, in_condct=st)
    else:
        c1 = np.cosh(np.pi*(X-2*ind*Sx)/(2*(Sy)))
        c2 = np.cos(np.pi*Ya/(2*(Sy)))
        c22 = np.cos(np.pi*y_p/(2*(Sy)))
        c3 = np.cosh(np.pi*ind*Sx/(Sy))
        c4 = np.cos(np.pi*radio/(2*(Sy)))
        Vm = np.abs(f_value * np.log((c1-c2)/(c1+c2))/np.log((c3-c4)/(c3+c4)))
        Vm2 = np.abs(f_value * np.log((c1-c22)/(c1+c22))/np.log((c3-c4)/(c3+c4)))
    #Vm[f_point] = f_value
    Vm = Dist_ValA(Vm, Vm, f_value, px, py, in_condct=st, copia=copiado)
    Vm[f_point] = f_value
    Vm[-1,:] = 0 # borde inferior
    Vm2[f_point] = f_value
    Vm2[-1,:] = 0 # borde inferior
    return Vm,Vm2
'''
def convergencia(rho_actual,rho_anterior, tol, epsilon=10**(-8)):
    # Calcula la diferencia relativa
    diferencia_relativa = np.abs(rho_actual - rho_anterior) / np.maximum.reduce([np.abs(rho_actual), np.abs(rho_anterior), np.ones_like(rho_actual)*epsilon])
    # Verifica si todas las diferencias relativas son menores a la tolerancia
    condicion = np.all(diferencia_relativa < tol)
    # Norma infinita de las diferencias relativas
    max_diferencia = np.linalg.norm(diferencia_relativa, ord=np.inf)
    return condicion, max_diferencia
'''
# Función para verificar si los valores tienen una tendencia horizontal
def verificar_convergencia(valores, ventana=110, umbral_variacion=6e-7):
    """
    Detecta si los valores oscilan en torno a un eje horizontal y si la variación es mínima.
    
    Parámetros:
        valores (list): Lista de valores.
        ventana (int): Número de valores recientes para evaluar el promedio.
        umbral_variacion (float): Umbral de variación aceptable entre promedios consecutivos.
        
    Retorna:
        bool: True si los valores han alcanzado una tendencia horizontal constante, False si no.
        float: Promedio actual de los últimos valores (None si no se cumplen las condiciones).
        float: Variación actual entre el último valor y el promedio actual (None si no se cumplen las condiciones).
    """
    if len(valores) < ventana:
        return None, None, None  # No se puede evaluar sin suficientes valores
    
    # Obtener los últimos 'ventana' valores
    valores_recientes = np.array(valores[-ventana:])
    promedio_actual = np.mean(valores_recientes)
    
    # Obtener el último valor
    ultimo_valor = valores[-1]
    
    # Calcular la variación actual (diferencia entre el último valor y el promedio actual)
    variacion_actual = np.abs(ultimo_valor - promedio_actual)
    
    # Inicializar dev_i con un valor por defecto
    dev_i = None
    
    # Comparar con la variación anterior (si existe)
    if hasattr(verificar_convergencia, 'ultima_variacion'):
        variacion_anterior = verificar_convergencia.ultima_variacion
        # Si la diferencia entre las dos variaciones es menor que el umbral, se considera convergencia
        dev_i = np.abs(variacion_actual - variacion_anterior)
        if dev_i < umbral_variacion:
            return True, promedio_actual, dev_i
    
    # Guardar la variación actual para la próxima iteración
    verificar_convergencia.ultima_variacion = variacion_actual
    
    return False, promedio_actual, dev_i

# Función personalizada para formatear números
def formatear_numero(valor):
    if np.abs(valor) < 1e-2 or np.abs(valor) > 1e3:
        return f"{valor:.3e}"  # Notación científica para valores extremos
    else:
        return f"{valor:.4f}"  # Notación decimal con 4 decimales

def convergencia(rho_actual, rho_anterior, tol):
    # Calcula la diferencia absoluta
    diferencia_absoluta = np.abs(rho_actual - rho_anterior)
    # Norma infinita de la diferencia absoluta
    max_diferencia = np.linalg.norm(diferencia_absoluta, ord=np.inf)
    # Verifica si la diferencia máxima es menor a la tolerancia
    condicion = max_diferencia < tol
    return condicion, max_diferencia

def evaluar_estado(dist_rho1, rhop0, diff_list, tol, px, py, ventana=110, umbral_variacion=1e-6, rho_hist=None):
    """
    Evalúa las condiciones de convergencia y estabilidad para el algoritmo iterativo.

    Parámetros:
        dist_rho1 (np.ndarray): Red actual.
        rhop0 (np.ndarray): Red de la iteración anterior.
        diff_list (list): Lista de diferencias relativas históricas.
        tol (float): Tolerancia para la convergencia de la red.
        px, py (int): Coordenadas del nodo a excluir en la máscara.
        ventana (int): Tamaño de la ventana para verificar tendencia.
        umbral_variacion (float): Umbral de variación para tendencia.
        rho_hist (list): Historial de redes iterativas (opcional).

    Retorna:
        bool: Si se alcanzó convergencia.
        dict: Información sobre el estado (promedio, desviación, diferencia relativa, etc.).
    """
    # Verificar convergencia por diferencia relativa
    condicion, diff = convergencia(dist_rho1, rhop0, tol)
    # Agregar la red actual al historial (si corresponde)
    if rho_hist is not False:
        rho_hist.append(dist_rho1.copy())
    # Agregar la diferencia relativa a la lista histórica
    if diff_list is not False:
        diff_list.append(diff)
    # Verificar tendencia al promedio
    tendencia, promedio, dev_dif = verificar_convergencia(diff_list, ventana, umbral_variacion)
    # Verificar que todos los nodos internos sean diferentes de cero excepto el nodo (px, py)
    mascara = np.ones_like(dist_rho1[1:-1, 1:-1], dtype=bool)  # Máscara con todos los valores True
    mascara[py - 1, px - 1] = False  # Excluir el nodo específico
    todos = np.all(dist_rho1[1:-1, 1:-1][mascara] != 0)
    # Evaluar condiciones finales de convergencia
    convergencia_total = (condicion and todos) or (tendencia and todos)
    # Construir información del estado
    estado = {
        "condicion": condicion,
        "diff": diff,
        "tendencia": tendencia,
        "promedio": promedio,
        "dev_dif": dev_dif,
        "todos": todos
    }
    return convergencia_total, estado


'''
def convergencia(rho_actual, rho_anterior, tol, estabilidad=15, delta_promedio=1e-9):
    """
    Verifica la convergencia basada en la diferencia máxima y la estabilidad del promedio de las últimas diferencias.
    
    Parámetros:
        rho_actual (ndarray): Valores actuales de la red.
        rho_anterior (ndarray): Valores anteriores de la red.
        tol (float): Tolerancia para la diferencia máxima.
        estabilidad (int): Número de diferencias a considerar para verificar constancia.
        delta_promedio (float): Tolerancia para evaluar si el promedio de las diferencias es constante.
    
    Retorna:
        condicion (bool): Verdadero si se cumple el criterio de convergencia.
        max_diferencia (float): Máxima diferencia absoluta entre `rho_actual` y `rho_anterior`.
    """
    # Calcula la diferencia absoluta
    diferencia_absoluta = np.abs(rho_actual - rho_anterior)
    # Norma infinita de la diferencia absoluta
    max_diferencia = np.linalg.norm(diferencia_absoluta, ord=np.inf)
    # Lista para almacenar diferencias históricas (debe mantenerse fuera de esta función)
    if not hasattr(convergencia, "hist_diffs"):
        convergencia.hist_diffs = []
    # Agrega la nueva diferencia a la lista
    convergencia.hist_diffs.append(max_diferencia)
    # Mantiene solo las últimas `estabilidad` diferencias
    if len(convergencia.hist_diffs) > estabilidad:
        convergencia.hist_diffs.pop(0)
    # Verifica si hay suficientes valores para evaluar la estabilidad
    if len(convergencia.hist_diffs) < estabilidad:
        return False, max_diferencia
    # Calcula el promedio de las diferencias recientes
    promedio_actual = np.mean(convergencia.hist_diffs)
    # Verifica si el promedio es constante dentro del umbral delta_promedio
    dif_mean = np.abs(promedio_actual - np.mean(convergencia.hist_diffs[:-1]))
    print(f'diferencia de los promedios es: {dif_mean}')
    if np.abs(promedio_actual - np.mean(convergencia.hist_diffs[:-1])) < delta_promedio:
        return True, max_diferencia
    # Verifica si la diferencia máxima es menor a la tolerancia
    condicion = max_diferencia < tol
    return condicion, max_diferencia
'''

def dev_1sentido(Ex, Ey, ep0, rhoi_1, rhoj_1, dx, dy):
    alpha = Ey*ep0/(2*dy) + Ex*ep0/(2*dx)
    beta = ep0*(Ey*rhoi_1/dy + Ex*rhoj_1/dx)
    dis = alpha**2 + beta+0j
    G = np.abs(-alpha + np.sqrt(dis))
    #G = np.sqrt(dis)
    #G[np.isnan(G)] = 0
    return G
'''
def dev_central(Ex, Ey,ep0, rhodx, rhoix, rhoiy, rhosy, dx, dy):
    d_rho_dx = (rhodx - rhoix) / (2 * dx)
    d_rho_dy = (rhoiy - rhosy) / (2 * dy)
    G = np.sqrt(-ep0 * (Ex * d_rho_dx + Ey * d_rho_dy)+0j)
    return G
'''
def dev_central(a,b,c, sig=1):
    G = np.sqrt(b**2-4*a*c+0j)
    if sig==1:
        return (-b + G)/(2*a)
    elif sig==-1:
        return (-b - G)/(2*a)
def dev_triple_di(Ex, Ey, ep0, rhoys, rhoyi, rhox, dx, dy):
    alpha  = Ex*ep0/(2*dx)
    beta = -ep0*(Ey*rhoyi/(2*dy) - (Ex*rhox/dx + Ey*rhoys/(2*dy)))
    G = np.abs(-alpha + np.sqrt(alpha**2 + beta+0j))
    return G
def dev_triple_si(Ex, Ey, ep0, rhoxd, rhoxi, rhoy, dx, dy):
    alpha  = Ey*ep0/(2*dy)
    beta = -ep0*(Ex*rhoxd/(2*dx) - (Ey*rhoy/dy + Ex*rhoxi/(2*dx)))
    G = np.abs(-alpha + np.sqrt(alpha**2 + beta +0j))
    return G


## Ocupa el negativo del gradiente del potencial calculado
def calcular_campo_electrico(V, dx, dy):
    # Inicializar las matrices para los componentes del campo eléctrico
    Ex = np.zeros_like(V)
    Ey = np.zeros_like(V)
    
    # Calcular las derivadas centrales (excepto en los bordes)
    Ey[1:-1, :] = (V[2:, :] - V[:-2, :]) / (2 * dy)
    Ex[:, 1:-1] = -(V[:, 2:] - V[:, :-2]) / (2 * dx)
    
    # Bordes: Aproximaciones hacia adelante y hacia atrás
    Ey[0, :] = (V[1, :] - V[0, :]) / dy  # Borde superior
    Ey[-1, :] = (V[-1, :] - V[-2, :]) / dy  # Borde inferior
    Ex[:, 0] = -(V[:, 1] - V[:, 0]) / dx  # Borde izquierdo
    Ex[:, -1] = -(V[:, -1] - V[:, -2]) / dx  # Borde derecho
    
    # Calcular el módulo del campo eléctrico en cada punto
    E_magnitud = np.sqrt(Ex**2 + Ey**2)
    
    return Ex, Ey, E_magnitud


# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
#def algoritmo_rho_v(V, rho_ini, dx, dy, windx, windy, max_iter_rho, Jplate, rho_A, visu, met):
def algoritmo_rho_V(V, rho_ini, rho_b, max_iter_rho, dx, dy, mov, wx,wy, val, pos, condi='si', copia='si',
                     rho_h=False, visu=15, nombre=None, muestra=False, fi=30):
    py1 ,px1 = pos
    rho1 = rho_ini.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
    rho1 = rho1.astype(complex)
    rho1 = Dist_ValA(rho1, rho_b, val, px1, py1, in_condct=condi,copia=copia)
    # En un comienzo se tiene el campo eléctrico electrostático
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    Ewx = Exxi + wx/mov
    Ewy = Eyyi + wy/mov
    a = 1/epsilon0
    b = 0
    # Se define a la malla rho_b con las condiciones de borde
    #rho_b[fixed_point] = rho_A
    #rhop_0 =impone_BC(rhop_0, rho_b, val1, px1, py1, in_cond=condi, copia=copia) # rhop_0 sería una red con condición de borde unicamente en la posición del conductor
    # COndición de desplazamiento de la densidad de carga debido al viento
    #rho_b[-1,:] = np.abs(Jplate/(mov*np.sqrt((Eyyi[-1,:]+windy[-1,:]/mov)**2+(Exxi[-1,:]+windx[-1,:]/mov)**2+0j))) # se postula que la dirección de E será solamente vertical
    #rho_b[-1,:] = np.abs(-(epsilon0/dy**2)*(V[-2,:]-2*V[-3,:]))
    #rho_b[-1,:] = 0

    #rho1[0,:]=0
    #rho1[-1,:]=0
    #rho1[:,0]=0
    #rho1[:,-1]=0
    RHO1 = rho1.copy()
    RHO1 =RHO1.astype(complex)
    difer = []
    if rho_h is not None:
        rho_historial = [rho1.copy()]
    else:
        rho_historial = None
    rho1, rhist, dif_list = update_rho_vectorized(max_iter_rho, dx, dy, Ewx, Ewy, RHO1, rho_b, val, px1, py1, a, b, sig=1, condi = condi,
                     copia=copia, rho_hist=rho_historial, nombre=nombre, est_gra=muestra, fi=fi)
    rho1 = extrapolate_borders(X,Y,rho1,num_layers=30)
    '''
    for iteration in range(max_iter_rho):
        rho_0 = rho1.copy()
        rho1 = update_rho_vectorized(RHO1, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_b, met = met)
        #rho1  = np.abs(rho1)
        # Criterio de convergencia
        # Define un pequeño valor para evitar divisiones por cero
         # Guardar estado actual para la animación
        rho_historial.append(rho1.copy())
        condicion,diff = convergencia(rho1, rho_0, 0.2)
        #print(r'Tolerancia $\rho_1$ = '+str(diff))
        if condicion:
            print(f"Convergencia para rho alcanzada en la iteración {iteration}")
            print(f'Dif relativa rho: {diff}')
            break
        #print(f'Dif relativa rho: {diff}')
        difer.append(diff)
        RHO1 = rho1
    #print(r'Tolerancia rho_1 = '+str(diff))
    #print(r'última iteración = '+str(iteration))
    difer_global.append(difer)
    '''
    if muestra is not False:
        visualizar_evolucion(X, Y, rhist, dif_list, titulo_base="Evolución densidad de carga", figura_num=visu, pausa=0.005)
    rhist =  False
    return np.real(rho1)


def visualizar_evolucion(X, Y, rho1_historial, list_dif, titulo_base="Evolución densidad de carga", figura_num=15, pausa=0.1):
    """
    Función para animar la evolución de una matriz a lo largo de iteraciones.

    Parámetros:
    - X, Y: Coordenadas de la malla.
    - rho1_historial: Lista de matrices que representan la evolución de la densidad de carga.
    - titulo_base: Título base para el gráfico.
    - pausa: Tiempo de pausa entre cada cuadro de la animación (en segundos).
    """
    if rho1_historial is not False:
        plt.figure(figura_num)  # Especificar el número de la figura
        fig, axis = plt.subplots(num=figura_num)  # Asegura que se use esa figura
        # Configurar el primer frame
        pcm = axis.pcolormesh(X, Y, np.abs(rho1_historial[0]), cmap='viridis', shading='auto', norm=LogNorm())
        cbar = plt.colorbar(pcm, ax=axis)
        cbar.set_label(r'Densidad de carga $C/m^3$', fontsize=11)
        # Iterar sobre el historial para actualizar la animación
        for i, rho in enumerate(rho1_historial):
            pcm.set_array(np.real(rho).ravel())
            if i==0:
                axis.set_title(f"{titulo_base} - Iteración: {i + 1}")
            else:
                axis.set_title(f"{titulo_base} - Iteración: {i + 1}, dif_absoluta: {formatear_numero(list_dif[i-1])}")
            plt.draw()
            plt.pause(pausa)
        # Suponiendo que tienes una lista diff_list
        plt.figure(figura_num+1)
        plt.plot(list_dif, label="Diferencia máxima por iteración")
        #plt.axhline(y=tol, color='r', linestyle='--', label="Tolerancia")
        plt.xlabel("Iteración")
        plt.ylabel("Diferencia absoluta (norma infinita)")
        plt.title("Evolución de la diferencia absoluta")
        plt.legend()
        plt.grid(True)
        plt.show()

# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
#def update_rho_vectorized(rhoini, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_bound, met=2):
def update_rho_vectorized(iteraciones, dx, dy, Ewx, Ewy, rho_iterado, rho_b, val, px, py, A, B, sig=1, condi = 'si',
                     copia='si', rho_hist=False, nombre=None, est_gra=False, fi=30):
    #rho_new = np.copy(rho)
    '''
    rho = rhoini.copy()
    rho  = rho.astype(complex)
    #rho[fixed_point] = rho_bound[fixed_point]  # Mantener la condición en el punto central
    #rho[posy_conductor+1, posx_conductor] = rho_bound[posy_conductor+1, posx_conductor]
    #rho[posy_conductor-1, posx_conductor] = rho_bound[posy_conductor-1, posx_conductor]
    #rho[posy_conductor, posx_conductor+1] = rho_bound[posy_conductor, posx_conductor+1]
    #rho[posy_conductor, posx_conductor-1] = rho_bound[posy_conductor, posx_conductor-1]
    #rho[-1, :] = rho_bound[-1,:]  # Mantener la condición en el borde inferior
    #Ewx = Ex + wx/mov
    #Ewy = Ey + wy/mov
    cx =fixed_point[1] # índice posición en x
    cy =fixed_point[0] # índice posición en y
    rhoN = rho.copy()
    rhoN = rhoN.astype(complex)
    if met==2:
        rhoN[cy:, cx:] = dev_1sentido(Ewx[cy:, cx:], Ewy[cy:, cx:], epsilon0, rho[cy-1:-1, cx:], rho[cy:, cx-1:-1], dx, dy) # 1er cuadrante
        rhoN[:cy+1, cx:] = dev_1sentido(Ewx[:cy+1, cx:], Ewy[:cy+1, cx:], epsilon0, rho[1:cy+2, cx:], rho[:cy+1, cx-1:-1], dx, dy) # 2do cuadrante
        rhoN[:cy+1, :cx+1] = dev_1sentido(Ewx[:cy+1, :cx+1], Ewy[:cy+1, :cx+1], epsilon0, rho[1:cy+2, :cx+1], rho[:cy+1, 1:cx+2], dx, dy) # 3er cuadrante
        rhoN[cy:, :cx+1] = dev_1sentido(Ewx[cy:, :cx+1], Ewy[cy:, :cx+1], epsilon0, rho[cy-1:-1, :cx+1], rho[cy:, 1:cx+2], dx, dy) # 4to cuadrante
    elif met == 0:
        # Calcular las derivadas en el interior (diferhorencias centrales) usando slicing
        # Ecuación diferencial discretizada en el interior
        rhoN[1:-1, 1:-1] = dev_central(Ewx[1:-1, 1:-1], Ewy[1:-1, 1:-1], epsilon0, rho[1:-1, 2:], rho[1:-1, :-2], rho[2:, 1:-1], rho[:-2, 1:-1], dx, dy)
        #rho[1:-1, 1:-1] = np.where(G >= 0, np.sqrt(G), 0)
        # Condiciones en los bordes
        #Gi = -epsilon0 * (Ewx[1:-1, 0] * (rho[1:-1, 1] - rho[1:-1, 0]) / dx + Ewy[1:-1, 0] * (rho[2:, 0] - rho[:-2, 0]) / (2 * dy))
        #rhoN[1:-1,0] = dev_triple_di(Ewx[1:-1,0], Ewy[1:-1,0], epsilon0, rho[:-2, 0], rho[2:, 0], rho[1:-1, 1], dx, dy) # borde izquierdo
        #Gd = -epsilon0 * (Ewx[1:-1, -1] * (rho[1:-1, -1] - rho[1:-1, -2]) / dx + Ewy[1:-1, -1] * (rho[2:, -1] - rho[:-2, -1]) / (2 * dy))
        #rhoN[1:-1,-1] = dev_triple_di(Ewx[1:-1,-1], Ewy[1:-1,-1], epsilon0, rho[:-2, -1], rho[2:, -1], rho[1:-1, -2], dx, dy) # borde derecho
        #Gs = -epsilon0 * (Ewx[0, 1:-1] * (rho[0, 2:] - rho[0, :-2]) / (2 * dx) + Ewy[0, 1:-1] * (rho[1, 1:-1] - rho[0, 1:-1]) / dy)
        #rhoN[0, 1:-1] = dev_triple_si(Ewx[0, 1:-1], Ewy[0, 1:-1], epsilon0, rho[0, 2:], rho[0, :-2], rho[1, 1:-1], dx, dy) # borde superior
        #rhoN[-1, 1:-1] = dev_triple_si(Ewx[-1, 1:-1], Ewy[-1, 1:-1], epsilon0, rho[-1, 2:], rho[-1, :-2], rho[-2, 1:-1], dx, dy) # borde inferior
        #Geis = -epsilon0 * (Ewx[0, 0] * (rho[0, 1] - rho[0, 0]) / dx + Ewy[0, 0] * (rho[1, 0] - rho[0, 0]) / dy)
        #rhoN[0, 0] = dev_1sentido(Ewx[0, 0], Ewy[0, 0], epsilon0, rho[1, 0], rho[0, 1], dx, dy) # esquina superior izquierda
        #Geds = -epsilon0 * (Ewx[0, -1] * (rho[0, -1] - rho[0, -2]) / dx + Ewy[0, -1] * (rho[1, -1] - rho[0, -1]) / dy)
        #rhoN[0, -1] = dev_1sentido(Ewx[0, -1], Ewy[0, -1], epsilon0, rho[1, -1], rho[0, -2], dx, dy) # esquina superior derecha
        #Geii = -epsilon0 * (Ewx[-1, 0] * (rho[-1, 1] - rho[-1, 0]) / dx + Ewy[-1, 0] * (rho[-1, 0] - rho[-2, 0]) / dy)
        #rhoN[-1, 0] = dev_1sentido(Ewx[-1, 0], Ewy[-1, 0], epsilon0, rho[-2, 0], rho[-1, 1], dx, dy) # esquina inferior izquierda
        #Gedi = -epsilon0 * (Ewx[-1, -1] * (rho[-1, -1] - rho[-1, -2]) / dx + Ewy[-1, -1] * (rho[-1, -1] - rho[-2, -1]) / dy)
        #rhoN[-1, -1] = dev_1sentido(Ewx[-1, -1], Ewy[-1, -1], epsilon0, rho[-2, -1], rho[-1, -2], dx, dy) # esquina inferior derecha
    # Aplicar condiciones de borde
    #rhoN[fixed_point] = rho_bound[fixed_point]  # Mantener la condición en el punto central
    rhoN[posy_conductor+1, posx_conductor] = rho_bound[posy_conductor+1, posx_conductor]
    rhoN[posy_conductor-1, posx_conductor] = rho_bound[posy_conductor-1, posx_conductor]
    rhoN[posy_conductor, posx_conductor+1] = rho_bound[posy_conductor, posx_conductor+1]
    rhoN[posy_conductor, posx_conductor-1] = rho_bound[posy_conductor, posx_conductor-1]
    rhoN[-1, :] = rho_bound[-1,:]  # Mantener la condición en el borde inferior
    #rhoN = np.abs(rhoN)
    #rhoN = extrapolate_borders(X, Y, rhoN)
    
    return rhoN
    '''
    dist_rho1 =np.zeros_like(rho_iterado) # asegura que los bordes sean nulos 
    dist_rho1 = dist_rho1.astype(complex)
    diff_list = []
    promedios = []
    desviaciones = []
    todos = False
    for i in range(iteraciones):
        if i == 0:
            C = delta_Vp(rho_iterado, Ewx, Ewy, dx, dy)
            rhop0 = rho_iterado.copy()
            dist_rho1[1:-1,1:-1] = dev_central(A, B, C, sig=sig)
            dist_rho1 = Dist_ValA(dist_rho1, rho_b, val, px, py, in_condct=condi, copia=copia)
        else:
            rhop0 = dist_rho1.copy()
            C = delta_Vp(dist_rho1, Ewx, Ewy, dx, dy)
            dist_rho1[1:-1,1:-1] = dev_central(A, B, C, sig=sig)
            dist_rho1 = Dist_ValA(dist_rho1, rho_b, val, px, py, in_condct=condi, copia=copia)
        #dist_rho1 = extrapolate_borders(X,Y, dist_rho1)
        #condicion,diff = convergencia(dist_rho1, rhop0,1e-8)
        # Evaluar condiciones de estado
        convergencia_total, estado = evaluar_estado(dist_rho1, rhop0, diff_list, tol=9e-07, px=px, py=py, rho_hist=rho_hist)
        # Registrar promedios y desviaciones
        if estado["promedio"] is not None and estado["dev_dif"] is not None:
            promedios.append(estado["promedio"])
            desviaciones.append(estado["dev_dif"])
        if est_gra is not False:
            if estado["promedio"] is not None and estado["dev_dif"] is not None:
                print(f'La tendencia al promedio es: {estado["tendencia"]}')
                print(f'Promedio actual: {estado["promedio"]}')
                print(f'Desviación del promedio actual: {estado["dev_dif"]}')
                print('....................')
        if convergencia_total:
            print(f'Diferencia relativa {nombre}: {estado["diff"]}')
            print(f'Promedio actual: {estado["promedio"]}')
            print(f'Desviación del promedio actual: {estado["dev_dif"]}')
            print(f'Término del algoritmo en iteración: {i}')
            break
        if i == iteraciones - 1:
            print(f'No hay convergencia de {nombre}: {estado["diff"]}')
    if est_gra is not False:
        plt.figure(fi)
        plt.plot(promedios, label=f'Promedio iteración {i}')
        plt.plot(desviaciones, label=f'Desviación iteración {i}')
        plt.xlabel('Iteraciones')
        plt.ylabel('Valor')
        plt.title('Evolución del Promedio y Desviación')
        plt.legend()
        plt.grid(True)
        plt.show()
    return dist_rho1, rho_hist, diff_list

def delta_Vp(rho_grid, ewx, ewy, dx, dy):
    '''
    Calcula el coeficiente constante de la ecuación de modelamiento de la densidad de carga espacial
    fac: factor (+1 o-1) el cual es asignado segun el tipo de conductor tratado.
    Deja  los bordes sin modificación
    '''
    d_rho_dx = (rho_grid[1:-1, 2:] - rho_grid[1:-1, :-2]) / (2 * dx)
    d_rho_dy = (rho_grid[2:, 1:-1] - rho_grid[:-2, 1:-1]) / (2 * dy)
    Coef = (d_rho_dx*ewx[1:-1,1:-1] + d_rho_dy*ewy[1:-1,1:-1])
    return Coef # Coef es de tamaño (n-2)x(m-2)

print('Densidad de carga preliminar alrededor del conductor calculado')
# Mostrar la matriz en una gráfica con un mapa de colores


#### Calculo del potencial con cálculo vectorizado
#################################3

# Inicializar la función u y el lado derecho f en la malla
def funcion_f(rho):
    return -1*rho.copy()/epsilon0  # f(x, y) = rho/epsilon0

def apply_laplacian_asymmetric(u, dx, dy):
    laplacian = np.zeros_like(u)
    
    # Calcular los coeficientes para las diferencias finitas
    dx2 = dx ** 2
    dy2 = dy ** 2

    # Aplicar la discretización asimétrica
    laplacian[1:-1, 1:-1] = (
        (u[:-2, 1:-1] - 2 * u[1:-1, 1:-1] + u[2:, 1:-1]) / dy2 +  # Término en x
        (u[1:-1, :-2] - 2 * u[1:-1, 1:-1] + u[1:-1, 2:]) / dx2    # Término en y
    )
    # Entrega red con las mismas dimensiones que u con los bordes con valores nulos
    return laplacian

# Método Jacobi para resolver el sistema de ecuaciones
def jacobi_step(u, f_rhs, dx, dy):
    dx2 = dx**2
    dy2 = dy**2
    unew = np.copy(u)
    unew[1:-1, 1:-1] = (1/(2*(dx2+dy2))) * (
        dy2*(u[0:-2, 1:-1] + u[2:, 1:-1]) + dx2*(u[1:-1, 0:-2] + u[1:-1, 2:])
        - f_rhs[1:-1, 1:-1]*dx2*dy2
    )
    return unew

# Definir las condiciones de borde usando una matriz 2D
def apply_boundary_conditions_2D(u, V_boundary, Exx=None, Eyy=None, dx=None, dy=None):
    """
    Aplica las condiciones de borde a la malla u utilizando la matriz V_boundary.
    También aplica condiciones de Neumann si se proporcionan Exx, Eyy, dx y dy.
    """
    # Borde izquierdo (Neumann si Exx está definido)
    if Exx is not None and dx is not None:
        #u[:, 0] = u[:, 1] + dx * (Exx[:, 0]+windx[:,0]/mov)
        u[:-1, 0] = u[:, 1] + dx * (Exx[:, 0])
    else:
        u[:-1, 0] = V_boundary[:-1, 0]  # Dirichlet

    # Borde derecho (Neumann si Exx está definido)
    if Exx is not None and dx is not None:
        #u[:, -1] = u[:, -2] - dx * (Exx[:, -1]+windx[:,-1]/mov)
        u[:-1, -1] = u[:, -2] - dx * (Exx[:, -1])
    else:
        u[:-1, -1] = V_boundary[:-1, -1]  # Dirichlet

    # Borde superior (Neumann si Eyy está definido)
    if Eyy is not None and dy is not None:
        #u[0, :] = u[1, :] + dy * (Eyy[0, :]+windy[0,:]/mov)
        u[0, :] = u[1, :] + dy * (Eyy[0, :])
    else:
        u[0, :] = V_boundary[0, :]  # Dirichlet

    # Borde inferior (siempre Dirichlet)
    u[-1, :] = V_boundary[-1, :]
    u[fixed_point] = 0

    return u


# Función modificada de v_cycle
def update_v(u, f_rhs, dx, dy, n_cycles=1, fixed_point=None, fixed_value=None, V_boundary=None, Exx=None, Eyy=None):
    if u.shape[0] <= 3:  # Resolver directamente cuando la grilla es muy pequeña
        return jacobi_step(u, f_rhs, dx, dy)

    for _ in range(n_cycles): 
        # Pre-suavizado
        u = jacobi_step(u, f_rhs, dx, dy)

        # Imponer condiciones de borde
        if V_boundary is not None:
            u = apply_boundary_conditions_2D(u, V_boundary, Exx=None, Eyy=None, dx=dx, dy=dy)
        
        # Imponer condición fija en el punto (si está definida)
        if fixed_point is not None and fixed_value is not None:
            u[fixed_point] = fixed_value
    
    return u

def densidad_voltaje(Volt, rho, Rho_val, pos, cond='si', cop='no'):
    py, px = pos
    # Calcula la densidad de cargas usando la ecuación de Poisson según el potencial ya calculado
    Rho_new = np.zeros_like(rho)
    Rho_new = -epsilon0*apply_laplacian_asymmetric(Volt, dx, dy)
    #Rho_new = -epsilon0*laplace(Volt)
    Rho_new = Dist_ValA(Rho_new, rho, Rho_val, px,py, in_condct=cond, copia=cop)
    Rho_new  = extrapolate_borders(X,Y, Rho_new)
    '''
    Rho_new[0,:] = rho[0,:]
    Rho_new[-1,:] = rho[-1,:]
    Rho_new[:,0] = rho[:,0]
    Rho_new[:,-1] = rho[:,-1]
    '''
    return Rho_new


###############################################################
###############################################################
# Función para interpolar los bordes sin considerar los valores nulos (ceros)
def interpolate_border(X, Y, grid, border, num_layers=4):
    if border == "bottom":
        coords_x, coords_y, values = X[-(num_layers+1):-1, 1:-1], Y[-(num_layers+1):-1, 1:-1], grid[-(num_layers+1):-1, 1:-1]
        target_x, target_y = X[-1, 1:-1], Y[-1, 1:-1]
    elif border == "top":
        coords_x, coords_y, values = X[1:num_layers+1, 1:-1], Y[1:num_layers+1, 1:-1], grid[1:num_layers+1, 1:-1]
        target_x, target_y = X[0, 1:-1], Y[0, 1:-1]
    elif border == "left":
        coords_x, coords_y, values = X[1:-1, 1:num_layers+1], Y[1:-1, 1:num_layers+1], grid[1:-1, 1:num_layers+1]
        target_x, target_y = X[1:-1, 0], Y[1:-1, 0]
    elif border == "right":
        coords_x, coords_y, values = X[1:-1, -(num_layers+1):-1], Y[1:-1, -(num_layers+1):-1], grid[1:-1, -(num_layers+1):-1]
        target_x, target_y = X[1:-1, -1], Y[1:-1, -1]
    else:
        raise ValueError("El borde debe ser 'top', 'bottom', 'left' o 'right'.")
    
    # Aplanar las matrices para usar griddata
    coords = np.column_stack((coords_x.flatten(), coords_y.flatten()))
    values_flat = values.flatten()

    # Eliminar puntos con valores NaN
    valid_mask = ~np.isnan(values_flat)
    coords = coords[valid_mask]
    values_flat = values_flat[valid_mask]

    # Usar griddata para interpolar
    target_coords = np.column_stack((target_x.flatten(), target_y.flatten()))
    interpolated_values = griddata(coords, values_flat, target_coords, method='nearest')

    return interpolated_values
# Función para tratar las esquinas con un promedio simple de los valores adyacentes
def corner_averages(X, Y, grid):
    # Esquinas: top-left, top-right, bottom-left, bottom-right
    # Cada esquina se calcula como la media de los tres puntos más cercanos
    top_left = np.mean([grid[1, 1], grid[0, 1], grid[1, 0]])  # top-left corner
    top_right = np.mean([grid[1, -2], grid[0, -2], grid[1, -1]])  # top-right corner
    bottom_left = np.mean([grid[-2, 1], grid[-1, 1], grid[-2, 0]])  # bottom-left corner
    bottom_right = np.mean([grid[-2, -2], grid[-1, -2], grid[-2, -1]])  # bottom-right corner
    # Devolver los valores de las esquinas como escalares
    return top_left, top_right, bottom_left, bottom_right
# Función para rellenar la red con los valores extrapolados

def extrapolate_borders(X, Y, grid, num_layers=4):
    # Interpolar para cada uno de los bordes (sin las esquinas)
    #top_interpolation = interpolate_border(X, Y, grid, "top", num_layers)
    bottom_interpolation = interpolate_border(X, Y, grid, "bottom", num_layers)
    #left_interpolation = interpolate_border(X, Y, grid, "left", num_layers)
    #right_interpolation = interpolate_border(X, Y, grid, "right", num_layers)

    # Rellenar la red con los valores extrapolados en los bordes
    #grid[0, 1:-1] = top_interpolation  # Llenar la fila superior
    grid[-1, 1:-1] = bottom_interpolation  # Llenar la fila inferior
    #grid[1:-1, 0] = left_interpolation  # Llenar la columna izquierda
    #grid[1:-1, -1] = right_interpolation  # Llenar la columna derecha
    #print("Shape of grid:", grid.shape)

    # Aplicar el tratamiento de las esquinas (solo después de llenar los bordes)
    top_left, top_right, bottom_left, bottom_right = corner_averages(X, Y, grid)

    # Asignar los valores de las esquinas
    grid[0, 0] = top_left
    grid[0, -1] = top_right
    grid[-1, 0] = bottom_left
    grid[-1, -1] = bottom_right

    return grid
################################################################

# Resolver usando FMG
# Algoritmo que obtien la red de potencial eléctrico en base a los valores de Rho
def algoritmo_V_rho(V, rho1, dx, dy, fixed_point, fixed_value, max_iter):
    f_rhs = funcion_f(rho1)
    #V_b =  Vm.copy()
    V_b = np.zeros_like(V)
    for iteration in range(max_iter):
        Vold = V.copy()
        V = update_v(V, f_rhs, dx,dy,fixed_point=fixed_point, fixed_value=fixed_value, V_boundary=V_b, Exx=None,Eyy=None)
        condicion, diff  = convergencia(V, Vold, 1e-6)
        if condicion:
            print(f"Convergencia alcanzada para V en la iteración {iteration}")
            break
    return V
'''
# ALGORITMO PRINCIPAL
# Se calcula primero el potencial incial en base a CSM
## Cálculo de potencial eléctrico inicial en base a la carga calculada por CSM
Vmi,Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, R, 0, carga=Q, st='noV')
##### Resolución ecuación de continuidad
Jp = 1988*10**(-8) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
# Condiciones de borde
rho_i = rhoA(Sx, Jp, mov, R, m, delta)
rho_inicial = np.zeros((nodosy, nodosx), dtype=complex)
rho_boundary = np.zeros_like(rho_inicial)
con_condct = 'si'
copiado = 'no'
rho_boundary = Dist_ValA(rho_boundary, rho_inicial, rho_i, posx_conductor, posy_conductor, in_condct=con_condct, copia=copiado)
visualizacion = 20
histl = True
Exxini, Eyyini, Em = calcular_campo_electrico(Vmi, dx, dy) # campo electrostático inicial
Vm = np.zeros_like(Vmi)
difer_global = []
vis = 15
it_global = 20
for n in range(it_global):
    Volder = Vm.copy()
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    if n==0:
        print('Primer rho')
        rho_n = algoritmo_rho_V(Vmi, rho_inicial, rho_boundary, max_iter_rho, dx, dy, mov, windx, windy, rho_i, fixed_point,
                                  condi=con_condct, copia=copiado, rho_h=histl, visu=visualizacion,  nombre='rho', muestra=None)
    else:
        print('----------------------')
        print(f'Comienza nueva actualización de rho número {n+1}')
        #rho_n = algoritmo_rho_v(Vm+Vmi, rho_n, dx, dy, windx, windy, max_iter_rho, Jp, rho_i, vis, met=0)
        #rho_n = algoritmo_rho_V(Vmi+Vm, rho_inicial, rho_boundary, max_iter_rho, dx, dy, mov, windx, windy, rho_i, fixed_point,
        #                          condi=con_condct, copia=copiado, rho_h=None, visu=visualizacion, nombre='rho', muestra=None)
        #rho_n = algoritmo_rho_V(Vmi+Vm, rho_n, rho_boundary, max_iter_rho, dx, dy, mov, windx, windy, rho_i, fixed_point,
        #                          condi=con_condct, copia=copiado, rho_h=None, visu=visualizacion, nombre='rho', muestra=None)
        rho_n = densidad_voltaje(Vm, rho_boundary, rho_i,fixed_point, cond=con_condct, cop=copiado)
    Vm = algoritmo_V_rho(Vm, rho_n, dx, dy, fixed_point, 0, max_iter)
    condicion, max_diff = convergencia(Vm, Volder, 3.1e-9)
    if condicion:
        print(f"Convergencia alcanzada para V en la iteración {n}")
        print(r'Diferencia relativa V y Vold: '+str(max_diff))
        break
    else:
        print(r'Diferencia relativa V y Vold: '+str(max_diff))
    

print(r'Diferencia relativa V y Vold: '+str(max_diff))
Di = np.zeros((len(difer_global), 2))
Di[:,0] = [np.mean(difer_global[i]) for i in range(len(difer_global))]
Di[:,1] = [np.std(difer_global[i]) for i in range(len(difer_global))]
#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
conver = 'Incompleto'
print('Potencial calculado')
##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica
Vol_def = Vm + Vmi
Edefx, Edefy, Edef = calcular_campo_electrico(Vol_def, dx, dy)
J = rho_n*mov*np.sqrt((Edefx+(windx/mov))**2 + (Edefy+(windy/mov))**2)
Jy = rho_n*mov*(Edefy+(windy/mov))
pxl,pyl = encuentra_nodos(0, l)
Ei = Edef[pyl,:] # Magnitud Campo eléctrico a nivel de piso
Ji = J[pyl,:] # Densidad de corriente a nivel de piso
Jave = np.mean(Ji)
print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:])/(10**(-9)))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')
if Jave/Jp < 1e-5:
    conver = 'Completo'
#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')
print(f'El estado de algoritmo está {conver}')
'''


def calcular_potencial_inicial(fixed_point, fixed_value, X, Y, R, Q):
    """Calcula el potencial inicial basado en el método CSM."""
    Vmi, Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, R, 0, carga=Q, st='noV')
    return Vmi

def inicializar_densidad(Sx, Jp, mov, R, m, delta, nodosy, nodosx, posx_conductor, posy_conductor, con_condct='si', copiado='no'):
    """Inicializa la densidad de carga y las condiciones de borde."""
    rho_i = rhoA(Sx, Jp, mov, R, m, delta)
    rho_inicial = np.zeros((nodosy, nodosx), dtype=complex)
    rho_boundary = np.zeros_like(rho_inicial)
    rho_boundary = Dist_ValA(rho_boundary, rho_inicial, rho_i, posx_conductor, posy_conductor, in_condct=con_condct, copia=copiado)
    return rho_inicial, rho_boundary, rho_i

def iterar_potencial(Vmi, rho_inicial, rho_boundary, rho_i, fixed_point, dx, dy, mov, windx, windy, max_iter_rho, max_iter,
                      it_global, visu=15, con_condct='si', copiado='no', must = False, fi=30):
    """Itera para calcular el potencial eléctrico y la densidad de carga."""
    Vm = np.zeros_like(Vmi)
    difer_global = []
    print('Campo eléctrico  electrostático libre iones calculado')
    for n in range(it_global):
        Volder = Vm.copy()
        if n == 0:
            print('Primer rho')
            rho_n = algoritmo_rho_V(
                Vmi, rho_inicial, rho_boundary, max_iter_rho, dx, dy, mov, windx, windy, rho_i, fixed_point,
                condi=con_condct, copia=copiado, rho_h=False, visu=visu, nombre='rho', muestra=must, fi=fi
            )
        else:
            print(f'Comienza nueva actualización de rho número {n + 1}')
            rho_n = densidad_voltaje(Vm, rho_boundary, rho_i, fixed_point, cond=con_condct, cop=copiado)
        Vm = algoritmo_V_rho(Vm, rho_n, dx, dy, fixed_point, 0, max_iter)
        condicion, max_diff = convergencia(Vm, Volder, 3.1e-9)
        if condicion:
            print(f"Convergencia alcanzada para V en la iteración {n}")
            break
        else:
            print(f'Diferencia relativa V y Vold: {max_diff}')
            print('*****************************')
    if n==it_global-1:
        print(f'se alcanzó el máximo de iteraciones: {n}')
    return Vm, rho_n

def calcular_resultados_finales(Vm, Vmi, rho_n, dx, dy, mov, windx, windy, l):
    """Calcula el campo eléctrico, densidad de corriente y verifica la convergencia."""
    Vol_def = Vm + Vmi
    Exxini, Eyyini, Eini =calcular_campo_electrico(Vmi, dx, dy)
    Edefx, Edefy, Edef = calcular_campo_electrico(Vol_def, dx, dy)
    Ei = Edef[encuentra_nodos(0, l)[1],  :] # Magnitud campo eléctrico nivel de piso
    J = rho_n * mov * np.sqrt((Edefx + (windx / mov))**2 + (Edefy + (windy / mov))**2)
    Jave = np.mean(J[-1,:]) # Promedio densidad de corriente a nivel  de piso
    Ji = J[encuentra_nodos(0, l)[1], :]  # Densidad de corriente a nivel de l
    Campo_ini = [Exxini, Eyyini, Eini]
    Campo_fin = [Edefx, Edefy, Edef]
    return Jave, Ji, Campo_ini, Campo_fin, Ei

def ejecutar_algoritmo(fixed_point, fixed_value, X, Y, R, Q, Sx, mov, m, delta, nodosy, nodosx, posx_conductor, posy_conductor,
                        dx, dy, windx, windy, max_iter_rho, max_iter, it_global, l, visualizacion, Jp_inicial,
                          tolerancia=[1-1e-5,1+1e-5], condct='si', copy='no', Muestra=False):
    """Ejecuta el algoritmo completo ajustando Jp hasta cumplir la condición de convergencia."""
    Jp = Jp_inicial
    convergencia_lograda = False
    cont = 0
    fi = 30
    while not convergencia_lograda and cont < 10:
        print('----------------------------------')
        print(f'Intentando con Jp = {Jp}')
        # Paso 1: Calcular potencial inicial
        Vmi = calcular_potencial_inicial(fixed_point, fixed_value, X, Y, R, Q)
        # Paso 2: Inicializar densidades
        rho_inicial, rho_boundary, rho_i = inicializar_densidad(Sx, Jp, mov, R, m, delta, nodosy, nodosx, posx_conductor, posy_conductor)
        # Paso 3: Iterar para calcular Vm y rho
        Vm, rho_n = iterar_potencial(Vmi, rho_inicial, rho_boundary, rho_i, fixed_point, dx, dy, mov, windx, windy, max_iter_rho, max_iter,
                                      it_global, visu=visualizacion, con_condct=condct, copiado=copy, must=Muestra, fi = fi)
        #fi += 1
        visualizacion += 2
        # Paso 4: Calcular resultados finales
        Vol_def = Vm+Vmi
        Jave, Ji, Campo_ini, Campo_fin, Ei = calcular_resultados_finales(Vm, Vmi, rho_n, dx, dy, mov, windx, windy, l)
        print(f'Jp promedio calculado: {Jave * 1e9} nA/m^2')
        print(f'Jp promedio propuesto: {Jp * 1e9} nA/m^2')
        cont += 1
        if (Jave / Jp) <= tolerancia[1] and (Jave / Jp) >= tolerancia[0]:
            convergencia_lograda = True
            print('Convergencia alcanzada!')
        else:
            resto = Jave - Jp
            Jp += resto  # Ajustar Jp
    if cont==10:
        print('Hubieron {10} iteraciones sin conseguir el mismo Jp')
    else:
        print('Algoritmo completado con éxito.')
    return Campo_ini, Vmi, rho_n, Vm, Vol_def, Campo_fin, Ei, Ji, Jave

######################################
############ PANEL DE CONTROL
######################################

# Parámetros constructivos
#Di = 20 * 10**(-2) # (m) diametro equivalente del conductor
numero_cond =  3
separacion = 20 #  cm
es_cm = True
conversion = 0.5067  # Conversión de MCM a mm²
area_cond = 140 # mm^2
es_mcm = False
Req, R = calculo_radio_eq(numero_cond, area_cond, separacion, conversion=1, es_mcm=es_mcm, es_cm=es_cm)
Req /= 10
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol = 400000 # (V) voltaje de la línea
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
mov = 7*10**(-4) # (m^2/Vs)
m = 0.7 # (AD) factor de rugosidad
#m = [0.4, 0.7, 1]
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0= 303  # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
#Y_coor = [16, 19, 22]
y_coor = 22 # (m) coordenada y ubicación conductor
x_coor = 0 # (m) coordenada x ubicación conductor
l = 1 # (m) distancia desde el suelo a la altura de interés de influencia de campo
fixed_value = Vol

# Discretización
Sx = 15 # (m) media longitud del plano de tierra 
Sy = 35 # (m) altura del área de estudio respecto de tierra
nodosx, nodosy = False, False # número de nodos en ambas dimensiones
# Definir coordenadas de la malla en x y y
x,y, nodosx, nodosy, Sx, Sy = discretiza(Req,x_coor,y_coor, nodox=nodosx, nodoy=nodosy, sx=Sx, sy=Sy)  
X,Y =np.meshgrid(x,y)
dx = np.abs(x[1] - x[0]) # (m) distancia física entre nodos en cada columna
dy = np.abs(y[1] - y[0]) # (m) distancia física entre nodos en cada fila
posx_conductor, posy_conductor = encuentra_nodos(x_coor, y_coor)
fixed_point = (posy_conductor, posx_conductor)

# Ajuste de viento
wndx = 0 # m/s
#wndx = [0, 5, 10]
wndy = 0 # m/s
modo = 'uniforme'
coef_drag = 0.1
wndx1 = np.ones((nodosy, nodosx)) * wndx
wndy1 = np.ones((nodosy, nodosx)) * wndy
windx, windy = windDist(wndx1, wndy1, l, Y, coef_drag, uni=modo) # Se calcula la distribución de viento

# Parámetros de iteración y convergencia
max_iter_rho = 400 # Número máximo de iteraciones para el cálculo de la densidad de carga
max_iter = 230 # Número máximo de iteraciones para el cálculo del potencial en el espacio
it_global = 10 # Número máximo de iteraciones para el cálculo del potencial definitivo
Jp_inicial = 198e-8 # Valor inicial dado, se espera que se recupere este valor al momento de conseguir la convergencia
TolDev = 1e-2 # Tolerancia de convergencia de Jp
Tolerancia = [1-TolDev, 1+TolDev] # Rango de tolerancia de Jp
visualizacion = 15 #  número de la ventana donde se muestra la evolución de la densidad de carga
in_condct = 'si' # Condición de que se ubica la densidad de carga alrededor del conductor o solo en un nodo central
copiado = 'no' # Condición de que los valores sean copiados o no al momento de mantener las condiciones de borde en las mallas
histl = False # Si se desea crear la lista de evolución de la densidad de carga
mostrar = False # Permite mostrar los datos de la convergencia de rho inicial

###########################################################
##### ALGORITMO PRINCIPAL Y EJECUCIÓN
###########################################################

coordenada = [(x_coor,y_coor)]
coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
largo = len(coordenada)
D= np.zeros((largo,largo))
D_im= np.zeros((largo,largo))
for i in range(len(coordenada)):
    for j in range(len(coordenada)):
        if i != j:
            D[i][j]=mod(coordenada[i], coordenada[j])
            D_im[i][j]=mod(coordenada[i], coordenada_im[j])

## Coeficientes de potencial
P = np.zeros((largo,largo)) # matriz coeficientes
for i in range(largo):
    for j in range(largo):
        if j==i:
            P[i][j] = K*np.log(2*h[i]/Req)
        else:
            if largo !=1 :
                P[i][j] = K*np.log(D_im[i][j]/D[i][j])
#V = [Vol,Vol] # voltajes sobre los conductores reales
V = [Vol]
Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores
Campo_ini, Vmi, rho_n, Vm, Vol_def, Campo_fin, Ei, Ji, Jave = ejecutar_algoritmo(fixed_point, fixed_value, X, Y, Req, Q, Sx, mov, m, delta, nodosy,
                                                                            nodosx, posx_conductor, posy_conductor, dx, dy, windx, windy, max_iter_rho,
                                                                                max_iter, it_global, l, visualizacion, Jp_inicial=Jp_inicial,
                                                                                tolerancia=Tolerancia, condct=in_condct, copy=copiado, Muestra=mostrar)

Exxini, Eyyini, Em = Campo_ini
Edefx, Edefy, Edef = Campo_fin
'''
E_camp = []
J_camp = []
for ine  in range(len(wndx)):
    P0 =101.3 # (kPa) Presión del aire a nivel de mar
    T0= 303  # (Kelvin) Temperatura de 25°C
    Pr =  90 # Presión del aire
    Tr= 290 # (Kelvin) Temperatura del sitio
    delta = Pr*T0/(P0*Tr) # () densidad del aire
    #Y_coor = [16, 19, 22]
    y_coor = 22 # (m) coordenada y ubicación conductor
    x_coor = 0 # (m) coordenada x ubicación conductor
    l = 1 # (m) distancia desde el suelo a la altura de interés de influencia de campo
    fixed_value = Vol

    # Discretización
    Sx = 15 # (m) media longitud del plano de tierra 
    Sy = 35 # (m) altura del área de estudio respecto de tierra
    nodosx, nodosy = None, None # número de nodos en ambas dimensiones
    # Definir coordenadas de la malla en x y y
    x,y, nodosx, nodosy, Sx, Sy = discretiza(R,x_coor,y_coor, nodox=nodosx, nodoy=nodosy, sx=Sx, sy=Sy)  
    X,Y =np.meshgrid(x,y)
    dx = np.abs(x[1] - x[0]) # (m) distancia física entre nodos en cada columna
    dy = np.abs(y[1] - y[0]) # (m) distancia física entre nodos en cada fila
    posx_conductor, posy_conductor = encuentra_nodos(x_coor, y_coor)
    fixed_point = (posy_conductor, posx_conductor)

    # Ajuste de viento
    wndx = 0 # m/s
    wndy = 0 # m/s
    modo = 'uniforme'
    coef_drag = 0.1
    wndx1 = np.ones((nodosy, nodosx)) * wndx[ine]
    wndy1 = np.ones((nodosy, nodosx)) * wndy
    windx, windy = windDist(wndx1, wndy1, l, Y, coef_drag, uni=modo) # Se calcula la distribución de viento

    # Parámetros de iteración y convergencia
    max_iter_rho = 400 # Número máximo de iteraciones para el cálculo de la densidad de carga
    max_iter = 230 # Número máximo de iteraciones para el cálculo del potencial en el espacio
    it_global = 10 # Número máximo de iteraciones para el cálculo del potencial definitivo
    Jp_inicial = 198e-8 # Valor inicial dado, se espera que se recupere este valor al momento de conseguir la convergencia
    TolDev = 1e-2 # Tolerancia de convergencia de Jp
    Tolerancia = [1-TolDev, 1+TolDev] # Rango de tolerancia de Jp
    visualizacion = 15 #  número de la ventana donde se muestra la evolución de la densidad de carga
    in_condct = 'si' # Condición de que se ubica la densidad de carga alrededor del conductor o solo en un nodo central
    copiado = 'no' # Condición de que los valores sean copiados o no al momento de mantener las condiciones de borde en las mallas
    histl = False # Si se desea crear la lista de evolución de la densidad de carga
    mostrar = False # Permite mostrar los datos de la convergencia de rho inicial

    ###########################################################
    ##### ALGORITMO PRINCIPAL Y EJECUCIÓN
    ###########################################################

    coordenada = [(x_coor,y_coor)]
    coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
    h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
    w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
    largo = len(coordenada)
    D= np.zeros((largo,largo))
    D_im= np.zeros((largo,largo))
    for i in range(len(coordenada)):
        for j in range(len(coordenada)):
            if i != j:
                D[i][j]=mod(coordenada[i], coordenada[j])
                D_im[i][j]=mod(coordenada[i], coordenada_im[j])

    ## Coeficientes de potencial
    P = np.zeros((largo,largo)) # matriz coeficientes
    for i in range(largo):
        for j in range(largo):
            if j==i:
                P[i][j] = K*np.log(2*h[i]/R)
            else:
                if largo !=1 :
                    P[i][j] = K*np.log(D_im[i][j]/D[i][j])
    #V = [Vol,Vol] # voltajes sobre los conductores reales
    V = [Vol]
    Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores

    Campo_ini, Vmi, rho_n, Vm, Vol_def, Campo_fin, Ei, Ji, Jave = ejecutar_algoritmo(fixed_point, fixed_value, X, Y, R, Q, Sx, mov, m, delta, nodosy,
                                                                            nodosx, posx_conductor, posy_conductor, dx, dy, windx, windy, max_iter_rho,
                                                                                max_iter, it_global, l, visualizacion, Jp_inicial=Jp_inicial,
                                                                                tolerancia=Tolerancia, condct=in_condct, copy=copiado, Muestra=mostrar)

    Exxini, Eyyini, Em = Campo_ini
    Edefx, Edefy, Edef = Campo_fin
    E_camp.append(Ei)
    J_camp.append(Ji)
'''

##########################
######## GRAFICOS
#########################
# Lista gráficos
def grafE(num):
    plt.figure(num)
    #mod = np.sqrt(U**2+Ww**2)
    plt.quiver(X, Y, Exxini, Eyyini, Em, cmap='plasma', scale_units='xy')
    #plt.imshow(E, cmap='viridis', interpolation='none')
    # Agregar la barra de colores
    cbar = plt.colorbar()
    cbar.set_label(r'Magnitud campo eléctrico $kV/m$', fontsize=11)
    # Obtén los ticks actuales
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    # Mostrar el gráfico
    plt.title('Campo electrostático libre de iones', fontsize=15)
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    plt.tight_layout()

def grafV(num):
    plt.figure(num)
    #plt.figure(figsize=(6, 6))
    #plt.contourf(X, Y, Vm, levels=200, cmap='plasma')
    plt.pcolormesh(X, Y, Vmi, cmap='plasma', shading='auto',norm=LogNorm())
    cbar = plt.colorbar()
    cbar.set_label(r'Potencial $kV$')
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    plt.title('Potencial electrostático', fontsize=15)
    plt.tight_layout()

def grafRho(num):
    plt.figure(num)
    plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto',norm=LogNorm())
    #plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto')
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    plt.title('Densidad de carga final',fontsize=15)
    # Añadir una barra de colores para mostrar la escala
    cbar = plt.colorbar()
    cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)

def grafVf(nm):
    plt.figure(nm)
    #plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
    plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto',norm=LogNorm())
    plt.title('Potencial iónico',fontsize=15) 
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    # Añadir una barra de colores para mostrar la escala
    cbar = plt.colorbar()
    cbar.set_label(r'Potencial iónico $kV$')
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    plt.tight_layout
    #plt.show()

def grafVdef(num):
    plt.figure(num)
    #plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
    plt.pcolormesh(X, Y, Vol_def, cmap='plasma', shading='auto',norm=LogNorm())
    plt.title('Potencial definitivo',fontsize=15) 
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    # Añadir una barra de colores para mostrar la escala
    cbar = plt.colorbar()
    cbar.set_label(r'Potencial definitivo $kV$')
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    plt.tight_layout
    #plt.show()

def grafEf(nm):
    plt.figure(nm)
    # Graficar el campo vectorial 2D
    # Usar Edefx, Edefy para las componentes del campo eléctrico
    # Edef para la magnitud, para mapear colores en las flechas
    Q=plt.quiver(X, Y, Exxini, Eyyini, Em, cmap='plasma', scale_units='xy', angles='xy')
    #Q = plt.quiver(X, Y, Edefx, Edefy, Edef, cmap='plasma', scale_units='xy', angles='xy', scale=1)
    # Agregar la barra de colores
    cbar = plt.colorbar(Q)
    cbar.set_label(r'Magnitud campo eléctrico $(kV/m)$', fontsize=11)
    # Obtén los ticks actuales
    ticks = cbar.get_ticks()
    # Cambiar las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    # Mostrar etiquetas y título
    plt.xlabel(r'Distancia horizontal (m)', fontsize=11)
    plt.ylabel(r'Distancia vertical (m)', fontsize=11)
    plt.title(r'Magnitud de campo eléctrico a nivel de suelo $(V/m)$', fontsize=13)
    # Ajustar la disposición para evitar que los textos se solapen
    plt.tight_layout()
    plt.show()


def grafJ1(num):
    plt.figure(num)
    plt.plot(x[30:-30], Ji[30:-30]*(10**9))
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
    plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
    plt.grid(True)

def grafE1(num):
    plt.figure(num)
    plt.plot(x[30:-30], Ei[30:-30]/1000)
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
    plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV'])
    plt.grid(True)



gra = ['Eele', 'Vele', 'Rhof', 'Vf', 'Vdef', 'Ef', 'E1', 'J1']

def show_plot(graf, ven=1):
    ven = 10**(ven-1)
    for i in graf:
        if i == 'Eele':
            grafE(ven)
        elif i == 'Vele':
            grafV(ven+1)
        elif i == 'Rhof':
            grafRho(ven+2)
            '''
            plt.figure(4)
            plt.pcolormesh(X, Y, rho_alternativo, cmap='viridis', shading='auto')
            #plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto')
            plt.xlabel('Distancia horizontal (m)',fontsize=11)
            plt.ylabel('Distancia vertical (m)',fontsize=11)
            plt.title('Densidad de carga final alternativo',fontsize=15)
            # Añadir una barra de colores para mostrar la escala
            cbar = plt.colorbar()
            cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
            '''
        elif i == 'Vf':
            grafVf(ven+4)
        elif i == 'Vdef':
            grafVdef(ven+5)
        elif i == 'Ef':
            grafEf(ven+6)
        elif i == 'J1':
            grafJ1(ven+7)
        elif i == 'E1':
            grafE1(ven+8)
            fig = plt.figure(figsize=(12, 12))  # Tamaño ajustado para los subgráficos
            # Subgráfico 1: Densidad de carga iónica
            ax1 = fig.add_subplot(221, projection='3d')
            surf1 = ax1.plot_surface(X, Y, rho_n * 10**6, cmap='viridis', edgecolor='none')
            fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)
            ax1.set_title(r'Densidad de carga iónica')
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_zlabel(r'$\rho (\mu C/m^3)$')
            # Subgráfico 2: Potencial electrostático
            ax2 = fig.add_subplot(222, projection='3d')
            surf2 = ax2.plot_surface(X, Y, Vmi / 1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)
            ax2.set_title('Potencial electrostático')
            ax2.set_xlabel('X')
            ax2.set_ylabel('Y')
            ax2.set_zlabel('V(kV)')
            # Subgráfico 3: Potencial iónico
            ax3 = fig.add_subplot(223, projection='3d')
            surf3 = ax3.plot_surface(X, Y, Vm / 1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf3, ax=ax3, shrink=0.5, aspect=10)
            ax3.set_title('Potencial iónico')
            ax3.set_xlabel('X')
            ax3.set_ylabel('Y')
            ax3.set_zlabel('V(kV)')
            # Subgráfico 4: Potencial definitivo
            ax4 = fig.add_subplot(224, projection='3d')
            surf4 = ax4.plot_surface(X, Y, Vol_def / 1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf4, ax=ax4, shrink=0.5, aspect=10)
            ax4.set_title('Potencial definitivo')
            ax4.set_xlabel('X')
            ax4.set_ylabel('Y')
            ax4.set_zlabel('V(kV)')
            # Ajustar diseño
            plt.tight_layout()
        # Renderiza cada gráfico y continúa al siguiente
        plt.pause(0.1)  # Pausa breve para que el gráfico se renderice
    plt.show()
    


'''
    plt.figure(10)
    #plt.figure(figsize=(6, 6))
    #plt.contourf(X, Y, Vm, levels=200, cmap='plasma')
    plt.pcolormesh(X, Y, Vmi2, cmap='plasma', shading='auto',norm=LogNorm())
    cbar = plt.colorbar()
    cbar.set_label(r'Potencial $kV$')
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    plt.title('Potencial electrostático', fontsize=15)
    plt.tight_layout()
'''
    

    ########## 
show_plot(gra)
'''
plt.figure()
for j  in range(len(E_camp)):
    Ei = E_camp[j]
    plt.plot(x[30:-30], Ei[30:-30]/1000, label=f'$w_x$ = {wndx[j]} m/s, $|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV')
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
#plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l), fontsize=13)
plt.tight_layout()
plt.legend()
plt.grid(True)

plt.figure()
for i  in range(len(J_camp)):
    Ji = J_camp[i]
    plt.plot(x[30:-30], Ji[30:-30]*(10**9), label=f'$w_x$ = {wndx[i]} m/s, $J_p$ = {str(np.round(np.mean(Ji)*(10**9),3))} $nA/m^2$')
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
#plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l), fontsize=13)
plt.tight_layout()
plt.legend()
plt.grid(True)
#threading.Thread(target=show_plot).start()
plt.show()
'''