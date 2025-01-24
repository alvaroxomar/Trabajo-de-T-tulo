if __name__ == "__main__":
        ##############################
    ###### LIBRERIAS IMPORTADAS
    ##############################
    import numpy as np
    import matplotlib.pyplot as plt
    import math as ma
    from tkinter import messagebox
    import tkinter as tk
    from matplotlib.colors import LogNorm
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.interpolate import griddata
    from scipy.ndimage import laplace
    import json
    import sys
    import pandas as pd
    import os
    ########## FUNCIONES ###############
    def pressure_height(altura,tempe):
        P0  = 101325 # Pa
        M = 0.029 # Kg/mol
        g =  9.8 # m/s^2
        R0 = 8.314 # J/mol K
        P = P0*np.e**(-g*M*altura/(R0*tempe))
        return P/1000 # en kPa

    def mod(z1,z2):
        return np.sqrt((z1[0]-z2[0])**2 + (z1[1]-z2[1])**2)
    
    # Función personalizada para formatear números
    def formatear_numero(valor):
        if np.abs(valor) < 1e-2 or np.abs(valor) > 1e3:
            return f"{valor:.3e}"  # Notación científica para valores extremos
        else:
            return f"{valor:.4f}"  # Notación decimal con 4 decimales

    def obtener_parametro(param, tipo):
        """
        Verifica si un parámetro está vacío y lo asigna como None o lo convierte al tipo especificado.

        :param param: Valor del parámetro a verificar.
        :param tipo: Tipo al que se debe convertir (float o int).
        :return: Valor convertido o None si está vacío.
        """
        if param == "" or param is None:
            return None
        return tipo(param)
    
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
    
    def distancia_equivalente_2(sepa, numero):
        """
        Calcula la distancia geométrica media (D_eq) a partir de una lista de lados.
        Parámetros:
            sepa: Separación subconductores vecinos
        Retorna:
            Distancia geométrica media (D_eq).
        """
        return sepa/(2*np.sin(np.pi/numero))

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
        #if numero==1:
        #    r_eq = radio
        #else:
        #    r_eq = distancia * (radio / distancia)**(1 / numero)
        #return 1.09 * r_eq if numero == 4 else r_eq
        return (radio*numero*distancia**(numero-1))**(1/numero)
    
    def calculo_radio_eq(numero, area_sub, sepa, conversion=1, es_mcm=False, es_cm=False):
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
        #lados = calcula_lados(numero, sepa)
        #distancia = distancia_equivalente(lados)
        distancia_e = distancia_equivalente_2(sepa, numero)
        radio = convierte_a_radio(area_sub, conversion, es_mcm)
        if es_cm:
            radio /= 10  # Convertir de mm a cm si es necesario
        radio_equi = radio_eq(radio, numero, distancia_e)
        if radio_equi==0:
            radio_equi = radio
        return radio_equi, radio
    
        # Encontrar los índices más cercanos
    def encuentra_nodos(x, y, x0,y0):
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

    def malla(ancho, Sx, Sy, nodox = None, nodoy = None):
        '''
        Discretiza la malla y define a los  arreglos de distancia que definen la malla
        Si nodox y nodoy no son dados, asume que se calculan a partir de 'ancho'
        Si son dados, solo discretiza el largo x e y
        '''
        if nodox is None and nodoy is None:
            x = np.arange(-Sx, Sx+ancho, ancho)
            y = np.arange(Sy, 0-ancho, -ancho)
            nodox = len(x)
            nodoy = len(y)
        else:
            x = np.linspace(-Sx, Sx, nodox)
            y = np.linspace(Sy, 0, nodoy)
        return x, y, nodox, nodoy

    def discretiza(min, nodox=None, nodoy=None, sx=None, sy=None):
        if (sx is None and sy is None) and (nodox is not None and nodoy is not None):
            Sx, Sy = espacio(min, min, nodox, nodoy)
            x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=nodox, nodoy=nodoy)
            distx=np.abs(x_coor)-Sx
            nodos_sx = int(distx/min) + 1 # Cantidad de nodos que faltan para que que la posicion x del conductor esté dentro del ancho del espacio de estudio
            disty=np.abs(y_coor)-Sy
            nodos_sy = int(disty/min) + 1 # Cantidad de nodos que faltan para que que la altura del conductor esté dentro de la altura del espacio de estudio
            assert np.abs(x_coor) < Sx, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sx+nodosx}, o bien ubica a menor distancia x_coor, la dif x-Sx es {distx}'
            assert np.abs(y_coor) < Sy, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {2*(nodos_sy+nodosy)}, o bien ubica a menor distancia y_coor, la dif y-Sy es {disty}'
        elif (sx is not None and sy is not None) and (nodox is None and nodoy is None):
            Sx, Sy = sx, sy
            x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=None, nodoy=None)
        elif (sx is not None and sy is not None) and (nodox is not None and nodoy is not None):
            Sx, Sy = sx, sy
            x, y, nodosx, nodosy = malla(min, sx, sy, nodox=nodox,nodoy=nodoy)
        else:
            print('Los parámetros no están bien ajustados en la función')
        return x, y, nodosx, nodosy, Sx, Sy
    
    def windDist(wndx, wndy, hr, y, alpha, uni='uniforme'):
        # Input: wndx: valocidad media en x, wndy: velocidad media en y, hr: altura de referencia
        # y: Matriz de alturas, alpha: coeficiente de rugosidad del piso
        if uni=='gradiente':
            Wx = wndx*(y/hr)**alpha
            Wy = wndy
        elif uni=='uniforme':
            Wx = wndx
            Wy = wndy
        else:
            print(f"el código de modo de viento {uni} no es aceptado (uniforme/gradiente)")
        return Wx, Wy
    
    def val_rhoA(val):
        r_s=val*np.cos((np.pi-np.pi)/2)
        r_n=val*np.cos((np.pi-0)/2)
        r_e=val*np.cos((np.pi-np.pi/2)/2)
        r_o=val*np.cos((np.pi-np.pi/2)/2)
        return r_s, r_n, r_e, r_o

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
    
    def Hay_corona(Vol, Sep, r, delta, m, g0):
        '''
        Vol: en kV
        Sep: en cm
        r: en cm
        '''
        Ev= grad_sup(g0, m, delta, r)
        ev = Vol_crit(Ev, Sep, r)
        if ev <= np.abs(Vol):
            print(f"{Vol} kV >= {ev} kV, hay efecto corona")
            return True
        else:
            print(f"{Vol} kV < {ev} kV, no hay efecto corona")
            return False
        
    def grad_sup(g0, m, delta, r):
        Ev = g0*delta*m*(1+0.301/np.sqrt(delta*r))
        return  Ev
    
    def Vol_crit(Ev, Sep, r):
        ev = Ev*r*np.log(Sep/r)
        return ev

    def rhoA(Sx, Jp, b, a,  m_ind, d):
        # Cálculo condición de borde densidad de carga conductor energizado
        # Depende del valor de Jp
        # Asume que existe efecto corona
        rho_i = (Sx*Jp)*10**(-3)/(np.pi*b*(100*a)*m_ind*(30*d + 9*np.sqrt(d/(100*a)))) # radio está en cm y rho_i en (c/m^3)
        return rho_i
    
    def dev_1sentido(Ex, Ey, ep0, rhoi_1, rhoj_1, dx, dy):
        alpha = Ey*ep0/(2*dy) + Ex*ep0/(2*dx)
        beta = ep0*(Ey*rhoi_1/dy + Ex*rhoj_1/dx)
        dis = alpha**2 + beta+0j
        G = np.abs(-alpha + np.sqrt(dis))
        #G = np.sqrt(dis)
        #G[np.isnan(G)] = 0
        return G
    
    def dev_central(a,b,c, sig=1):
        G = np.sqrt(b**2-4*a*c+0j)
        if sig==1:
            return (-b + G)/(2*a)
        elif sig==-1:
            return (-b - G)/(2*a)
    
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
    
    ############ FUNCIONES PRINCIPALES ######################
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
    
    def f1(vol, vol_co):
        # Usado solo para caso de conductores fasciculados
        return 1.1339 - 0.1678*(vol/vol_co) + 0.03*(vol/vol_co)**2

    ## Ocupa el negativo del gradiente del potencial calculado
    def calcular_campo_electrico(V, Vo, dx, dy, posicion_cond, Evv, Vco, kapzov_borde=False, bundle=False):
        print(f"posicion_cond dentro de la funcion calcular_campo es {posicion_cond}")
        posy, posx = posicion_cond
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
        if kapzov_borde:
            if bundle:
                print(f"Hay conductores fasciculados")
                print(f"el voltaje aplicado es {Vo}")
                Evv *= f1(Vo, Vco)
                print(f"El nuevo gradiente superficial depende de {Vo/1000} kV y es {Evv/(10**5)} kV/cm")
            # Calcular el módulo del campo eléctrico en cada punto
            E_magnitud = np.sqrt(Ex**2 + Ey**2)
            chi_sur = E_magnitud[posy+1, posx]/Evv
            chi_norte = E_magnitud[posy-1, posx]/Evv
            chi_este = E_magnitud[posy, posx+1]/Evv
            chi_oeste = E_magnitud[posy, posx-1]/Evv
            Ex[posy+1, posx] /= chi_sur
            Ex[posy-1, posx] /= chi_norte
            Ex[posy, posx+1] /= chi_este
            Ex[posy, posx-1] /= chi_oeste
            Ey[posy+1, posx] /= chi_sur
            Ey[posy-1, posx] /= chi_norte
            Ey[posy, posx+1] /= chi_este
            Ey[posy, posx-1] /= chi_oeste
        # Calcular el módulo del campo eléctrico en cada punto
        E_magnitud = np.sqrt(Ex**2 + Ey**2)
        return Ex, Ey, E_magnitud
    
    def algoritmo_rho_V(V, Vo, Vco, rho_ini, rho_b, max_iter_rho, dx, dy, mov, wx, wy, val, pos, Evv, condi='si', copia='si',
                        rho_h=False, visu=15, nombre=None, muestra=False, fi=30, is_bundled=False):
        py1 ,px1 = pos
        rho1 = rho_ini.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
        rho1 = rho1.astype(complex)
        rho1 = Dist_ValA(rho1, rho_b, val, px1, py1, in_condct=condi, copia=copia)
        # En un comienzo se tiene el campo eléctrico electrostático
        Exxi, Eyyi, Em = calcular_campo_electrico(V, Vo, dx, dy, pos, Evv, Vco, kapzov_borde=True, bundle=is_bundled)
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

    def update_rho_vectorized(iteraciones, dx, dy, Ewx, Ewy, rho_iterado, rho_b, val, px, py, A, B, sig=1, condi = 'si',
                        copia='si', rho_hist=False, nombre=None, est_gra=False, fi=30):
        #rho_new = np.copy(rho)
        dist_rho1 =np.zeros_like(rho_iterado) # asegura que los bordes sean nulos 
        dist_rho1 = dist_rho1.astype(complex) #np.complex64
        diff_list = []
        promedios = []
        desviaciones = []
        todos = False
        for i in range(iteraciones):
            if i == 0:
                C = delta_Vp(rho_iterado, Ewx, Ewy, dx, dy)
                rhop0 = rho_iterado.copy()
                dist_rho1[1:-1,1:-1] = dev_central(A, B, C, sig=sig)
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

    def calcular_potencial_inicial(fixed_point, fixed_value, X, Y, R, Q):
        """Calcula el potencial inicial basado en el método CSM."""
        Vmi, Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, R, 0, carga=Q, st='noV')
        return Vmi

    def inicializar_densidad(Sx, Jp, mov, R, m, delta, Vol, yco, g0, nodosy, nodosx, posx_conductor, posy_conductor, con_condct='si', copiado='no'):
        """Inicializa la densidad de carga y las condiciones de borde."""
        Corona = Hay_corona(Vol/1000, yco*100, R*100, delta, m, g0)
        if Corona:
            rho_i = rhoA(Sx, Jp, mov, R, m, delta)
        else:
            rho_i = 0
        rho_inicial = np.zeros((nodosy, nodosx), dtype=complex)
        rho_boundary = np.zeros_like(rho_inicial)
        rho_boundary = Dist_ValA(rho_boundary, rho_inicial, rho_i, posx_conductor, posy_conductor, in_condct=con_condct, copia=copiado)
        return rho_inicial, rho_boundary, rho_i

    def iterar_potencial(Vmi, Vo, Vco, rho_inicial, rho_boundary, rho_i, fixed_point, dx, dy, mov, windx, windy, max_iter_rho, max_iter, Evv,
                        it_global, visu=15, con_condct='si', copiado='no', must = False, fi=30, is_bundled=False):
        """Itera para calcular el potencial eléctrico y la densidad de carga."""
        Vm = np.zeros_like(Vmi)
        difer_global = []
        print('Campo eléctrico  electrostático libre iones calculado')
        for n in range(it_global):
            Volder = Vm.copy()
            if n == 0:
                print('Primer rho')
                rho_n = algoritmo_rho_V(
                    Vmi, Vo, Vco, rho_inicial, rho_boundary, max_iter_rho, dx, dy, mov, windx, windy, rho_i, fixed_point, Evv,
                    condi=con_condct, copia=copiado, rho_h=False, visu=visu, nombre='rho', muestra=must, fi=fi, is_bundled=is_bundled
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

    def calcular_resultados_finales(Vm, Vo, Vco, Vmi, rho_n, dx, dy, mov, windx, windy, l, posicion_cond, Evv, is_bundled=False):
        """Calcula el campo eléctrico, densidad de corriente y verifica la convergencia."""
        Vol_def = Vm + Vmi
        Exxini, Eyyini, Eini =calcular_campo_electrico(Vmi, Vo, dx, dy, posicion_cond, Evv, Vco, kapzov_borde=False, bundle=False)
        Edefx, Edefy, Edef = calcular_campo_electrico(Vol_def, Vo, dx, dy, posicion_cond, Evv, Vco, kapzov_borde=True, bundle=is_bundled)
        Ei = Edef[encuentra_nodos(x, y, 0, l)[1],  :] # Magnitud campo eléctrico nivel de piso
        J = rho_n * mov * np.sqrt((Edefx + (windx / mov))**2 + (Edefy + (windy / mov))**2)
        Jave = np.mean(J[-1,:]) # Promedio densidad de corriente a nivel  de piso
        Ji = J[encuentra_nodos(x, y, 0, l)[1], :]  # Densidad de corriente a nivel de l
        Campo_ini = [Exxini, Eyyini, Eini]
        Campo_fin = [Edefx, Edefy, Edef]
        return Jave, Ji, Campo_ini, Campo_fin, Ei

    def ejecutar_algoritmo(fixed_point, fixed_value, Vco, X, Y, R, Q, Sx, mov, m, delta, g0, nodosy, nodosx, posx_conductor, posy_conductor, yco,
                            dx, dy, windx, windy, max_iter_rho, Evv, max_iter, it_global, l, visualizacion, Jp_inicial,
                            tolerancia=[1-1e-5,1+1e-5], condct='si', copy='no', Muestra=False, is_bundled=False):
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
            rho_inicial, rho_boundary, rho_i = inicializar_densidad(Sx, Jp, mov, R, m, delta, fixed_value, yco, g0, nodosy, nodosx, posx_conductor, posy_conductor)
            # Paso 3: Iterar para calcular Vm y rho
            Vm, rho_n = iterar_potencial(Vmi, fixed_value, Vco, rho_inicial, rho_boundary, rho_i, fixed_point, dx, dy, mov, windx, windy, max_iter_rho, max_iter, Evv,
                                        it_global, visu=visualizacion, con_condct=condct, copiado=copy, must=Muestra, fi = fi, is_bundled=is_bundled)
            #fi += 1
            visualizacion += 2
            # Paso 4: Calcular resultados finales
            Vol_def = Vm+Vmi
            Jave, Ji, Campo_ini, Campo_fin, Ei = calcular_resultados_finales(Vm, fixed_value, Vco, Vmi, rho_n, dx, dy, mov, windx, windy, l, fixed_point, Evv, is_bundled=is_bundled)
            print(f'Jp promedio calculado: {Jave * 1e9} nA/m^2')
            print(f'Jp promedio propuesto: {Jp * 1e9} nA/m^2')
            cont += 1
            if Jp != float(0):
                if (Jave / Jp) <= tolerancia[1] and (Jave / Jp) >= tolerancia[0]:
                    convergencia_lograda = True
                    print('Convergencia alcanzada!')
                else:
                    resto = Jave - Jp
                    Jp += resto  # Ajustar Jp
            else:
                convergencia_lograda = True
                print('No hay efecto corona por lo que no hay corriente de iones a nivel de suelo')
        if cont==10:
            print('Hubieron {10} iteraciones sin conseguir el mismo Jp')
        else:
            print('Algoritmo completado con éxito.')
        return Campo_ini, Vmi, rho_n, Vm, Vol_def, Campo_fin, Ei, Ji, Jave

    ####### Almacenamiento de datos en archivos xlsx y csv

    def guardar_datos_multihoja(distancias, campo_electrico, densidad_corriente, ruta=".", nombre_csv="datos.csv", nombre_excel="datos.xlsx", guarda=False):
        """
        Guarda listas en archivos CSV y Excel. Los datos se organizan en múltiples hojas en el archivo Excel.

        Args:
            distancias (list or np.array): Lista de distancias físicas.
            campo_electrico (list or np.array): Lista de valores del campo eléctrico.
            densidad_corriente (list or np.array): Lista de valores de densidad de corriente.
            x_coor (float): Posición en X.
            otros_datos (dict): Diccionario con otras hojas y datos, formato {"NombreHoja": DataFrame}.
            ruta (str): Ruta donde se guardarán los archivos.
            nombre_csv (str): Nombre del archivo CSV a generar.
            nombre_excel (str): Nombre del archivo Excel a generar.
            guarda (bool): Indicador para guardar o no los archivos.
        """
        if guarda:
            # Asegurarse de que la ruta existe
            if not os.path.exists(ruta):
                os.makedirs(ruta)

            # Crear DataFrame principal (Campo Eléctrico y Corriente)
            data_principal = {
                "Lateral Distance x[m]": distancias,
                "|E| [kV/m]": campo_electrico,
                "J [nA/m^2]": densidad_corriente}
            df_principal = pd.DataFrame(data_principal)

            # Crear las rutas completas para los archivos
            ruta_csv = os.path.join(ruta, nombre_csv)
            ruta_excel = os.path.join(ruta, nombre_excel)

            # Guardar como CSV (solo datos principales)
            if not os.path.exists(ruta_csv):
                df_principal.to_csv(ruta_csv, index=False)
                print(f"Datos principales guardados como '{ruta_csv}'.")
            else:
                print(f"Archivo CSV existente. No se sobrescribió: '{ruta_csv}'.")

            # Preparar los datos extra para cada hoja
            datos_conductores = {
                "Tipo": [tipo],
                "Nombre": [nombre],
                "Área (mm²)": [area],
                "Diámetro (mm)": [diametro],
                "Subconductores": [cantidad],
                "Separación (cm)": [sep],
                "Voltaje (kV)": [Vol/1000],
                "Posición en x (m)": [x_coor],
                "Posición en y (m)": [y_coor],
                "factor conductor": [m],
                "Jp inicial (nA/m^2)": [Jp_inicial],
                "Radio equivalente (cm)": [Req*100],
                "Gradiente superficial crítico (kV/cm)": Evv,
                "Potencial crítico corona (kV)": evv
            }
            df_conductores = pd.DataFrame(datos_conductores)

            datos_ambientales = {
                "Temperatura (°C)": [params["Temperatura (°C)"]],
                "Presión (kPa)": [Pres],
                "Viento x (m/s)": [viento_x],
                "Viento y (m/s)": [viento_y],
                "Rugosidad terreno": [rug],
                "Movilidad iónica (m2/kVs)": [mov*1000],
            }
            df_ambientales = pd.DataFrame(datos_ambientales)

            datos_dimensionamiento = {
                "Ancho (m)": [Sx],
                "Altura (m)": [Sy],
                "nodos x": [nodosx],
                "nodos y": [nodosy],
                "Medición (m)": [l],
                "dx (m)": [dx],
                "dy (m)": [dy]
            }
            df_dimensionamiento = pd.DataFrame(datos_dimensionamiento)

            datos_iteraciones = {
                "Max iter rho": [max_iter_rho],
                "Max iter V": [max_iter],
                "Max iter Gob": [it_global],
            }
            df_iteraciones = pd.DataFrame(datos_iteraciones)

            # Guardar como Excel con múltiples hojas
            with pd.ExcelWriter(ruta_excel, engine="openpyxl") as writer:
                # Hoja principal (Campo Eléctrico y Corriente)
                df_principal.to_excel(writer, index=False, sheet_name="Campo y Corriente")
                print("Hoja 'Campo y Corriente' añadida al archivo Excel.")

                # Hoja de Características de los conductores
                df_conductores.to_excel(writer, index=False, sheet_name="Conductores")
                print("Hoja 'Conductores' añadida al archivo Excel.")

                # Hoja de Características ambientales
                df_ambientales.to_excel(writer, index=False, sheet_name="Ambientales")
                print("Hoja 'Ambientales' añadida al archivo Excel.")

                # Hoja de Dimensionamiento y discretización
                df_dimensionamiento.to_excel(writer, index=False, sheet_name="Dim y discr")
                print("Hoja 'Dim y discr' añadida al archivo Excel.")

                # Hoja de Iteraciones
                df_iteraciones.to_excel(writer, index=False, sheet_name="Iteraciones")
                print("Hoja 'Iteraciones' añadida al archivo Excel.")
            print(f"Datos guardados como '{ruta_excel}'.")

    def guarda_en_carpeta(Vol, ycoor, nx, ny, x, Ei, Ji, ruta, guarda=False):
        """
        Genera los nombres de archivo y llama a la función guardar_datos para guardar los datos.

        Args:
            Vol (float): Voltaje aplicado.
            ycoor (float): Coordenada Y.
            nx (int): Número de nodos en X.
            ny (int): Número de nodos en Y.
            x (list or np.array): Lista de distancias físicas.
            Ei (list or np.array): Lista de valores del campo eléctrico.
            Ji (list or np.array): Lista de valores de densidad de corriente.
            guarda (bool): Indicador para guardar o no los archivos.
        """
        
        #carpeta = os.makedirs(ruta_destino)
        nom_csv = f"modelo_{Vol/1000}_{ycoor}_{nx}_{ny}.csv"
        nom_xlsx = f"modelo_{Vol/1000}_{ycoor}_{nx}_{ny}.xlsx"
        #guardar_datos(x, Ei, Ji, ruta=ruta_destino, nombre_csv=nom_csv, nombre_excel=nom_xlsx, guarda=guarda)
        guardar_datos_multihoja(x, Ei, Ji, ruta=ruta, nombre_csv=nom_csv, nombre_excel=nom_xlsx, guarda=guarda)
    
    ############### PARÁMETROS ###################
    # Capturar el diccionario JSON desde los argumentos
    if len(sys.argv) < 2:
        print("Se requiere un diccionario JSON como argumento.")
        sys.exit(1)

    # Cargar el diccionario desde el JSON
    try:
        params = json.loads(sys.argv[1])
    except (IndexError, json.JSONDecodeError) as e:
        print(f"Error al leer los parámetros: {e}")
    # Desestructurar el diccionario en variables específicas
    
    # Extraer parámetros necesarios
    #area = float(params["Área (mm²)"])
    tipo, nombre, areas = params["Conductor"]
    area =float(areas[0])
    diametro = float(areas[1])
    cantidad = int(params["Subconductores"])
    sep = float(params["Separación (cm)"])
    is_bundled = False
    if sep > 0:
        is_bundled = True
    es_cm = True # si separacion está en cm entonces es True
    es_mcm = False # si  aárea está en mm2 entonces es False
    conversion = 0.5067  # Conversión de MCM a mm²
    Vol = float(params["Voltaje (kV)"])*1000 # en Volt
    x_coor = float(params["Posición en x (m)"])
    y_coor = float(params["Posición en y (m)"])
    # Obtener los parámetros
    Sx = obtener_parametro(params["Ancho (m)"], float)  # (m) Media longitud del plano de tierra
    Sy = obtener_parametro(params["Altura (m)"], float)  # (m) Altura del área de estudio respecto de tierra
    nodosx = obtener_parametro(params["nodos x"], int)  # Número de nodos en x
    nodosy = obtener_parametro(params["nodos y"], int)  # Número de nodos en y
    mov = float(params["Movilidad iónica (m2/kVs)"])/1000 # m2/Vs
    Tem = float(params["Temperatura (°C)"]) + 273 # kelvin
    #Pres = float(params["Presión (kPa)"]) # kPa
    altura = float(params["Altitud (m)"]) # m
    Pres = pressure_height(altura, Tem) # kPa
    print(f"la altura es {altura} y la presion es {Pres} kPa")
    viento_x = float(params["Viento x (m/s)"])
    viento_y = float(params["Viento y (m/s)"])
    modo = str(params["Modo (str)"])
    rug =  float(params["Rugosidad terreno"])
    m = float(params["Estado sup cond 1"])
    l = float(params["Medición (m)"])
    gra = params["graficos"]
    print(f"los graficos seleccionados son {gra}")

    print(str(area))
    print(str(sep))
    Req, R = calculo_radio_eq(cantidad, area, sep, conversion=conversion, es_mcm=es_mcm, es_cm=es_cm) # están en cm
    print(f"radio eq {Req} y radio subconductor {R}")
    Req /= 100 # en m
    R /= 100
    print(f'{Req},{R} en m')
    
    max_iter_rho = int(params["Max iter rho"])
    max_iter = int(params["Max iter V"])
    it_global = int(params["Max iter Gob"])

    Jp_inicial = float(params["Jp (nA/m^2)"]) # en nA/m^2
    
    # Parámetros constructivos
    fixed_value = Vol
    P0 =101.3 # (kPa) Presión del aire a nivel de mar
    T0= 298.15  # (Kelvin) Temperatura de 25°C  + 273.15
    delta = Pres*T0/(P0*Tem) # () densidad del aire
    g0 = 29.8 # kV/cm
    Evv = grad_sup(g0, m, delta, Req*100) # kV/cm
    evv = Vol_crit(Evv, y_coor*100, Req*100) # kV
    print(f"potencial critico corona es {evv} kV y gradiente superficial critico es {Evv} kV/cm")
    epsilon0 = (1 / (36 * np.pi)) * 10**(-9)  # (F/m) permitividad del vacío
    K = 1 / (2 * np.pi * epsilon0)  # factor de multiplicación
    TolDev = 1e-2 # Tolerancia de convergencia de Jp
    #TolDev =  float(params["Tol Jp"])
    Tolerancia = [1-TolDev, 1+TolDev] # Rango de tolerancia de Jp
    visualizacion = 15 #  número de la ventana donde se muestra la evolución de la densidad de carga
    #in_condct = 'si' # Condición de que se ubica la densidad de carga alrededor del conductor o solo en un nodo central
    copiado = 'no' # Condición de que los valores sean copiados o no al momento de mantener las condiciones de borde en las mallas
    histl = False # Si se desea crear la lista de evolución de la densidad de carga
    mostrar = False # Permite mostrar los datos de la convergencia de rho inicial
    guarda = params["guardar"]
    print(f"guardado: {guarda}")
    in_condct = params["Interior conductor"]
    print(f"enconductor es: {in_condct}")
    # Definir coordenadas de la malla en x y y
    x,y, nodosx, nodosy, Sx, Sy = discretiza(Req, nodox=nodosx, nodoy=nodosy, sx=Sx, sy=Sy)  
    X,Y =np.meshgrid(x,y)
    dx = np.abs(x[1] - x[0]) # (m) distancia física entre nodos en cada columna
    dy = np.abs(y[1] - y[0]) # (m) distancia física entre nodos en cada fila
    posx_conductor, posy_conductor = encuentra_nodos(x, y, x_coor, y_coor)
    fixed_point = (posy_conductor, posx_conductor)
    print(f"malla X {X}")
    print("________________")
    print(f"malla Y {Y}")
    print(f"nodos x: {nodosx}, nodos y: {nodosy}")
    print(f"Sx: {Sx}, Sy: {Sy}")
    print(f"dx: {dx}, dy: {dy}")

    # Ajuste de viento
    wndx1 = np.ones((nodosy, nodosx)) * viento_x
    wndy1 = np.ones((nodosy, nodosx)) * viento_y
    windx, windy = windDist(wndx1, wndy1, l, Y, rug, uni=modo) # Se calcula la distribución de viento
    print(f"viento x:  {windx}")
    print("__________")
    print(f"viento y:  {windy}")

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
                P[i][j] = K*np.log(2*h[i]/(Req))
            else:
                if largo !=1 :
                    P[i][j] = K*np.log(D_im[i][j]/D[i][j])
    #V = [Vol,Vol] # voltajes sobre los conductores reales
    V = [Vol]
    Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores

    Campo_ini, Vmi, rho_n, Vm, Vol_def, Campo_fin, Ei, Ji, Jave = ejecutar_algoritmo(fixed_point, fixed_value, evv*1000, X, Y, Req, Q, Sx, mov, m, delta, g0, nodosy,
                                                                            nodosx, posx_conductor, posy_conductor, y_coor, dx, dy, windx, windy, max_iter_rho, Evv*10**5,
                                                                                max_iter, it_global, l, visualizacion, Jp_inicial=Jp_inicial,
                                                                                tolerancia=Tolerancia, condct=in_condct, copy=copiado, Muestra=mostrar, is_bundled=is_bundled)

    Exxini, Eyyini, Em = Campo_ini
    Edefx, Edefy, Edef = Campo_fin
    
    ##########################
    ######## GRAFICOS
    #########################
    # Lista gráficos
    def guarda_graficos(nombre, ruta, guarda=False):
        if guarda:
            # Asegúrate de que la carpeta existe
            if not os.path.exists(ruta):
                os.makedirs(ruta)
            # Construye la ruta de forma segura
            ruta_carpeta = os.path.join(ruta, f"{nombre}.png")
            print(ruta_carpeta)
            try:
                plt.savefig(ruta_carpeta)
                print(f"Imagen guardada en: {ruta_carpeta}")
            except PermissionError:
                print(f"No se pudo guardar el archivo en {ruta_carpeta}. Verifica permisos.")
            except Exception as e:
                print(f"Error al guardar la imagen: {e}")

    def grafE(num, ruta, mostrar=False, guarda=False):
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
        guarda_graficos("Campo_electrostatico", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()
       

    def grafV(num, ruta, mostrar=False, guarda=False):
        plt.figure(num)
        #plt.figure(figsize=(6, 6))
        #plt.contourf(X, Y, Vm, levels=200, cmap='plasma')
        c = plt.pcolormesh(X, Y, Vmi/1000, cmap='plasma', shading='auto', norm=LogNorm(vmax=np.max(Vmi)/1000))
        # Barra de color
        cbar = plt.colorbar(c)
        cbar.set_label(r'Potencial $kV$')
        ticks = cbar.get_ticks()
        # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{tick:.1f}' for tick in ticks]) 
        plt.xlabel('Distancia horizontal (m)',fontsize=11)
        plt.ylabel('Distancia vertical (m)',fontsize=11)
        plt.title('Potencial electrostático', fontsize=15)
        plt.tight_layout()
        guarda_graficos("Potencial_electrostatico", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()
       

    def grafRho(num, ruta, mostrar=False, guarda=False):
        plt.figure(num)
        
        # Verificar si rho_n es una matriz de solo ceros
        if not np.all(rho_n == 0):  # Si no todos los valores son cero
            plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto', norm=LogNorm())
            # Añadir una barra de colores para mostrar la escala
            cbar = plt.colorbar()
            cbar.set_label(r'Densidad de carga $C/m^3$', fontsize=11)
        else:
            plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto')
        
        plt.xlabel('Distancia horizontal (m)', fontsize=11)
        plt.ylabel('Distancia vertical (m)', fontsize=11)
        plt.title('Densidad de carga final', fontsize=15)
        plt.tight_layout()
        
        guarda_graficos("Densidad_carga", ruta, guarda=guarda)
        
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()

        

    def grafVf(nm, ruta, mostrar=False, guarda=False):
        plt.figure(nm)
        
        # Verificar si Vm es una matriz de solo ceros
        if not np.all(Vm == 0):  # Si no todos los valores son cero
            plt.pcolormesh(X, Y, Vm/1000, cmap='plasma', shading='auto', norm=LogNorm())
            # Añadir una barra de colores para mostrar la escala
            cbar = plt.colorbar()
            cbar.set_label(r'Potencial iónico $kV$')
        else:
            plt.pcolormesh(X, Y, Vm/1000, cmap='viridis', shading='auto')
        
        plt.title('Potencial iónico', fontsize=15)
        plt.xlabel('Distancia horizontal (m)', fontsize=11)
        plt.ylabel('Distancia vertical (m)', fontsize=11)
        
        # Si se dibujó la barra de colores, ajustar los ticks
        if np.any(Vm != 0):
            ticks = cbar.get_ticks()
            # Cambiar las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
            cbar.set_ticks(ticks)
            cbar.set_ticklabels([f'{tick:.1f}' for tick in ticks])
        
        plt.tight_layout()
        guarda_graficos("Potencial_ionico", ruta, guarda=guarda)
        
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()

        

    def grafVdef(num, ruta, mostrar=False, guarda=False):
        plt.figure(num)
        #plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
        plt.pcolormesh(X, Y, Vol_def/1000, cmap='plasma', shading='auto',norm=LogNorm())
        plt.title('Potencial definitivo',fontsize=15) 
        plt.xlabel('Distancia horizontal (m)',fontsize=11)
        plt.ylabel('Distancia vertical (m)',fontsize=11)
        # Añadir una barra de colores para mostrar la escala
        cbar = plt.colorbar()
        cbar.set_label(r'Potencial definitivo $kV$')
        ticks = cbar.get_ticks()
        # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{tick:.1f}' for tick in ticks]) 
        plt.tight_layout()
        guarda_graficos("Potencial_definitivo", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()
        

    def grafEf(nm, ruta, mostrar=False, guarda=False):
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
        guarda_graficos("Campo_definitivo", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()


    def grafJ1(num, ruta, mostrar=False, guarda=False):
        plt.figure(num)
        plt.plot(x[30:-30], Ji[30:-30]*(10**9))
        plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
        plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
        plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m, $w_x=$'+str(viento_x), fontsize=13)
        plt.tight_layout()
        plt.legend([f'$J_p$ = {str(formatear_numero(np.round(Jave*(10**9),3)))} $nA/m^2$'])
        plt.grid(True)
        guarda_graficos("Corriente_nivel_piso", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()

    def grafE1(num, ruta, mostrar=False, guarda=False):
        plt.figure(num)
        plt.plot(x[30:-30], Ei[30:-30]/1000)
        plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
        plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
        plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(viento_x), fontsize=13)
        plt.tight_layout()
        plt.legend([f'$|E|_a$ = {str(formatear_numero(np.round(np.mean(Ei/1000),3)))} kV'])
        plt.grid(True)
        guarda_graficos("Campo_nivel_piso", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()

    def grafSP(num, ruta, mostrar=False, guarda=False):
        fig = plt.figure(11)
        #fig = plt.figure(figsize=(12, 12))  # Tamaño ajustado para los subgráficos
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
        guarda_graficos("Graficos_3d", ruta, guarda=guarda)
        # Mostrar la figura si mostrar=True
        if mostrar:
            plt.show(block=False)  # Muestra la figura sin detener la ejecución
        else:
            plt.close()
        
    
    # Crear una ventana emergente para preguntar al usuario
    def preguntar_nuevo_modelo():
        ventana = tk.Tk()
        ventana.title("Nuevo modelo")
        ventana.geometry("300x150")
        
        label = tk.Label(ventana, text="¿Desea iniciar un nuevo modelo?", font=("Arial", 12))
        label.pack(pady=10)
        
        # Acción al presionar el botón
        def iniciar_nuevo_modelo():
            plt.close("all")  # Cierra todas las figuras de Matplotlib
            ventana.destroy()  # Cierra la ventana emergente
            print("Nuevo modelo iniciado.")  # Aquí puedes iniciar el nuevo modelo
        
        boton = tk.Button(ventana, text="Iniciar nuevo modelo", command=iniciar_nuevo_modelo, bg="lightblue")
        boton.pack(pady=10)
        
        # Permitir que la ventana se cierre sin bloquear las figuras
        ventana.protocol("WM_DELETE_WINDOW", ventana.destroy)
        ventana.mainloop()

    def show_plot(graf, ruta, guarda=False):
        """
        Controla el guardado y la muestra de gráficos según las opciones.
        
        Args:
            graf (list): Lista de claves que se desean mostrar.
            ruta (str): Ruta donde se guardan las imágenes.
            guarda (bool): Si True, se guardan todas las figuras.
        """
        # Diccionario con las funciones de generación de gráficos
        graficos = {
            'Eele': lambda: grafE(1, ruta, mostrar='Eele' in graf, guarda=guarda),
            'Vele': lambda: grafV(2, ruta, mostrar='Vele' in graf, guarda=guarda),
            'Rhof': lambda: grafRho(3, ruta, mostrar='Rhof' in graf, guarda=guarda),
            'Vf': lambda: grafVf(4, ruta, mostrar='Vf' in graf, guarda=guarda),
            'Vdef': lambda: grafVdef(5, ruta, mostrar='Vdef' in graf, guarda=guarda),
            'Ef': lambda: grafEf(6, ruta, mostrar='Ef' in graf, guarda=guarda),
            'J1': lambda: grafJ1(7, ruta, mostrar='J1' in graf, guarda=guarda),
            'E1': lambda: grafE1(8, ruta, mostrar='E1' in graf, guarda=guarda),
            'SPv': lambda: grafSP(9, ruta, mostrar='SPv' in graf, guarda=guarda),
        }

        # Iterar sobre todas las claves y ejecutar las funciones
        for key, func in graficos.items():
            func()  # Genera y guarda/muestra según corresponda
        plt.show(block=False)  # No bloquea la ejecución
        # Llamar a la ventana emergente
        preguntar_nuevo_modelo()


    carpeta = f"modelo_{Vol/1000}_{cantidad}_{y_coor}_{nodosx}_{nodosy}"
    ruta_destino = f"C:\\Users\\HITES\\Desktop\\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\{carpeta}"
    #print(ruta_destino)
    guarda_en_carpeta(Vol, y_coor, nodosx, nodosy, x, Ei, Ji, ruta_destino, guarda=guarda)
    show_plot(gra, ruta_destino, guarda=guarda)
    
 