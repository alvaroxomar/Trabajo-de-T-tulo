import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.colorbar as cbar
from matplotlib.colors import Normalize
from scipy.interpolate import griddata

'''
Implementa  un algoritmoque primero resuelva en toda la malla la distribucion de las densidad de carga de uno de los polos
y luego, usando esa distribución, resolver para la densidad del otro polo.


AUN QUEDA PENDIENTE DECIDIR CUAL METODO USAR PARA HACER CONVERGER LAS DENSIDADESDE CARGA
'''

plt.close('all')


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

def discretiza(min, x_coor, y_coor, nodox=None, nodoy=None, sx=None, sy=None):
    if (sx is None and sy is None) and (nodox is not None and nodoy is not None):
        Sx, Sy = espacio(min, min, nodox, nodoy)
        x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=nodox, nodoy=nodoy)
        distx=np.abs(x_coor)-Sx
        nodos_sx = int(distx/min) + 1
        disty=np.abs(y_coor)-Sy
        nodos_sy = int(disty/min) + 1
        assert np.abs(x_coor) < Sx, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sx+nodosx}, o bien ubica a menor distancia x_coor, la dif x-Sx es {distx}'
        assert np.abs(y_coor) < Sy, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sy+nodosy}, o bien ubica a menor distancia y_coor, la dif x-Sx es {disty}'
    elif (sx is not None and sy is not None) and (nodox is None and nodoy is None):
        Sx, Sy = sx, sy
        x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=None, nodoy=None)
    else:
        print('Los parámetros no están bien ajustados en la función')
    return x, y, nodosx, nodosy, Sx, Sy






def mod(z1,z2):
    return np.sqrt((z1[0]-z2[0])**2 + (z1[1]-z2[1])**2)


def calcular_campo_electrico_inicial(nodosx, nodosy, x, y, w, h, Q, K):
    # Inicializar matrices para las componentes y magnitud del campo eléctrico
    E = np.zeros((nodosy, nodosx))
    Exx = np.zeros((nodosy, nodosx))
    Eyy = np.zeros((nodosy, nodosx))
    # Calcular el campo eléctrico en cada punto de la malla
    for i in range(nodosy):
        for j in range(nodosx):
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




print('Campo eléctrico  electrostático libre iones calculado')


## Resolución para densidad de carga en el espacio
'''
def calcular_potencial_con_borde_vectorizado(E_x, E_y, dx, dy):
    """
    Calcula el potencial eléctrico V en una malla 2D a partir de las componentes E_x y E_y del campo eléctrico,
    estableciendo el potencial como nulo en todo el borde inferior, de forma vectorizada.
    
    Parámetros:
    - E_x: matriz 2D de las componentes del campo eléctrico en la dirección x.
    - E_y: matriz 2D de las componentes del campo eléctrico en la dirección y.
    - dx: paso en la dirección x de la malla.
    - dy: paso en la dirección y de la malla.
    
    Retorna:
    - V: matriz 2D del potencial eléctrico en cada punto de la malla.
    """
    # Tamaño de la malla
    nodosx, nodosy = E_x.shape
    
    # Inicializar el potencial V con ceros en el borde inferior
    V = np.zeros((nodosx, nodosy))
    
    # Integración acumulada en la dirección x (horizontal)
    # Empezamos desde el borde izquierdo en cada fila, acumulando el cambio debido a E_x
    V[:, 1:] = -np.cumsum(E_x[:, :-1] * dx, axis=1)
    
    # Integración acumulada en la dirección y (vertical) desde el borde inferior
    # Acumulamos el cambio debido a E_y en cada columna
    V[1:, :] += -np.cumsum(E_y[:-1, :] * dy, axis=0)
    
    return V
'''
# Función personalizada para formatear números
def formatear_numero(valor):
    if np.abs(valor) < 1e-2 or np.abs(valor) > 1e3:
        return f"{valor:.3e}"  # Notación científica para valores extremos
    else:
        return f"{valor:.4f}"  # Notación decimal con 4 decimales

#def rhoA(Sx, Jp, b, a,  m_ind, d,I2l):
#rho_ip = (I2l)*10**(-3)/(np.pi*b[0]*(100*a[0])*m_ind[0]*(30*d + 9*np.sqrt(d/(100*a[0])))) # radio está en cm
    #rho_in = -(Sx*Jp-I2l/2)*10**(-3)/(np.pi*b[1]*(100*a[1])*m_ind[1]*(30*d + 9*np.sqrt(d/(100*a[1])))) # radio está en cm


def visualizar_evolucion(X, Y, rho1_historial, list_dif, titulo_base="Evolución densidad de carga", figura_num=15, pausa=0.005):
    """
    Función para animar la evolución de una matriz a lo largo de iteraciones.

    Parámetros:
    - X, Y: Coordenadas de la malla.
    - rho1_historial: Lista de matrices que representan la evolución de la densidad de carga.
    - titulo_base: Título base para el gráfico.
    - pausa: Tiempo de pausa entre cada cuadro de la animación (en segundos).
    """
    #assert len(rho1_historial)>=0, 'No se tiene una lista de mapas para graficar'
    if rho1_historial is not False:
        plt.figure(figura_num)  # Especificar el número de la figura
        fig, axis = plt.subplots(num=figura_num)  # Asegura que se use esa figura
        # Configurar el primer frame
        vmin, vmax = -np.max(np.real(rho1_historial)), np.max(np.real(rho1_historial))
        pcm = axis.pcolormesh(X, Y, np.real(rho1_historial[0]), cmap='coolwarm', shading='auto',
                            norm=SymLogNorm(linthresh=1e-10, linscale=1, vmin=vmin, vmax=vmax))
        #pcm = axis.pcolormesh(X, Y, np.real(rho1_historial[0]), cmap='viridis', shading='auto', norm=LogNorm())
        cbar = plt.colorbar(pcm, ax=axis)
        cbar.set_label(r'Densidad de carga $C/m^3$', fontsize=11)
        for i, rho in enumerate(rho1_historial):
            pcm.set_array(np.real(rho).ravel())
            if i==0:
                axis.set_title(f"Evolución densidad de {titulo_base} - Iteración: {i + 1}")
            else:
                axis.set_title(f"Evolución densidad de {titulo_base} - Iteración: {i + 1}, dif_absoluta: {formatear_numero(list_dif[i-1])}")
                # Iterar sobre el historial para actualizar la animación
            plt.draw()
            plt.pause(pausa)
        # Suponiendo que tienes una lista diff_list
        plt.figure(figura_num+1)
        plt.plot(list_dif, label="Diferencia máxima por iteración")
        #plt.axhline(y=tol, color='r', linestyle='--', label="Tolerancia")
        plt.xlabel("Iteración")
        plt.ylabel("Diferencia absoluta (norma infinita)")
        plt.title(f"Evolución de la diferencia absoluta densidad de {titulo_base}")
        plt.legend()
        plt.grid(True)
        plt.show()


def rhoA(Ey,D, ep0, V0_p,V0_s, V, Ecrit_p,Ecrit_s, R):    
    '''
    Ey: Campo  de carga libre en medio camino entre los dos conductores
    Ecrit_p,_s:  Campo gradiente corona para conductor positivo y negativo
    V0_p,_s: Voltajes de gradiente corona  en los conductores positivosy negativos, estos incluyen al factor de rugosidad
    '''
    rho_ip = Ey*D*8*ep0*V0_p*(V[0]-V0_p)/(Ecrit_p*R[0]*D**2*V[0]*(5-4*V0_p/V[0]))
    rho_in = -Ey*D*8*ep0*V0_s*(V[1]-V0_s)/(Ecrit_s*R[1]*D**2*V[1]*(5-4*V0_p/V[1]))
    return rho_ip, rho_in # Entrega condiciones de borde del condcutor positivo y negativo respectivamente
'''
y_p = np.zeros_like(Y)
indices_sup = (np.abs(Y[:,0] - int(Sy-y_coor))).argmin()
indices_inf = (np.abs(Y[:,0] - int(y_coor))).argmin()
y_p[:(nodosy-indices_sup),:] = Y[indices_sup:,:]
y_p[(nodosy-indices_sup):,:] = Y[(indices_inf-len(y_p[(nodosy-indices_sup):,0])):indices_inf,:]
'''
#y_f = np.zeros_like(Y)
#y_f[indices_inf:,:] = Y[(indices_inf-(nodosy-2)):indices_inf,:]
#Ya = (np.flip(Y)-(Sy-y_coor))
def potencial_electrostático(f_point, f_value, X, Y, radio, ind, carga=None, st = 'noV', copiado='no'):
    Vm = np.zeros((nodosy,nodosx))
    Vm2 = Vm.copy()
    ## Condiciones de borde potencial inicial
    Vm[-1,:] = 0 # borde inferior
    for g in range(len(f_point)):
        py, px = f_point[g] # ubica el voltaje en el nodo de posición central del conductor
        Vm[py,px] = f_value[g]
        Vm2[f_point[g]] = f_value[g]
        Vm2[-1,:] = 0 # borde inferior
        Vm = impone_BC(Vm, Vm, f_value[g], px, py, in_condct=st, copia=copiado) # por defecto, ubica al voltaje en 5 nodos,uno central y 4 que lo  rodean
    if carga is not None:
        #print(carga)
        for z in range(len(carga)):
            # Calcular la distancia entre cada punto de la malla y el punto (w, h)
            Mod_coor = np.sqrt((X - w[z])**2 + (Y - h[z])**2)
            Mod_coor2 = np.sqrt((X - w[z])**2 + (Y + h[z])**2)
            #print(Mod_coor)
            # Calcular el potencial eléctrico en cada nodo de la malla
            #Vm += carga[z] * 0.5 * K*(1/Mod_coor - 1/Mod_coor2)
            Vm += carga[z]*K*(np.log(20/Mod_coor) - np.log(20/Mod_coor2))
            #Vm = calcular_potencial_con_borde_vectorizado(Exx, Eyy,dx, dy)
    else:
        c1 = np.cosh(np.pi*(X-2*ind*Sx)/(2*(Sy)))
        c2 = np.cos(np.pi*Ya/(2*(Sy)))
        c22 = np.cos(np.pi*y_p/(2*(Sy)))
        c3 = np.cosh(np.pi*ind*Sx/(Sy))
        c4 = np.cos(np.pi*radio/(2*(Sy)))
        Vm = np.abs(f_value * np.log((c1-c2)/(c1+c2))/np.log((c3-c4)/(c3+c4)))
        Vm2 = np.abs(f_value * np.log((c1-c22)/(c1+c22))/np.log((c3-c4)/(c3+c4)))
    for g in range(len(f_point)):
        py, px= f_point[g]
        Vm[f_point[g]] = f_value[g]
        Vm2[f_point[g]] = f_value[g]
        Vm2[-1,:] = 0 # borde inferior
        Vm = impone_BC(Vm, Vm, f_value[g], px, py, in_condct=st, copia=copiado)
    Vm[-1,:] = 0 # borde inferior
    return Vm,Vm2

def resuelve2(a, b, c, d, g, f):
    A2 = d*(d-g*b/a)
    B2 = f*(2*d-g*b/a)+2*c*g**2
    C2 = f**2
    Gr = (-B2 + np.sqrt(B2**2-4*A2*C2+0j))/(2*A2)
    ysol = np.sqrt(Gr)
    Gy = (b*ysol)**2 - 4*a*c
    xsol = (-b*ysol + np.sqrt(Gy+0j))/(2*a)
    return xsol,ysol


'''
# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
def update_rho_vectorized(rhoinip, rhoinin, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_boundp, rho_boundn, met=0):
    # Aquí hace falta copiar los rhop y rhon por lo que debo integrarlos en los parámetros de estas función.
    #rho_new = np.copy(rho)
    rhop = rhoinip.copy()
    rhon = rhoinin.copy()
    rhop = rhop.astype(complex)
    rhon = rhon.astype(complex)
    rhop[fixed_point[0]] = rho_boundp[fixed_point[0]]  # Mantener la condición en el o los puntos centrales
    rhop[-1, :] = rho_boundp[-1,:]  # Mantener la condición en el borde inferior
    rhon[fixed_point[1]] = rho_boundn[fixed_point[1]]  # Mantener la condición en el o los puntos centrales
    rhon[-1, :] = rho_boundn[-1,:]  # Mantener la condición en el borde inferior
    # Calcular las derivadas en el interior (diferencias centrales) usando slicing
    d_rho_dxp = (rhop[1:-1, 2:] - rhop[1:-1, :-2]) / (2 * dx)
    d_rho_dyp = (rhop[2:, 1:-1] - rhop[:-2, 1:-1]) / (2 * dy)
    d_rho_dxn = (rhon[1:-1, 2:] - rhon[1:-1, :-2]) / (2 * dx)
    d_rho_dyn = (rhon[2:, 1:-1] - rhon[:-2, 1:-1]) / (2 * dy)
    Ewxp = Ex*movp + wx
    Ewyp = Ey*movp + wy
    Ewxn = Ex*movn - wx
    Ewyn = Ey*movn - wy
    a = movp/epsilon0
    b = Rep/ele - movp/epsilon0
    c = (d_rho_dxp*Ewxp[1:-1,1:-1] + d_rho_dyp*Ewyp[1:-1,1:-1])
    d = movn/epsilon0
    g = Rep/ele - movn/epsilon0
    f = -(d_rho_dxn*Ewxn[1:-1,1:-1] + d_rho_dyn*Ewyn[1:-1,1:-1])
    if met==0:
        rp, rn = resuelve2(a,b,c,d,g,f)
        rhop[1:-1,1:-1] = rp
        rhon[1:-1,1:-1] = rn
    elif met==1:
        rp_p, rp_n = resuelve1(a, b , c)
    # Ecuación diferencial discretizada en el interior
    G = -epsilon0 * (Ewx[1:-1, 1:-1] * d_rho_dx + Ewy[1:-1, 1:-1] * d_rho_dy)
    D1 = np.sqrt(np.abs(G))
    rho[1:-1, 1:-1] = np.where(G >= 0, D1, np.sqrt(np.abs(G)))
    #rho[1:-1, 1:-1] = np.where(G >= 0, np.sqrt(G), 0)
    
    # Condiciones en los bordes
    Gi = -epsilon0 * (Ewx[1:-1, 0] * (rho[1:-1, 1] - rho[1:-1, 0]) / dx + Ewy[1:-1, 0] * (rho[2:, 0] - rho[:-2, 0]) / (2 * dy))
    Gd = -epsilon0 * (Ewx[1:-1, -1] * (rho[1:-1, -1] - rho[1:-1, -2]) / dx + Ewy[1:-1, -1] * (rho[2:, -1] - rho[:-2, -1]) / (2 * dy))
    Gs = -epsilon0 * (Ewx[0, 1:-1] * (rho[0, 2:] - rho[0, :-2]) / (2 * dx) + Ewy[0, 1:-1] * (rho[1, 1:-1] - rho[0, 1:-1]) / dy)
    Geis = -epsilon0 * (Ewx[0, 0] * (rho[0, 1] - rho[0, 0]) / dx + Ewy[0, 0] * (rho[1, 0] - rho[0, 0]) / dy)
    Geds = -epsilon0 * (Ewx[0, -1] * (rho[0, -1] - rho[0, -2]) / dx + Ewy[0, -1] * (rho[1, -1] - rho[0, -1]) / dy)
    Geii = -epsilon0 * (Ewx[-1, 0] * (rho[-1, 1] - rho[-1, 0]) / dx + Ewy[-1, 0] * (rho[-1, 0] - rho[-2, 0]) / dy)
    Gedi = -epsilon0 * (Ewx[-1, -1] * (rho[-1, -1] - rho[-1, -2]) / dx + Ewy[-1, -1] * (rho[-1, -1] - rho[-2, -1]) / dy)
    D2 = np.sqrt(np.abs(Gi))
    D3 = np.sqrt(np.abs(Gd))
    D4 = np.sqrt(np.abs(Gs))
    D5 = np.sqrt(np.abs(Geis))
    D6 = np.sqrt(np.abs(Geds))
    D7 = np.sqrt(np.abs(Geii))
    D8 = np.sqrt(np.abs(Gedi))
    rho[1:-1, 0] = np.where(Gi >= 0, D2, np.sqrt(np.abs(Gi)))
    #rho[1:-1, 0] = np.where(Gi >= 0, np.sqrt(Gi), 0) # Borde izquierdo
    rho[1:-1, -1] = np.where(Gd >= 0, D3, np.sqrt(np.abs(Gd)))
    #rho[1:-1, -1] = np.where(Gd >= 0, np.sqrt(Gd), 0) # Borde derecho
    rho[0, 1:-1] = np.where(Gs >= 0, D4, np.sqrt(np.abs(Gs)))
    #rho[0, 1:-1] = np.where(Gs >= 0, np.sqrt(Gs), 0) # Borde superior
    rho[0, 0] = np.where(Geis >= 0, D5, np.sqrt(np.abs(Geis))) # Esquina superior izquierda
    rho[0, -1] = np.where(Geds >= 0, D6, np.sqrt(np.abs(Geds))) # Esquina superior derecha
    rho[-1, 0] = np.where(Geis >= 0, D7, np.sqrt(np.abs(Geii)))  # Esquina inferior izquierda
    rho[-1, -1] = np.where(Geds >= 0, D8, np.sqrt(np.abs(Gedi)))  # Esquina inferior derecha

    # Aplicar condiciones de borde
    rhop[fixed_point[0]] = rho_boundp[fixed_point[0]]  # Mantener la condición en el o los puntos centrales
    rhop[-1, :] = rho_boundp[-1,:]  # Mantener la condición en el borde inferior
    rhon[fixed_point[1]] = rho_boundn[fixed_point[1]]  # Mantener la condición en el o los puntos centrales
    rhon[-1, :] = rho_boundn[-1,:]  # Mantener la condición en el borde inferior
    return rhop, rhon
'''

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

'''
def convergencia(ro1,ro0, tol):
    # Define un pequeño valor para evitar divisiones por cero
    epsilon = 1e-10
    epsilon_array = np.full_like(ro1, epsilon)
    # Calcula la diferencia relativa usando el máximo entre |rho1|, |rho_0| y epsilon
    diferencia_relativa = np.abs(ro1 - ro0) / np.maximum.reduce([np.abs(ro1), np.abs(ro0), epsilon_array])
    # Verifica si todas las diferencias relativas son menores al 10%
    condicion = np.all(diferencia_relativa < tol)
    diff = np.linalg.norm(diferencia_relativa,ord=np.inf)
    return condicion, diff
'''
def convergencia(rho_actual, rho_anterior, tol):
    # Calcula la diferencia absoluta
    diferencia_absoluta = np.abs(rho_actual - rho_anterior)
    # Norma infinita de la diferencia absoluta
    max_diferencia = np.linalg.norm(diferencia_absoluta, ord=np.inf)
    # Verifica si la diferencia máxima es menor a la tolerancia
    condicion = max_diferencia < tol
    return condicion, max_diferencia

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

def evaluar_estado(dist_rho1, rhop0, diff_list, tol, px, py, ventana=110, umbral_variacion=1e-6, rho_hist=False):
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
# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
'''
def algoritmo_rho_v(V, rho_ini, dx, dy, windx, windy,max_iter_rho, Jplatep, Jplaten, rho_A):
    conv = 0
    rho_bp =np.zeros((nodosx,nodosy), dtype=complex)
    rho_bn =np.zeros((nodosx,nodosy), dtype=complex)
    # En un comienzo se tiene el campo eléctrico electrostático
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    # Se define a la malla rho_b con las condiciones de borde
    rho_bp[fixed_point[0]] = rho_A[0]
    rho_bn[fixed_point[1]] = rho_A[1]
    # COndición de desplazamiento de la densidad de carga debido al viento
    rho_bp[-1,:] = Jplatep/(movp*np.sqrt(np.abs((Eyyi[-1,:]+windy[-1,:]/movp)**2))) # se postula que la dirección de E será solamente vertical
    rho_bn[-1,:] = Jplaten/(movn*np.sqrt(np.abs((Eyyi[-1,:]+windy[-1,:]/movn)**2)))
    rho1 = rho_ini.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
    rho1 = rho1.astype(complex)
    for iteration in range(max_iter_rho):
        rho0 = rho1.copy()
        if iteration == 0:
            rho1p, rho1n = update_rho_vectorized(rho_bp, rho_bn, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_bp, rho_bn)
        else:
            rho1p, rho1n = update_rho_vectorized(rho1p, rho1n, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_bp, rho_bn)
        rho1 = rho1p-rho1n # rho1py rho1n son usados en sus valores absolutos
        # Criterio de convergencia
        condicion, diff = convergencia(rho1, rho0, 0.01)
        #print(r'Tolerancia $\rho_1$ = '+str(diff))
        if condicion:
            #print(f"Convergencia para rho alcanzada en la iteración {iteration}")
            break
    #print(r'Tolerancia rho_1 = '+str(diff))
    return np.real(rho1), np.real(rho1p), np.real(rho1n)
'''
################################################################
################################################################
def algoritmo_rho(V, rho_inip, rho_inin, rho_b, max_iter_rho, dx, dy, movp, movn, wx,wy, val_cond, pos, condi='si', copia='no',
                   rho_h=False, visu=15, muestra=False, fi=30):
    '''
    Algoritmo principal de resolución del sistema de ecuaciones de la densidad de cara positiva y negativa
    Comienza con la fijacion de la distribucion de rho- como valores nulos excepto en el conductor negativo
    Actuliza los valores de rho+ y una vez resuelto, actuliza los valores de rho- con rho+ fijo
    Entrega las partes reales de rho+, rho- y rho+-rho- 
    V: potencial existente, se usa para el cálculo del campo eleçtrico
    rho_ini: matriz con las dimensiones (ny, nx), sirve de referencia para crear rho+ y rho-
    val_cond: arreglo con las condiciones de borde de conductor positivo y negativo
    pos: arreglo de las posiciones de los condcutores
    '''
    val1 = val_cond[0]
    val2 = val_cond[1]
    py1, px1 = pos[0]
    py2, px2 = pos[1]
    rhop_0 = rho_inip.copy() # copia la  red con los valores de rho_p que ya existe de antes
    rhop_0= rhop_0.astype(complex)
    rhon_0 = rho_inin.copy() # copia la  red con los valores de rho_n que ya existe de antes
    rhon_0= rhon_0.astype(complex)
    rhop_0 =impone_BC(rhop_0, rho_b, val1, px1, py1, in_condct=condi, copia=copia) # rhop_0 sería una red con condición de borde unicamente en la posición del conductor
    rhon_0 =impone_BC(rhon_0, rho_b, val2, px2, py2, in_condct=condi, copia=copia) # rhon_0 sería una red con condición de borde unicamente en la posición del conductor
    print(f'el máximo de rho_p es: {rhop_0[(py1,px1)]}')
    '''
    plt.figure(50)
    mesh = plt.pcolormesh(X, Y, np.real(rhop_0))
    cbar = plt.colorbar(mesh)
    cbar.set_label('Densidad de carga pos ini')
    plt.draw()  # Dibuja el gráfico
    plt.figure(51)
    mesh = plt.pcolormesh(X, Y, np.real(rhon_0))
    cbar = plt.colorbar(mesh)
    cbar.set_label('Densidad de carga neg ini')
    plt.draw()  # Dibuja el gráfico
    plt.show()
    input("Presiona Enter para continuar...")
    '''
    print(f'el mínimo de rho_n es: {rhon_0[(py2,px2)]}')
    if rho_h is not False:
        rhop_hist = [rhop_0.copy()]
        rhon_hist = [rhon_0.copy()]
    else:
        rhop_hist = False
        rhon_hist = False
    Ex, Ey, Em = calcular_campo_electrico(V, dx, dy)
    Ewxp = Ex*movp + wx
    Ewyp = Ey*movp + wy
    Ewxn = Ex*movn - wx
    Ewyn = Ey*movn - wy
    a = movp/epsilon0
    #a = 1/epsilon0
    b = Rep/ele - movp/epsilon0
    d = movn/epsilon0
    #d = 1/epsilon0
    g = Rep/ele - movn/epsilon0
    rhop_updated, rphist, dif_listp = iteraciones_rho(max_iter_rho, dx, dy, Ewxp, Ewyp, rhop_0, rhon_0, rho_b, val1, px1, py1, a, b, fac_c=1,
                                            sig=1, condi=condi, copia=copia, rho_hist=rhop_hist, nombre='rho_p', est_gra=muestra, fi=fi)
    ext_rhop = extrapolate_borders(X, Y, rhop_updated.copy(), num_layers=30)
    rhon_updated, rnhist, dif_listn = iteraciones_rho(max_iter_rho, dx, dy, Ewxn, Ewyn, rhon_0, rhop_0, rho_b, val2, px2, py2, d, g, fac_c=-1,
                                            sig=-1, condi=condi, copia=copia, rho_hist=rhon_hist, nombre='rho_n', est_gra=muestra, fi=fi+1)
    ext_rhon = extrapolate_borders(X, Y, rhon_updated.copy(), num_layers=30)
    '''
    rhop_updated1, rphist, dif_listp = iteraciones_rho(max_iter_rho, dx, dy, Ewxp, Ewyp, rhop_updated, rhon_updated, rho_b, val1, px1, py1, a, b, fac_c=1,
                                            sig=1, condi=condi, copia=copia, rho_hist=rhop_hist, nombre='rho_p', est_gra=muestra, fi=fi)
    ext_rhop = extrapolate_borders(X, Y, rhop_updated1.copy(), num_layers=30)
    rhon_updated1, rnhist, dif_listn = iteraciones_rho(max_iter_rho, dx, dy, Ewxn, Ewyn, rhon_updated, rhop_updated, rho_b, val2, px2, py2, d, g, fac_c=-1,
                                            sig=-1, condi=condi, copia=copia, rho_hist=rhon_hist, nombre='rho_n', est_gra=muestra, fi=fi+1)
    ext_rhon = extrapolate_borders(X, Y, rhon_updated1.copy(), num_layers=30)
    #rhop_updated, rphist, dif_listp = iteraciones_rho(max_iter_rho, dx, dy, Ewxp, Ewyp, ext_rhop, ext_rhon, rho_b, val1, px1, py1, a, b, fac_c=1,
    #                                        sig=1, condi=condi, copia=copia, rho_hist=rhop_hist, nombre='rho_p', est_gra=muestra, fi=fi)
    #ext_rhop = extrapolate_borders(X, Y, rhop_updated.copy(), num_layers=30)
    '''
    print(f'confirmo maximo rhop es  {rhop_updated[(py1,px1)]}')
    print(f'confirmo minimi rhon es  {rhon_updated[(py2,px2)]}')
    rho_total = ext_rhop + ext_rhon
    visualizar_evolucion(X, Y, rnhist, dif_listp, titulo_base='carga negativa', figura_num=visu, pausa=0.005)
    visualizar_evolucion(X, Y, rphist, dif_listn, titulo_base='carga positiva', figura_num=visu+2, pausa=0.005)
    return np.real(ext_rhop), np.real(ext_rhon), np.real(rho_total)


def iteraciones_rho(iteraciones, dx, dy, Ewx, Ewy, rho_iterado, rho_fijo, rho_b, val, px, py, A, B, fac_c=1, sig=1, condi = 'si',
                     copia='si', rho_hist=False, nombre=None, est_gra=False, fi=30):
    '''
    Resuelve de manera iterativa la ecuación cuadrática para la densidad de carga asociado a un condcutor específico.
    Deja los bordes sin modificación.
    Entrega la red de rho actualizado con sus valores  complejos
    rho_iterado: es la red con los valores de rho que se desean actualizar iteradamente
    rho_fijo: el la red de valores de rho correspondiente al conductor opuesto el cual se mantiene fijo durante todas las iteraciones
    val: valor de la densidad de carga en el conductor iterado
    px,py: posiciones del conductor iterado
    fac_c: depende del tipo de conductor considerado
    '''
    dist_rho1 =np.zeros_like(rho_iterado)
    dist_rho1 = dist_rho1.astype(complex)
    diff_list = []
    promedios = []
    desviaciones = []
    todos = False
    for i in range(iteraciones):
        if i == 0:
            C = delta_Vp(fac_c, rho_iterado, Ewx, Ewy, dx, dy) # toma  la red que corresponde únicamente a un conductor
            rhop0 = rho_iterado.copy()
            dist_rho1[1:-1,1:-1] = resuelve1(A, B, C, rho_fijo[1:-1, 1:-1], sig=sig)
            dist_rho1 = impone_BC(dist_rho1, rho_b, val, px, py, in_condct=condi, copia=copia)
        else:
            rhop0 = dist_rho1.copy()
            C = delta_Vp(fac_c, dist_rho1, Ewx, Ewy, dx, dy)
            dist_rho1[1:-1,1:-1] = resuelve1(A, B, C, rho_fijo[1:-1, 1:-1], sig=sig)
            dist_rho1 = impone_BC(dist_rho1, rho_b, val, px, py, in_condct=condi, copia=copia)
        convergencia_total, estado = evaluar_estado(dist_rho1, rhop0, diff_list, tol=9e-07, px=px, py=py, rho_hist=rho_hist)
        #condicion,diff = convergencia(dist_rho1, rhop0,1e-6)
        if estado["promedio"] is not None and estado["dev_dif"] is not None:
            promedios.append(estado["promedio"])
            desviaciones.append(estado["dev_dif"])
        #if est_gra is not False:
            #if estado["promedio"] is not None and estado["dev_dif"] is not None:
            #    print(f'La tendencia al promedio es: {estado["tendencia"]}')
            #    print(f'Promedio actual: {estado["promedio"]}')
            #    print(f'Desviación del promedio actual: {estado["dev_dif"]}')
            #    print('....................')
        if convergencia_total:
            print(f'Diferencia relativa {nombre}: {estado["diff"]}')
            print(f'Promedio actual: {estado["promedio"]}')
            print(f'Desviación del promedio actual: {estado["dev_dif"]}')
            print(f'Término del algoritmo en iteración: {i}')
            break
        #if rho_hist is not None:
        #    rho_hist.append(dist_rho1.copy())
        #if condicion:
            #print(f'Convergencia alcanzada para {nombre} en la iteración: '+str(i))
         #   print(f'Diferencia relativa {nombre}: '+str(diff))
          #  break
        if i == iteraciones-1:
            print(f'No hay convergencia de {nombre}: {estado['diff']}')
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

def resuelve1(a, b, c, rho_fijo, sig=1):
    B = b*rho_fijo
    if sig==1:
        sol1 = (-B + np.sqrt(B**2-4*a*c+0j))/(2*a)
        return sol1
    elif sig==-1:
        sol2 = (-B - np.sqrt(B**2-4*a*c+0j))/(2*a)
        return sol2

def delta_Vp(fac, rho_grid, ewx,ewy, dx, dy):
    '''
    Calcula el coeficiente constante de la ecuación de modelamiento de la densidad de carga espacial
    fac: factor (+1 o-1) el cual es asignado segun el tipo de conductor tratado.
    Deja  los bordes sin modificación
    '''
    d_rho_dx = (rho_grid[1:-1, 2:] - rho_grid[1:-1, :-2]) / (2 * dx)
    d_rho_dy = (rho_grid[2:, 1:-1] - rho_grid[:-2, 1:-1]) / (2 * dy)
    Coef = fac*(d_rho_dx*ewx[1:-1,1:-1] + d_rho_dy*ewy[1:-1,1:-1])
    return Coef

def impone_BC(rho_red, rho_b, valor, posx, posy, in_condct='si', copia='no'):
    '''
    rho_red: Malla en el cual se impone la condicion del punto
    rho_b: Malla de referencia donde se almacena los valores de los conductores
    valor: Valor de las densidad de carga del conductor, es un único valor
    in_cond: condicional si se desea que se distribuya el valor en la periferia del conductor o bien en un único punto
    copia: condicional si se desea copiar los valores de rho_b o no
    '''
    #print(f'valor es: {valor}')
    #print(f'posx es: {posx}')
    #print(f'posy es: {posy}')
    if in_condct=='si':
        rho_red[posy, posx] = valor
    elif in_condct=='no':
        if copia=='no':
            rs, rn, re, ro = val_rhoA(valor)
        elif copia=='si':
            rs = rho_b[posy+1, posx]
            rn = rho_b[posy-1, posx]
            re = rho_b[posy, posx+1]
            ro = rho_b[posy, posx-1]
        rho_red[posy+1, posx] = rs
        rho_red[posy-1, posx] = rn
        rho_red[posy, posx+1] = re
        rho_red[posy, posx-1] = ro
    elif in_condct == 'noV':
        rho_red[posy+1, posx] = valor
        rho_red[posy-1, posx] = valor
        rho_red[posy, posx+1] = valor
        rho_red[posy, posx-1] = valor
    else:
        print('in_condct escogido es inválido')
    return rho_red

def val_rhoA(val):
    r_s=val*np.cos((np.pi-np.pi)/2)
    r_n=val*np.cos((np.pi-0)/2)
    r_e=val*np.cos((np.pi-np.pi/2)/2)
    r_o=val*np.cos((np.pi-np.pi/2)/2)
    return r_s, r_n, r_e, r_o

############################################################
############################################################
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
        u[:, 0] = u[:, 1] + dx * (Exx[:, 0])
    else:
        u[:, 0] = V_boundary[:, 0]  # Dirichlet

    # Borde derecho (Neumann si Exx está definido)
    if Exx is not None and dx is not None:
        #u[:, -1] = u[:, -2] - dx * (Exx[:, -1]+windx[:,-1]/mov)
        u[:, -1] = u[:, -2] - dx * (Exx[:, -1])
    else:
        u[:, -1] = V_boundary[:, -1]  # Dirichlet

    # Borde superior (Neumann si Eyy está definido)
    if Eyy is not None and dy is not None:
        #u[0, :] = u[1, :] + dy * (Eyy[0, :]+windy[0,:]/mov)
        u[0, :] = u[1, :] + dy * (Eyy[0, :])
    else:
        u[0, :] = V_boundary[0, :]  # Dirichlet

    # Borde inferior (siempre Dirichlet)
    u[-1, :] = V_boundary[-1, :]
    u[fixed_point[0]] = fixed_value[0]
    u[fixed_point[1]] = fixed_value[1]

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
            u = apply_boundary_conditions_2D(u, V_boundary, Exx=Exx, Eyy=Eyy, dx=dx, dy=dy)
        
        # Imponer condición fija en el punto (si está definida)
        if fixed_point is not None and fixed_value is not None:
            u[fixed_point[0]] = fixed_value[0]
            u[fixed_point[1]] = fixed_value[1]
    
    return u

# Resolver usando FMG
# Algoritmo que obtien la red de potencial eléctrico en base a los valores de Rho

def algoritmo_V_rho(V, rho1, dx, dy, fd_point, fd_value, max_iter):
    f_rhs = funcion_f(rho1)
    #V_b =  Vb.copy()
    #Vi = V.copy()
    V_b = np.zeros_like(V)
    for iteration in range(max_iter):
        Vold = V.copy()
        V = update_v(V, f_rhs, dx,dy,fixed_point=fd_point, fixed_value=[0,0], V_boundary=V_b, Exx=None,Eyy=None)
        condicion, diff = convergencia(V, Vold, 1e-6)
        if condicion:
            print(f"Convergencia alcanzada para V en la iteración {iteration}")
            break
    return V

def E_onset(m, delta, a):
    '''
    a: radio conductor (m)
    delta: densidad relativa (ad)
    m: factor de conductor (ad)
    '''
    g0 =29.8 # kV/cm
    return 100*g0*m*(delta + 0.301*np.sqrt(delta/(100*a))) # kV/m
def V_onset(m, Ecr, a):
    '''
    a: radio conductor (m)
    m: factor de conductor (ad)
    Ecr: Campo crítico corona (kV/m)
    '''
    return m*a*2*Ecr*np.log(2*y_coor1/a)
'''
Von = V_onset

Exxini, Eyyini, Em = calcular_campo_electrico(Vmi, dx, dy) # campo electrostático inicial

#rho_inicial[fixed_point[0]] = rho_ip
#rho_inicial[fixed_point[1]] = rho_in
#rho_inicial[-1,:] = Jpp/(movp*np.sqrt(np.abs((Eyyini[-1,:]+windy[-1,:]/movp)**2)))
#rho_inicial[-1,:] += Jpn/(movn*np.sqrt(np.abs((Eyyini[-1,:]-windy[-1,:]/movn)**2)))
#print(Vmi)
Vm = Vmi.copy()
it_global = 2
for n in range(it_global):
    Volder = Vm.copy()
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    if n==0:
        #rho_d, rho_p, rho_n = algoritmo_rho_v(Vmi, rho_inicial, dx, dy, windx, windy, max_iter_rho, Jpp, Jpn, [rho_ip, rho_in])
        rho_p, rho_n, rho_d = algoritmo_rho(Vmi, rho_inicial, rho_inicial, rho_boundary, max_iter_rho, dx, dy, movp, movn, windx, windy,
                                             [rho_ip, rho_in], fixed_point, condi=cn_cond, copia=copiado, rho_h=None)
    else:
        #rho_d, rho_p, rho_n = algoritmo_rho_v(Vm, rho_d0, dx, dy, windx, windy, max_iter_rho, Jpp, Jpn, [rho_ip, rho_in])
        rho_p, rho_n, rho_d = algoritmo_rho(Vmi+Vm, rho_p, rho_n, rho_boundary, max_iter_rho, dx, dy, movp, movn, windx, windy,
                                             [rho_ip, rho_in], fixed_point, condi=cn_cond, copia=copiado, rho_h=None)
    #rho_d0 = rho_p-rho_n # LOS rho_pn son puestos en sus valores absolutos
    Vm = algoritmo_V_rho(Vm, rho_d, dx, dy, fixed_point, fixed_value, max_iter, Vmi)
    condicion,diff = convergencia(Vm, Volder, 0.01)
    if condicion:
        print(f"Convergencia alcanzada para V en la iteración {n}")
        print(r'Diferencia relativa V y Vold: '+str(diff))
        break
    else:
        print(r'Diferencia relativa V y Vold: '+str(diff))
print(r'Diferencia relativa V y Vold: '+str(diff))

#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
print('Potencial calculado')
##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica

#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')
'''

def densidad_voltaje(Volt, rho, Rho_val1, Rho_val2, points, cond='si', cop='no'):
    py1, px1 = points[0]
    py2, px2 = points[1]
    '''
    print(f'pos conductor izquierdo {points[0]}')
    print(f'pos conductor derecho {points[1]}')
    '''
    if len(Volt) != 2:
        # Calcula la densidad de cargas usando la ecuación de Poisson según el potencial ya calculado
        Rho_new = np.zeros_like(rho)
        Rho_new = -epsilon0*apply_laplacian_asymmetric(Volt, dx, dy)
        #Rho_new = -epsilon0*laplace(Volt)
        Rho_new = impone_BC(Rho_new, rho, Rho_val1, px1, py1, in_condct=cond, copia=cop)
        Rho_new = impone_BC(Rho_new, rho, Rho_val2, px2, py2, in_condct=cond, copia=cop)
        Rho_new  = extrapolate_borders(X,Y, Rho_new)
        '''
        Rho_new[0,:] = rho[0,:]
        Rho_new[-1,:] = rho[-1,:]
        Rho_new[:,0] = rho[:,0]
        Rho_new[:,-1] = rho[:,-1]
        '''
        return Rho_new
    elif len(Volt)==2:
        Rho_newp = np.zeros_like(rho)
        Rho_newn = np.zeros_like(rho)
        Rho_newp = -epsilon0*apply_laplacian_asymmetric(Volt[0], dx, dy)
        Rho_newn = -epsilon0*apply_laplacian_asymmetric(Volt[1], dx, dy)
        #Rho_new = -epsilon0*laplace(Volt)
        Rho_newp = impone_BC(Rho_newp, rho, Rho_val1, px1, py1, in_condct=cond, copia=cop)
        Rho_newn = impone_BC(Rho_newn, rho, Rho_val2, px2, py2, in_condct=cond, copia=cop)
        Rho_newp  = extrapolate_borders(X,Y, Rho_newp)
        Rho_newn  = extrapolate_borders(X,Y, Rho_newn)
        return Rho_newp, Rho_newn

def calcular_potencial_inicial(fixed_point, fixed_value, X, Y, R, Q, st='noV'):
    """Calcula el potencial inicial basado en el método CSM."""
    Vmi, Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, R, 0, carga=Q, st=st)
    return Vmi

def inicializar_densidad(Ey,  diametro, R, m, delta, nodosy, nodosx, pos_conductor1, pos_conductor2, con_condct='si', copiado='no'):
    """Inicializa la densidad de carga y las condiciones de borde."""
    mp, mn = m # inidices de  rugosidad
    posy1, posx1= pos_conductor1 # en formato (y,x)
    posy2, posx2= pos_conductor2 # en formato (y,x)
    Ecritp= E_onset(mp, delta, R[0]) # cálcula gradiente superficial crítico polo positivo
    Ecritn= E_onset(mn, delta, R[1]) # cálcula gradiente superficial crítico polo negativo
    V0p = V_onset(mp, Ecritp, R[0]) # cálcula voltaje crítico polo positivo
    V0n = V_onset(mn, Ecritn, R[1]) # cálcula voltaje crítico polo negativo
    rho_ip, rho_in = rhoA(Ey, diametro, epsilon0, V0p, V0n, V, Ecritp, Ecritn, R)
    print(f'condicion rho_p: {rho_ip}')
    print(f'condición rho_n: {rho_in}')
    rho_inicial = np.zeros((nodosy, nodosx),dtype=complex) # Distribución inicial solo con valores nulos
    rho_boundary = np.zeros_like(rho_inicial)
    rho_boundary = impone_BC(rho_boundary,rho_inicial, rho_in, posx2, posy2, in_condct=con_condct, copia='no')# ubica densidad carga condc. negativo
    rho_boundary = impone_BC(rho_boundary,rho_inicial, rho_ip, posx1, posy1, in_condct=con_condct, copia='no')# ubica densidad carga condc. positivo
    # Graficar
    plt.figure(40)
    mesh = plt.pcolormesh(X, Y, np.real(rho_boundary))
    cbar = plt.colorbar(mesh)
    cbar.set_label('Densidad de carga')
    plt.draw()  # Dibuja el gráfico
    plt.show()
    #input("Presiona Enter para continuar...")
    return rho_inicial, rho_boundary, rho_ip, rho_in

def iterar_potencial(Vmi, rho_inicial, rho_boundary, rho_i, fixed_point, dx, dy, windx, windy, max_iter_rho, max_iter,
                      it_global, rho_h = False, visu=15, con_condct='si', copiado='no', must = False, fi=30):
    """Itera para calcular el potencial eléctrico y la densidad de carga."""
    Vm = np.zeros_like(Vmi)
    difer_global = []
    rho_ip, rho_in = rho_i
    '''
    print(f'rho_p  encond es: {rho_ip}')
    print(f'rho_n  encond es: {rho_in}')
    print('Campo eléctrico  electrostático libre iones calculado')
    '''
    for n in range(it_global):
        Volder = Vm.copy()
        if n == 0:
            print('Primer rho')
            rho_p, rho_n, rho_d = algoritmo_rho(Vmi, rho_inicial, rho_inicial, rho_boundary, max_iter_rho, dx, dy, movp, movn, windx, windy,
                                             [rho_ip, rho_in], fixed_point, condi=con_condct, copia=copiado, rho_h=rho_h, visu=visu, muestra=must, fi=fi)
            #Vp =algoritmo_V_rho(Vm, rho_p, dx, dy,fixed_point[0], 0, max_iter)
            #Vn =algoritmo_V_rho(Vm, rho_n, dx, dy,fixed_point[1], 0, max_iter)
        else:
            print(f'Comienza nueva actualización de rho número {n + 1}')
            #rho_p, rho_n, rho_d = algoritmo_rho(Vmi+Vm, rho_p, rho_n, rho_boundary, max_iter_rho, dx, dy, movp, movn, windx, windy,
            #                                 [rho_ip, rho_in], fixed_point, condi=con_condct, copia=copiado, rho_h=rho_h, visu=visu, muestra=must, fi=fi)
            #rho_p, rho_n = densidad_voltaje([Vp,Vn], rho_boundary, rho_ip, rho_in, fixed_point, cond=con_condct, cop=copiado)
            rho_d = densidad_voltaje(Vm, rho_boundary, rho_ip, rho_in, fixed_point, cond=con_condct, cop=copiado)
        Vm = algoritmo_V_rho(Vm, rho_d, dx, dy, fixed_point, 0, max_iter)
        #rho_d = densidad_voltaje(Vm, rho_boundary, rho_ip, rho_in, fixed_point, cond=con_condct, cop=copiado)
        condicion, max_diff = convergencia(Vm, Volder, 3.1e-9)
        if condicion:
            print(f"Convergencia alcanzada para V en la iteración {n}")
            break
        else:
            print(f'Diferencia relativa V y Vold: {max_diff}')
            print('*****************************')
        if n==it_global-1:
            print(f'se alcanzó el máximo de iteraciones: {n}')
    Vm = -Vm
    return Vm, rho_p, rho_n, rho_d


def calcular_resultados_finales(Vm, Vmi, rho_p, rho_n, dx, dy, mov, windx, windy, l):
    """Calcula el campo eléctrico, densidad de corriente y verifica la convergencia."""
    Vol_def = Vm + Vmi
    movp, movn = mov
    Exxini, Eyyini, Eini =calcular_campo_electrico(Vmi, dx, dy)
    Edefx, Edefy, Edef = calcular_campo_electrico(Vol_def, dx, dy)
    Ei = Edef[encuentra_nodos(0, l)[1],  :] # Magnitud campo eléctrico nivel de piso
    Jplus = rho_p*movp*np.sqrt((Edefx+(windx/movp))**2 + (Edefy+(windy/movp))**2)
    Jdiff = rho_n*movn*np.sqrt((Edefx-(windx/movn))**2 + (Edefy-(windy/movn))**2)
    Jtotal = np.sqrt((rho_p*movp*(Edefx+(windx/movp))+rho_n*movn*(Edefx+(windx/movn)))**2 + (rho_p*movp*(Edefy+(windy/movp))+rho_n*movn*(Edefy+(windy/movn)))**2)
    J = Jplus + Jdiff
    Jave = np.mean(J[-1,:]) # Promedio densidad de corriente a nivel  de piso
    Ji = J[encuentra_nodos(0, l)[1], :]  # Densidad de corriente a nivel de l
    Campo_ini = [Exxini, Eyyini, Eini]
    Campo_fin = [Edefx, Edefy, Edef]
    return Jave, Ji, Campo_ini, Campo_fin, Ei, Jtotal

#print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:])/(10**(-9)))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
#print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')

def ejecutar_algoritmo(fixed_point, diametro, fixed_value, X, Y, R, Q, mov, m, delta, nodosy, nodosx,
                        dx, dy, windx, windy, max_iter_rho, max_iter, it_global, l, visualizacion, rho_h=False,
                          condct='si', copy='no', Muestra=False, fi=30):
    """Ejecuta el algoritmo completo ajustando Jp hasta cumplir la condición de convergencia."""
    convergencia_lograda = False
    cont = 0
    fi = 30
    pos_conductor1 =  fixed_point[0]
    pos_conductor2 =  fixed_point[1]
    print('----------------------------------')
    # Paso 1: Calcular potencial inicial
    Vmi = calcular_potencial_inicial(fixed_point, fixed_value, X, Y, R, Q, st='noV') # Ubicar en 5 nodos con el voltaje de +- en cada conductor
    Exx,  Eyy,  E = calcular_campo_electrico(Vmi, dx, dy)
    Ey = E[nod_central[1],nod_central[0]] # Componente vertical del campo eléctrico en punto medio entre los conductores, necesario  para cb carga en conductor
    # Paso 2: Inicializar densidades
    rho_inicial, rho_boundary, rho_ip, rho_in = inicializar_densidad(Ey,  diametro, R, m, delta, nodosy, nodosx, pos_conductor1,
                                                             pos_conductor2, con_condct=condct, copiado=copy)
    # Paso 3: Iterar para calcular Vm y rho
    Vm, rho_p, rho_n, rho_d = iterar_potencial(Vmi, rho_inicial, rho_boundary, [rho_ip, rho_in], fixed_point, dx, dy, windx, windy,
                                                max_iter_rho, max_iter, it_global, rho_h=rho_h, visu=visualizacion, con_condct=condct, copiado=copy, must = Muestra, fi=fi)
    #fi += 1
    visualizacion += 2
    # Paso 4: Calcular resultados finales
    Vol_def = Vm+Vmi
    Jave, Ji, Campo_ini, Campo_fin, Ei , Jtot= calcular_resultados_finales(Vm, Vmi, rho_p, rho_n, dx, dy, mov, windx, windy, l)
    '''
    print(f'Jp promedio calculado: {Jave * 1e9} nA/m^2')
    #print(f'Jp promedio propuesto: {Jp * 1e9} nA/m^2')
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
    
    '''
    return Campo_ini, Vmi, rho_p, rho_n, rho_d, Vm, Vol_def, Campo_fin, Ei, Ji, Jave, Jtot


#############################
########### PANEL DE CONTROL
############################
## Parámetros genericos
Rp = 10 * 10**(-2) # (m) radio equivalente del conductor positivo
Rn = 10 * 10**(-2) # (m) radio equivalente del conductor negativo
R = [Rp, Rn]
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol1 = 400000 # (V) voltaje de la línea
Vol2 = -400000 # V
ele = 1.6021*10**(-19) # (C) carga del electrón
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
movp = 3*10**(-4) # (m^2/Vs) movilidad iones positivos
movn = 5*10**(-4) # (m^2/Vs) movilidad iones negativos
mov = [movp, movn]
mp = 0.85 # (AD) factor de rugosidad conductor positivo
mn = 0.85 # (AD) factor de rugosidad conductor negativo
m = [mp, mn]
Rep = 1.8*10**(-6) # (m^3/2) Coeficiente de recombinación iones
#Ren = # (m^3/2) Coeficiente de recombinación iones
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
x_coor1 = -6
y_coor1 = 20
x_coor2 = 6
y_coor2 = 20
dia = np.abs(x_coor1-x_coor2)
l=1 # (m) altura  de medición de campos
fixed_value1 = Vol1
fixed_value2 = Vol2

# Discretización
Sx = 15 # (m) media longitud del plano de tierra 
Sy = 35 # (m) altura del área de estudio respecto de tierra
nodosy, nodosx = None, None # número de nodos en ambas dimensiones
# Definir coordenadas de la malla en x y y
x,y, nodosx, nodosy, Sx, Sy = discretiza(Rp, x_coor2, y_coor2, nodox=nodosx, nodoy=nodosy, sx=Sx, sy=Sy)  
X,Y =np.meshgrid(x,y)
dx = np.abs(x[1] - x[0]) # (m) distancia física entre nodos en cada columna
dy = np.abs(y[1] - y[0]) # (m) distancia física entre nodos en cada fila
posx_conductor1, posy_conductor1 = encuentra_nodos(x_coor1, y_coor1)
posx_conductor2, posy_conductor2 = encuentra_nodos(x_coor2, y_coor2)
fixed_point1 = (posy_conductor1, posx_conductor1)
fixed_point2 = (posy_conductor2, posx_conductor2)
fixed_point = [fixed_point1, fixed_point2]
fixed_value = [fixed_value1, fixed_value2]
# Obtener las coordenadas físicas de los nodos
pos_central = ((x_coor1+x_coor2)/2, (y_coor1+y_coor2)/2)
nod_central = encuentra_nodos(pos_central[0], pos_central[1])
# Calcular el punto central

# Ajuste de viento
wndx = 0 # m/s
wndy = 0 # m/s
modo = 'uniforme'
coef_drag = 0.1
wndx1 = np.ones((nodosy, nodosx)) * wndx
wndy1 = np.ones((nodosy, nodosx)) * wndy
windx, windy = windDist(wndx1, wndy1, l, Y, coef_drag, uni=modo) # Se calcula la distribución de viento

# Parámetros de iteración y convergencia
max_iter_rho = 600 # Número máximo de iteraciones para el cálculo de la densidad de carga
max_iter = 400 # Número máximo de iteraciones para el cálculo del potencial en el espacio
it_global = 4 # Número máximo de iteraciones para el cálculo del potencial definitivo
in_condct = 'si' # Condición de que se ubica la densidad de carga alrededor del conductor o solo en un nodo central
copiado = 'no' # Condición de que los valores sean copiados o no al momento de mantener las condiciones de borde en las mallas
histl = False # Si se desea crear la lista de evolución de la densidad de carga
visualizacion = 15 #  número de la ventana donde se muestra la evolución de la densidad de carga
mostrar = False # Permite mostrar los datos de la convergencia de rho inicial
fi = 30 # número de la ventana donde se muestra la evolución de la diferencia promedio y la desviación



##############################################
############# ALGORITMO PRINCIPAL Y EJECUCIÓN
###############################################

coordenada = [(x_coor1,y_coor1),(x_coor2,y_coor2)]
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
        for z in range(2):
            if j==i:
                P[i][j] = K*np.log(2*h[i]/R[z])
            else:
                if largo !=1 :
                    P[i][j] = K*np.log(D_im[i][j]/D[i][j])
#V = [Vol,Vol] # voltajes sobre los conductores reales
V = [Vol1,  Vol2]
Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores




Campo_ini, Vmi, rho_p, rho_n, rho_d, Vm, Vol_def, Campo_fin, Ei, Ji, Jave, Jtot = ejecutar_algoritmo(fixed_point, dia,
                                                                                                fixed_value, X, Y, R,
                                                                                                  Q, mov, m, delta, nodosy,
                                                                                                    nodosx, dx, dy, windx,
                                                                                                      windy, max_iter_rho,
                                                                                                        max_iter, it_global,
                                                                                                          l, visualizacion, rho_h=histl,
                                                                                                            condct=in_condct, copy=copiado,
                                                                                                              Muestra=mostrar, fi=fi)
Exxini, Eyyini, Em = Campo_ini
Edefx, Edefy, Edef = Campo_fin
Jave = np.mean(Jtot[encuentra_nodos(0,l)[1],:])

######## GRAFICOS


'''
class AsymmetricLogNorm(Normalize):
    def __init__(self, linthresh, vmin=None, vmax=None, base=10):
        super().__init__(vmin, vmax)
        self.linthresh = linthresh
        self.base = base

    def __call__(self, value, clip=None):
        # Normalización asimétrica para valores negativos y positivos
        result = np.zeros_like(value)
        mask_pos = value > self.linthresh
        mask_neg = value < -self.linthresh
        mask_lin = ~mask_pos & ~mask_neg
        
        # Escalas logarítmicas
        result[mask_pos] = np.log(value[mask_pos] / self.linthresh) / np.log(self.base)
        result[mask_neg] = -np.log(-value[mask_neg] / self.linthresh) / np.log(self.base)
        
        # Escala lineal cerca de 0
        result[mask_lin] = value[mask_lin] / self.linthresh
        return result
'''



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
    # Define el límite para la transición a la escala lineal
    linthresh = 1e+2  # Ajusta este valor según tus datos
    # Usa SymLogNorm para manejar valores negativos y positivos
    norm = SymLogNorm(linthresh=linthresh, vmin=Vmi.min(), vmax=Vmi.max(), base=50)
    # Genera el gráfico
    mesh = plt.pcolormesh(X, Y, Vmi, cmap='plasma', shading='auto', norm=norm)
    # Añade la barra de colores
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Potencial $kV$')
    # Define manualmente los ticks del colorbar, incluyendo el mínimo, el máximo y valores intermedios
    min_val = Vmi.min()
    max_val = Vmi.max()
    ticks = [min_val, -linthresh, 0, linthresh, max_val]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.2f}' for tick in ticks])  # Divide por 1000 si es necesario
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Potencial electrostático', fontsize=15)
    plt.tight_layout()


def grafRhop(num):
    plt.figure(num)
    #linthresh = 1e-8  # Ajusta según la magnitud de rho_p
    min_linthresh = 1e-5  # Un valor mínimo positivo
    linthresh = max(np.abs(min(np.abs(np.min(rho_p)), np.max(rho_p)) / 5), min_linthresh)
    #linthresh = abs(np.min(rho_p)) * 1
    #norm = AsymmetricLogNorm(linthresh=linthresh, vmin=rho_p.min(), vmax=rho_p.max())
    norm = SymLogNorm(linthresh=linthresh, vmin=rho_p.min(), vmax=rho_p.max(), base=10)
    mesh = plt.pcolormesh(X, Y, rho_p, cmap='viridis', shading='auto', norm=norm)
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Densidad de carga positiva $C/m^3$', fontsize=11)
    ticks = [rho_p.min(), -linthresh, 0, linthresh, rho_p.max()]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick:.2e}' for tick in ticks])
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Densidad de carga positiva', fontsize=15)
    plt.tight_layout()

def grafRhon(num):
    plt.figure(num)
    #linthresh = 1e-8  # Ajusta según la magnitud de rho_n
    min_linthresh = 1e-5  # Un valor mínimo positivo
    linthresh = max(np.abs(min(np.abs(np.min(rho_n)), np.max(rho_n)) / 5), min_linthresh)
    norm = SymLogNorm(linthresh=linthresh, vmin=rho_n.min(), vmax=rho_n.max(), base=10)
    mesh = plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto', norm=norm)
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Densidad de carga negativa $C/m^3$', fontsize=11)
    ticks = [rho_n.min(), -linthresh, 0, linthresh, rho_n.max()]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick:.2e}' for tick in ticks])
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Densidad de carga negativa', fontsize=15)
    plt.tight_layout()

def grafRhoT(num):
    plt.figure(num)
    min_linthresh = 1e-9  # Un valor mínimo positivo
    linthresh = max(np.abs(min(np.abs(np.min(rho_d)), np.max(rho_d)) / 5), min_linthresh)
    norm = SymLogNorm(linthresh=linthresh, vmin=rho_d.min(), vmax=rho_d.max(), base=10)
    mesh = plt.pcolormesh(X, Y, rho_d, cmap='viridis', shading='auto', norm=norm)
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Densidad de carga total $C/m^3$', fontsize=11)
    ticks = [rho_d.min(), -linthresh, 0, linthresh, rho_d.max()]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick:.2e}' for tick in ticks])
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Densidad de carga final', fontsize=15)
    plt.tight_layout()

def grafVf(num):
    plt.figure(num)
    linthresh = 1e+2  # Ajusta según la magnitud de Vm
    norm = SymLogNorm(linthresh=linthresh, vmin=Vm.min(), vmax=Vm.max(), base=10)
    mesh = plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto', norm=norm)
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Potencial iónico $kV$')
    ticks = [Vm.min(), -linthresh, 0, linthresh, Vm.max()]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.2f}' for tick in ticks])
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Potencial iónico', fontsize=15)
    plt.tight_layout()

def grafVdef(num):
    plt.figure(num)
    linthresh = 1e+2  # Ajusta según la magnitud de Vol_def
    norm = SymLogNorm(linthresh=linthresh, vmin=Vol_def.min(), vmax=Vol_def.max(), base=10)
    mesh = plt.pcolormesh(X, Y, Vol_def, cmap='plasma', shading='auto', norm=norm)
    cbar = plt.colorbar(mesh)
    cbar.set_label(r'Potencial definitivo $kV$')
    ticks = [Vol_def.min(), -linthresh, 0, linthresh, Vol_def.max()]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.2f}' for tick in ticks])
    plt.xlabel('Distancia horizontal (m)', fontsize=11)
    plt.ylabel('Distancia vertical (m)', fontsize=11)
    plt.title('Potencial definitivo', fontsize=15)
    plt.tight_layout()

def grafEf(nm):
    plt.figure(nm)
    #plt.quiver(Xe, Ye, Edefx[1:-1, 1:-1], Edefy[1:-1, 1:-1], Edef[1:-1, 1:-1], cmap='plasma')
    plt.quiver(X, Y, Edefx, Edefy, Edef, cmap='plasma', scale_units='xy')
    # Agregar la barra de colores
    cbar = plt.colorbar()
    cbar.set_label(r'Magnitud campo eléctrico $kV/m$', fontsize=11)
    # Obtén los ticks actuales
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    # Mostrar el gráfico
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Distancia vertical (m)',fontsize=11)
    plt.title(r'Campo eléctrico definitivo (V/m)', fontsize=13)
    plt.tight_layout()

def grafJ1(num):
    plt.figure(num)
    #plt.plot(x[30:-30], Ji[30:-30]*(10**9))
    plt.plot(x[30:-30], Jtot[encuentra_nodos(0, l)[1],30:-30]*(10**9))
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
    plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
    plt.grid(True)

def grafE1(num):
    plt.figure(num)
    plt.plot(x[30:-30], -Edefy[encuentra_nodos(0, l)[1],30:-30]/1000)
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
    plt.title(r'Componente $E_y$ campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$E_y$ = {str(np.round(np.mean(Edefy[encuentra_nodos(0, l)[1],:]/1000),5))} kV'])
    plt.grid(True)

gra = ['Eele', 'Vele', 'Rhop', 'Rhon', 'Rhof', 'Vf', 'Vdef', 'Ef', 'E1', 'J1']

def show_plot(graf):
    for i in graf:
        if i == 'Eele':
            grafE(1)
        elif i == 'Vele':
            grafV(2)
        elif i == 'Rhop':
            grafRhop(3)
        elif i == 'Rhon':
            grafRhon(4)
        elif i == 'Rhof':
            grafRhoT(5)
        elif i == 'Vf':
            grafVf(6)
        elif i == 'Vdef':
            grafVdef(7)
        elif i == 'Ef':
            grafEf(8)
        elif i == 'J1':
            grafJ1(9)
            # Crear la figura
            #plt.figure(30)
            fig = plt.figure(figsize=(12, 12))  # Tamaño ajustado para los subgráficos

            # Subgráfico 1
            ax1 = fig.add_subplot(131, projection='3d')  # 1 fila, 2 columnas, gráfico 1
            surf1 = ax1.plot_surface(X, Y, rho_d*10**(6), cmap='viridis', edgecolor='none')
            fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)  # Barra de color
            ax1.set_title(r'Densidad de carga iónica total')
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_zlabel(r'$\rho (\mu C/m^3)$')
            # Subgráfico 2
            ax2 = fig.add_subplot(132, projection='3d')  # 1 fila, 2 columnas, gráfico 1
            surf2 = ax2.plot_surface(X, Y, rho_n*10**(6), cmap='viridis', edgecolor='none')
            fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)  # Barra de color
            ax2.set_title(r'Densidad de carga iónica negativa')
            ax2.set_xlabel('X')
            ax2.set_ylabel('Y')
            ax2.set_zlabel(r'$\rho (\mu C/m^3)$')
            # Subgráfico 3
            ax3 = fig.add_subplot(133, projection='3d')  # 1 fila, 2 columnas, gráfico 1
            surf3 = ax3.plot_surface(X, Y, rho_p*10**(6), cmap='viridis', edgecolor='none')
            fig.colorbar(surf3, ax=ax3, shrink=0.5, aspect=10)  # Barra de color
            ax3.set_title(r'Densidad de carga iónica positiva')
            ax3.set_xlabel('X')
            ax3.set_ylabel('Y')
            ax3.set_zlabel(r'$\rho (\mu C/m^3)$')
        elif i == 'E1':
            grafE1(11)
        # Crear la figura
            #plt.figure(21)
            fig = plt.figure(figsize=(12, 12))  # Tamaño ajustado para los subgráficos

            # Subgráfico 1
            ax1 = fig.add_subplot(131, projection='3d')  # 1 fila, 2 columnas, gráfico 1
            surf1 = ax1.plot_surface(X, Y, Vmi/1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)  # Barra de color

            ax1.set_title(r'Potencial electrostático')
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_zlabel(r'V(kV)')

            # Subgráfico 2
            ax2 = fig.add_subplot(132, projection='3d')  # 1 fila, 3 columnas, gráfico 2
            surf2 = ax2.plot_surface(X, Y, Vm/1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)  # Barra de color

            ax2.set_title('Potencial iónico')
            ax2.set_xlabel('X')
            ax2.set_ylabel('Y')
            ax2.set_zlabel('V(kV)')

            # Subgráfico 3
            ax3 = fig.add_subplot(133, projection='3d')  # 1 fila, 3 columnas, gráfico 2
            surf3 = ax3.plot_surface(X, Y, Vol_def/1000, cmap='plasma', edgecolor='none')
            fig.colorbar(surf3, ax=ax3, shrink=0.5, aspect=10)  # Barra de color

            ax3.set_title('Potencial definitivo')
            ax3.set_xlabel('X')
            ax3.set_ylabel('Y')
            ax3.set_zlabel('V(kV)')
            
            plt.tight_layout()
        plt.show()

show_plot(gra)
