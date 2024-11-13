import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')  # O usa 'Qt5Agg' o 'Agg' si no funciona
import matplotlib.pyplot as plt
import math as ma
#import threading
from matplotlib.colors import LogNorm 
plt.close('all')
################################
## Parámetros genericos
R = 6.35 * 10**(-2) # (m) radio equivalente del conductor
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol = 300000 # (V) voltaje de la línea
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
mov = 7*10**(-4) # (m^2/Vs)
m = 0.7 # (AD) factor de rugosidad
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
wndx = 10 # m/s
wndy = 0 # m/s
# Sx=3H
# Sy = 2H
Sx = 60 # (m) media longitud del plano de tierra 
Sy = 70 # (m) altura del área de estudio respecto de tierra
l = 1 # (m) distancia desde el suelo a la altura de interés de influencia de campo
## Definición coordenadas conductores caso bipolar
#coordenada = [(6,7), (-6,7)] # (m) coordenadas cargas en posición de los conductores
x_coor = 0
y_coor = 30
coordenada = [(x_coor,y_coor)]
coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
#Hor = np.linspace(-10,10,100) # (m) Ancho total del área de estudio
#Ver = 1 # (m) Altura donde interesa medir el campo eléctrico
Tol = 10**(-4)
con = 0
max_iter_rho = 300
max_iter = 800
nodosx = 100
nodosy = 100
def mod(z1,z2):
    return np.sqrt((z1[0]-z2[0])**2 + (z1[1]-z2[1])**2)

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


windx = np.ones((nodosx, nodosy)) * wndx
windy = np.ones((nodosx, nodosy)) * wndy
# Definir coordenadas de la malla en x y y
x = np.linspace(-Sx, Sx, nodosy)
y = np.linspace(Sy, 0, nodosx)
dx = np.abs(x[1] - x[0])
dy = np.abs(y[1] - y[0])
def calcular_campo_electrico_inicial(nodosx, nodosy, x, y, w, h, Q, K):
    # Inicializar matrices para las componentes y magnitud del campo eléctrico
    E = np.zeros((nodosx, nodosy))
    Exx = np.zeros((nodosx, nodosy))
    Eyy = np.zeros((nodosx, nodosy))
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

Exx, Eyy, E = calcular_campo_electrico_inicial(nodosx,nodosy, x, y, w, h, Q, K)

X,Y =np.meshgrid(x,y)
U = Exx
Ww = Eyy
print('Campo eléctrico  electrostático libre iones calculado')




posx_conductor = int((x_coor + Sx) / dx)  # Para x_conductor = 0
posy_conductor = int((Sy - y_coor) / dy)
fixed_point = (posy_conductor, posx_conductor)
fixed_value = Vol

## Resolución para densidad de carga en el espacio


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

def rhoA(Sx, Jp, b, a,  m_ind, d):
    # Cálculo condición de borde densidad de carga conductor energizado
    # Depende del valor de Jp
    # Asume que existe efecto corona
    rho_i = (Sx*Jp)*10**(-3)/(np.pi*b*(100*a)*m_ind*(30*d + 9*np.sqrt(d/(100*a)))) # radio está en cm
    return rho_i

y_p = np.zeros_like(Y)
indices_sup = (np.abs(Y[:,0] - int(Sy-y_coor))).argmin()
indices_inf = (np.abs(Y[:,0] - int(y_coor))).argmin()
y_p[:(nodosy-indices_sup),:] = Y[indices_sup:,:]
y_p[(nodosy-indices_sup):,:] = Y[(indices_inf-len(y_p[(nodosy-indices_sup):,0])):indices_inf,:]
#y_f = np.zeros_like(Y)
#y_f[indices_inf:,:] = Y[(indices_inf-(nodosy-2)):indices_inf,:]
Ya = (np.flip(Y)-(Sy-y_coor))

def potencial_electrostático(f_point, f_value, X, Y, radio, ind, carga=None):
    Vm = np.zeros((nodosx,nodosy))
    Vm2 = Vm.copy()
    ## Condiciones de borde potencial inicial
    Vm[f_point] = f_value
    Vm[-1,:] = 0 # borde inferior
    Vm2[f_point] = f_value
    Vm2[-1,:] = 0 # borde inferior
    if carga is not None:
        # Calcular la distancia entre cada punto de la malla y el punto (w, h)
        Mod_coor = np.sqrt((X - w)**2 + (Y - h)**2)
        Mod_coor2 = np.sqrt((X - w)**2 + (Y + h)**2)
        # Calcular el potencial eléctrico en cada nodo de la malla
        #Vm = carga * 0.5 * K*(1/Mod_coor - 1/Mod_coor2)
        Vm = carga*K*(np.log(1/Mod_coor) - np.log(1/Mod_coor2))
        #Vm = calcular_potencial_con_borde_vectorizado(Exx, Eyy,dx, dy)
    else:
        c1 = np.cosh(np.pi*(X-2*ind*Sx)/(2*(Sy)))
        c2 = np.cos(np.pi*Ya/(2*(Sy)))
        c22 = np.cos(np.pi*y_p/(2*(Sy)))
        c3 = np.cosh(np.pi*ind*Sx/(Sy))
        c4 = np.cos(np.pi*radio/(2*(Sy)))
        Vm = np.abs(f_value * np.log((c1-c2)/(c1+c2))/np.log((c3-c4)/(c3+c4)))
        Vm2 = np.abs(f_value * np.log((c1-c22)/(c1+c22))/np.log((c3-c4)/(c3+c4)))
    Vm[f_point] = f_value
    Vm[-1,:] = 0 # borde inferior
    Vm2[f_point] = f_value
    Vm2[-1,:] = 0 # borde inferior
    return Vm,Vm2

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

def dev_1sentido(Ex, Ey, ep0, rhoi_1, rhoj_1, dx, dy):
    alpha = Ey*ep0/(2*dy) + Ex*ep0/(2*dx)
    beta = ep0*(Ey*rhoi_1/dy + Ex*rhoj_1/dx)
    dis = np.abs(alpha**2 + beta)
    G = -alpha + np.sqrt(dis)
    #G = np.sqrt(dis)
    #G[np.isnan(G)] = 0
    return G
def dev_central(Ex, Ey,ep0, rhodx, rhoix, rhoiy, rhosy, dx, dy):
    d_rho_dx = (rhodx - rhoix) / (2 * dx)
    d_rho_dy = (rhoiy - rhosy) / (2 * dy)
    G = np.sqrt(np.abs(-ep0 * (Ex * d_rho_dx + Ey * d_rho_dy)))
    return G
def dev_triple_di(Ex, Ey, ep0, rhoys, rhoyi, rhox, dx, dy):
    alpha  = Ex*ep0/(2*dx)
    beta = -ep0*(Ey*rhoyi/(2*dy) - (Ex*rhox/dx + Ey*rhoys/(2*dy)))
    G = -alpha + np.sqrt(np.abs(alpha**2 + beta))
    return G
def dev_triple_si(Ex, Ey, ep0, rhoxd, rhoxi, rhoy, dx, dy):
    alpha  = Ey*ep0/(2*dy)
    beta = -ep0*(Ex*rhoxi/(2*dx) - (Ey*rhoy/dy + Ex*rhoxd/(2*dx)))
    G = -alpha + np.sqrt(np.abs(alpha**2 + beta))
    return G

# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
def update_rho_vectorized(rhoini, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_bound, met=2):
    #rho_new = np.copy(rho)
    rho = rhoini.copy()
    rho[fixed_point] = rho_bound[fixed_point]  # Mantener la condición en el punto central
    rho[-1, 1:-1] = rho_bound[-1,1:-1]  # Mantener la condición en el borde inferior
    Ewx = Ex + wx/mov
    Ewy = Ey + wy/mov
    cx =fixed_point[1] # índice posición en x
    cy =fixed_point[0] # índice posición en y
    rhoN = rho.copy()
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
        rhoN[1:-1,0] = dev_triple_di(Ewx[1:-1,0], Ewy[1:-1,0], epsilon0, rho[:-2, 0], rho[2:, 0], rho[1:-1, 1], dx, dy) # borde izquierdo
        #Gd = -epsilon0 * (Ewx[1:-1, -1] * (rho[1:-1, -1] - rho[1:-1, -2]) / dx + Ewy[1:-1, -1] * (rho[2:, -1] - rho[:-2, -1]) / (2 * dy))
        rhoN[1:-1,-1] = dev_triple_di(Ewx[1:-1,-1], Ewy[1:-1,-1], epsilon0, rho[:-2, -1], rho[2:, -1], rho[1:-1, -2], dx, dy) # borde derecho
        #Gs = -epsilon0 * (Ewx[0, 1:-1] * (rho[0, 2:] - rho[0, :-2]) / (2 * dx) + Ewy[0, 1:-1] * (rho[1, 1:-1] - rho[0, 1:-1]) / dy)
        rhoN[0, 1:-1] = dev_triple_si(Ewx[0, 1:-1], Ewy[0, 1:-1], epsilon0, rho[0, 2:], rho[0, :-2], rho[1, 1:-1], dx, dy) # borde superior
        #Geis = -epsilon0 * (Ewx[0, 0] * (rho[0, 1] - rho[0, 0]) / dx + Ewy[0, 0] * (rho[1, 0] - rho[0, 0]) / dy)
        rhoN[0, 0] = dev_1sentido(Ewx[0, 0], Ewy[0, 0], epsilon0, rho[1, 0], rho[0, 1], dx, dy) # esquina superior izquierda
        #Geds = -epsilon0 * (Ewx[0, -1] * (rho[0, -1] - rho[0, -2]) / dx + Ewy[0, -1] * (rho[1, -1] - rho[0, -1]) / dy)
        rhoN[0, -1] = dev_1sentido(Ewx[0, -1], Ewy[0, -1], epsilon0, rho[1, -1], rho[0, -2], dx, dy) # esquina superior derecha
        #Geii = -epsilon0 * (Ewx[-1, 0] * (rho[-1, 1] - rho[-1, 0]) / dx + Ewy[-1, 0] * (rho[-1, 0] - rho[-2, 0]) / dy)
        rhoN[-1, 0] = dev_1sentido(Ewx[-1, 0], Ewy[-1, 0], epsilon0, rho[-2, 0], rho[-1, 1], dx, dy) # esquina inferior izquierda
        #Gedi = -epsilon0 * (Ewx[-1, -1] * (rho[-1, -1] - rho[-1, -2]) / dx + Ewy[-1, -1] * (rho[-1, -1] - rho[-2, -1]) / dy)
        rhoN[-1, -1] = dev_1sentido(Ewx[-1, -1], Ewy[-1, -1], epsilon0, rho[-2, -1], rho[-1, -2], dx, dy) # esquina inferior derecha
    # Aplicar condiciones de borde
    rhoN[fixed_point] = rho_bound[fixed_point]  # Mantener la condición en el punto central
    rhoN[-1, 1:-1] = rho_bound[-1,1:-1]  # Mantener la condición en el borde inferior
    
    return rhoN

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
V_b = np.zeros_like(Vm)
V_b[-1,:] = 0
# Aplicar condiciones de borde tipo Neumann
# Borde izquierdo (columna 0)
V_b[:, 0] = V_b[:, 1] + dx * Exx[:, 0]
# Borde derecho (última columna)
V_b[:, -1] = V_b[:, -2] - dx * Exx[:, -1]
# Borde superior (fila 0)
V_b[0, :] = V_b[1, :] + dy * Eyy[0, :]
'''


# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
def algoritmo_rho_v(V, rho_ini, dx, dy, windx, windy, max_iter_rho, Jplate, rho_A,met):
    conv = 0
    rho_b =np.zeros((nodosx,nodosy))
    # En un comienzo se tiene el campo eléctrico electrostático
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    # Se define a la malla rho_b con las condiciones de borde
    rho_b[fixed_point] = rho_A
    # COndición de desplazamiento de la densidad de carga debido al viento
    rho_b[-1,1:-1] = Jplate/(mov*np.sqrt(np.abs((Eyyi[-1,1:-1]+windy[-1,1:-1]/mov)**2+(Exxini[-1,1:-1]+windx[-1,1:-1]/mov)**2))) # se postula que la dirección de E será solamente vertical
    rho1 = rho_ini.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
    '''
    plt.figure(12)
    plt.pcolormesh(X, Y, rho1, cmap='viridis', shading='auto',norm=LogNorm())
    plt.show()
    '''
    for iteration in range(max_iter_rho):
        rho_0 = rho1.copy()
        rho1 = update_rho_vectorized(rho1, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_b, met = met)
        # Criterio de convergencia
        # Define un pequeño valor para evitar divisiones por cero
        condicion,diff = convergencia(rho1, rho_0, 0.01)
        '''
        epsilon = 1e-10
        epsilon_array = np.full_like(rho1, epsilon)
        # Calcula la diferencia relativa usando el máximo entre |rho1|, |rho_0| y epsilon
        diferencia_relativa = np.abs(rho1 - rho_0) / np.maximum.reduce([np.abs(rho1), np.abs(rho_0), epsilon_array])
        # Verifica si todas las diferencias relativas son menores al 10%
        condicion = np.all(diferencia_relativa < 0.01)
        diff = np.linalg.norm(diferencia_relativa,ord=np.inf)
        '''
        #print(r'Tolerancia $\rho_1$ = '+str(diff))
        if condicion:
            print(f"Convergencia para rho alcanzada en la iteración {iteration}")
            break
    #print(r'Tolerancia rho_1 = '+str(diff))
    #print(r'última iteración = '+str(iteration))
    return rho1
    

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
        (u[0:-2, 1:-1] - 2 * u[1:-1, 1:-1] + u[2:, 1:-1]) / dy2 +  # Término en x
        (u[1:-1, 0:-2] - 2 * u[1:-1, 1:-1] + u[1:-1, 2:]) / dx2    # Término en y
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
    u[fixed_point] = fixed_value

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

def densidad_voltaje(Volt, rho):
    # Calcula la densidad de cargas usando la ecuación de Poisson según el potencial ya calculado
    Rho_new = np.zeros_like(rho)
    Rho_new = -epsilon0*apply_laplacian_asymmetric(Volt, dx, dy)
    Rho_new[fixed_point] = rho_i
    Rho_new[0,:] = rho[0,:]
    Rho_new[-1,:] = rho[-1,:]
    Rho_new[:,0] = rho[:,0]
    Rho_new[:,-1] = rho[:,-1]
    return Rho_new
# Resolver usando FMG
# Algoritmo que obtien la red de potencial eléctrico en base a los valores de Rho
def algoritmo_V_rho(V, rho1, dx, dy, fixed_point, fixed_value, max_iter):
    f_rhs = funcion_f(rho1)
    V_b =  Vm.copy()
    for iteration in range(max_iter):
        Vold = V.copy()
        V = update_v(V, f_rhs, dx,dy,fixed_point=fixed_point, fixed_value=fixed_value, V_boundary=V_b, Exx=Exx,Eyy=Eyy)
        condicion,diff = convergencia(V, Vold, 0.01)
        if condicion:
            #print(f"Convergencia alcanzada para V en la iteración {iteration}")
            break
    return V

# ALGORITMO PRINCIPAL
# Se calcula primero el potencial incial en base a CSM
## Cálculo de potencial eléctrico inicial en base a la carga calculada por CSM
Vmi,Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, R, 0, carga=Q)
##### Resolución ecuación de continuidad
Jp = 0.6744*10**(-8) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
# Condiciones de borde
rho_i = rhoA(Sx, Jp, mov, R, m, delta)
rho_inicial = np.zeros((nodosx, nodosy))
Exxini, Eyyini, Em = calcular_campo_electrico(Vmi, dx, dy) # campo electrostático inicial
rho_inicial[fixed_point] = rho_i
rho_inicial[-1,1:-1] = Jp/(mov*np.sqrt(np.abs((Eyyini[-1,1:-1]+windy[-1,1:-1]/mov)**2+(Exxini[-1,1:-1]+windx[-1,1:-1]/mov)**2)))
#print(Vmi)
Vm = Vmi.copy()
it_global = 260
for n in range(it_global):
    Volder = Vm.copy()
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    if n==0:
        rho_n = algoritmo_rho_v(Vmi, rho_inicial, dx, dy, windx, windy, max_iter_rho, Jp, rho_i,met=0)
    else:
        rho_n = algoritmo_rho_v(Vm, rho_n, dx, dy, windx, windy, max_iter_rho, Jp, rho_i, met=0)
    Vm = algoritmo_V_rho(Vm, rho_n, dx, dy, fixed_point, fixed_value, max_iter)
    condicion,diff = convergencia(Vm, Volder, 0.01)
    if n%20 == 0:
        print(r'Diferencia relativa V y Vold: '+str(diff))
    if condicion:
        print(f"Convergencia alcanzada para V en la iteración {n}")
        break
#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
print('Potencial calculado')
##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica
#Voldef = Vm + Vmi
Edefx, Edefy, Edef = calcular_campo_electrico(Vm, dx, dy)
J = rho_n*mov*np.sqrt((Edefx+(windx/mov))**2 + (Edefy+(windy/mov))**2)
Ei = Edef[int((Sy-l)/dy),:] # Magnitud Campo eléctrico a nivel de piso
Ji = J[int((Sy-l)/dy),:] # Densidad de corriente a nivel de piso
Jave = np.mean(Ji)
print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:])/(10**(-9)))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')
#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')


######## GRAFICOS
# Lista gráficos
def grafE():
    plt.figure(1)
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

def grafV():
    plt.figure(2)
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

def grafRho():
    plt.figure(3)
    plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto',norm=LogNorm())
    #plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto')
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    plt.title('Densidad de carga final',fontsize=15)
    # Añadir una barra de colores para mostrar la escala
    cbar = plt.colorbar()
    cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)

def grafVf():
    plt.figure(4)
    #plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
    plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto',norm=LogNorm())
    plt.title('Potencial',fontsize=15) 
    plt.xlabel('Distancia horizontal (m)',fontsize=11)
    plt.ylabel('Distancia vertical (m)',fontsize=11)
    # Añadir una barra de colores para mostrar la escala
    cbar = plt.colorbar()
    cbar.set_label(r'Potencial $kV$')
    ticks = cbar.get_ticks()
    # Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
    plt.tight_layout
    #plt.show()

def grafEf():
    plt.figure(5)
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
    plt.title(r'Magnitud de campo eléctrico a nivel de suelo (V/m)', fontsize=13)
    plt.tight_layout()

def grafJ1():
    plt.figure(6)
    plt.plot(x[10:-10], Ji[10:-10]*(10**9))
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
    plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
    plt.grid(True)

def grafE1():
    plt.figure(7)
    plt.plot(x[10:-10], Ei[10:-10]/1000)
    plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
    plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
    plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
    plt.tight_layout()
    plt.legend([f'$|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV'])
    plt.grid(True)



gra = ['Eele', 'Vele', 'Rhof', 'Vf', 'Ef', 'E1', 'J1']

def show_plot(graf):
    for i in graf:
        if i == 'Eele':
            grafE()
        elif i == 'Vele':
            grafV()
        elif i == 'Rhof':
            grafRho()
        elif i == 'Vf':
            grafVf()
        elif i == 'Ef':
            grafEf()
        elif i == 'J1':
            grafJ1()
        elif i == 'E1':
            grafE1()
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
#threading.Thread(target=show_plot).start()