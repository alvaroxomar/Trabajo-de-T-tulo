import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.colorbar as cbar
plt.close('all')
################################
## Parámetros genericos
Rp = 6.35 * 10**(-2) # (m) radio equivalente del conductor positivo
Rn = 6.35 * 10**(-2) # (m) radio equivalente del conductor negativo
R = [Rp, Rn]
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol1 = 300000 # (V) voltaje de la línea
Vol2 = -300000 # V
ele = -1.6021*10**(-19) # (C) carga del electrón
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
movp = 1.5*10**(-4) # (m^2/Vs) movilidad iones positivos
movn = 2*10**(-4) # (m^2/Vs) movilidad iones negativos
mp = 0.85 # (AD) factor de rugosidad conductor positivo
mn = 0.85 # (AD) factor de rugosidad conductor negativo
Rep = 1.8*10**(-6) # (m^3/2) Coeficiente de recombinación iones
#Ren = # (m^3/2) Coeficiente de recombinación iones
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
wndx = 0 # m/s
wndy = 0 # m/s

def windDist(wndx, wndy, hr, y, alpha, uni):
    # Input: wndx: valocidad media en x, wndy: velocidad media en y, hr: altura de referencia
    # y: Matriz de alturas, alpha: coeficiente de rugosidad del piso
    if uni==1:
        Wx = wndx*(y/hr)**alpha
        Wy = wndy
    elif uni==0:
        Wx = wndx
        Wy = wndy
    return Wx, Wy
# Sx=3H
# Sy = 2H
Sx = 60 # (m) media longitud del plano de tierra 
Sy = 70 # (m) altura del área de estudio respecto de tierra
l = 0.3 # (m) distancia desde el suelo a la altura de interés de influencia de campo
## Definición coordenadas conductores caso bipolar
#coordenada = [(6,7), (-6,7)] # (m) coordenadas cargas en posición de los conductores
x_coor1 = -20
y_coor1 = 30
x_coor2 = 20
y_coor2 = 30
coordenada = [(x_coor1,y_coor1),(x_coor2,y_coor2)]
coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
#Hor = np.linspace(-10,10,100) # (m) Ancho total del área de estudio
#Ver = 1 # (m) Altura donde interesa medir el campo eléctrico
Tol = 10**(-4)
con = 0
max_iter_rho = 800
max_iter = 1000
nodosx = 100
nodosy = 100
windx = np.ones((nodosx, nodosy)) * wndx
windy = np.ones((nodosx, nodosy)) * wndy
# Definir coordenadas de la malla en x y y
x = np.linspace(-Sx, Sx, nodosy)
y = np.linspace(Sy, 0, nodosx)
dx = np.abs(x[1] - x[0])
dy = np.abs(y[1] - y[0])

# Encontrar los índices más cercanos
def encuentra_nodos(x0,y0):
    indice_x0 = (np.abs(x - x0)).argmin()
    indice_y0 = (np.abs(y - y0)).argmin()
    return (indice_x0,indice_y0)

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
        for z in range(2):
            if j==i:
                P[i][j] = K*np.log(2*h[i]/R[z])
            else:
                if largo !=1 :
                    P[i][j] = K*np.log(D_im[i][j]/D[i][j])
#V = [Vol,Vol] # voltajes sobre los conductores reales
V = [Vol1,  Vol2]
Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores



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





#posx_conductor = int(Sx/dx)
#posy_conductor = int((Sy - coordenada[0][1])/dy)
posx_conductor1 = int((x_coor1 + Sx) / dx)  # Para x_conductor = 0
posy_conductor1 = int((Sy - y_coor1) / dy)
posx_conductor2 = int((x_coor2 + Sx) / dx)  # Para x_conductor = 0
posy_conductor2 = int((Sy - y_coor2) / dy)
fixed_point1 = (posy_conductor1, posx_conductor1)
fixed_point2 = (posy_conductor2, posx_conductor2)
fixed_value1 = Vol1
fixed_value2 = Vol2
fixed_point = [fixed_point1, fixed_point2]
fixed_value = [fixed_value1, fixed_value2]

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

#def rhoA(Sx, Jp, b, a,  m_ind, d,I2l):
#rho_ip = (I2l)*10**(-3)/(np.pi*b[0]*(100*a[0])*m_ind[0]*(30*d + 9*np.sqrt(d/(100*a[0])))) # radio está en cm
    #rho_in = -(Sx*Jp-I2l/2)*10**(-3)/(np.pi*b[1]*(100*a[1])*m_ind[1]*(30*d + 9*np.sqrt(d/(100*a[1])))) # radio está en cm
def rhoA(Ey,D, ep0, V0_p,V0_s, V, Ecrit_p,Ecrit_s, R):    
    '''
    Ey: Campo  de carga libre en medio camino entre los dos conductores
    Ecrit_p,_s:  Campo gradiente corona para conductor positivo y negativo
    V0_p,_s: Voltajes de gradientecorona  en los conductores positivosy negativos, estos incluyen al factor de rugosidad
    '''
    rho_ip = Ey*D*8*ep0*V0_p*(V[0]-V0_p)/(Ecrit_p*R[0]*D**2*V[0]*(5-4*V0_p/V[0]))
    rho_in = Ey*D*8*ep0*V0_s*(V[1]-V0_s)/(Ecrit_s*R[1]*D**2*V[1]*(5-4*V0_p/V[1]))
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
def potencial_electrostático(f_point, f_value, X, Y, radio, ind, carga=None):
    Vm = np.zeros((nodosx,nodosy))
    Vm2 = Vm.copy()
    ## Condiciones de borde potencial inicial
    for g in range(len(f_point)):
        Vm[f_point[g]] = f_value[g]
        Vm[-1,:] = 0 # borde inferior
        Vm2[f_point[g]] = f_value[g]
        Vm2[-1,:] = 0 # borde inferior
    if carga is not None:
        print(carga)
        for z in range(len(carga)):
            # Calcular la distancia entre cada punto de la malla y el punto (w, h)
            Mod_coor = np.sqrt((X - w[z])**2 + (Y - h[z])**2)
            Mod_coor2 = np.sqrt((X - w[z])**2 + (Y + h[z])**2)
            print(Mod_coor)
            # Calcular el potencial eléctrico en cada nodo de la malla
            #Vm += carga[z] * 0.5 * K*(1/Mod_coor - 1/Mod_coor2)
            Vm += carga[z]*K*(np.log(1/Mod_coor) - np.log(1/Mod_coor2))
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
        Vm[f_point[g]] = f_value[g]
        Vm[-1,:] = 0 # borde inferior
        Vm2[f_point[g]] = f_value[g]
        Vm2[-1,:] = 0 # borde inferior
    return Vm,Vm2

def resuelve2(a, b, c, d, g, f):
    A2 = d*(d-g*b/a)
    B2 = f*(2*d-g*b/a)+2*c*g**2
    C2 = f**2
    Gr = (-B2 + np.sqrt(B2**2-4*A2*C2+0j))/(2*A2)
    ysol = np.abs(np.sqrt(Gr))
    Gy = (b*ysol)**2 - 4*a*c
    xsol = np.abs((-b*ysol + np.sqrt(Gy+0j))/(2*a))
    return xsol,ysol

# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
def update_rho_vectorized(rhoinip, rhoinin, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_boundp, rho_boundn):
    # Aquí hace falta copiar los rhop y rhon por lo que debo integrarlos en los parámetros de estas función.
    #rho_new = np.copy(rho)
    rhop = rhoinip.copy()
    rhon = rhoinin.copy()
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
    rp, rn = resuelve2(a,b,c,d,g,f)
    rhop[1:-1,1:-1] = rp
    rhon[1:-1,1:-1] = rn
    '''
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
    '''
    # Aplicar condiciones de borde
    rhop[fixed_point[0]] = rho_boundp[fixed_point[0]]  # Mantener la condición en el o los puntos centrales
    rhop[-1, :] = rho_boundp[-1,:]  # Mantener la condición en el borde inferior
    rhon[fixed_point[1]] = rho_boundn[fixed_point[1]]  # Mantener la condición en el o los puntos centrales
    rhon[-1, :] = rho_boundn[-1,:]  # Mantener la condición en el borde inferior
    return rhop, rhon

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
# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
def algoritmo_rho_v(V, rho_ini, dx, dy, windx, windy,max_iter_rho, Jplatep, Jplaten, rho_A):
    conv = 0
    rho_bp =np.zeros((nodosx,nodosy))
    rho_bn =np.zeros((nodosx,nodosy))
    # En un comienzo se tiene el campo eléctrico electrostático
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    # Se define a la malla rho_b con las condiciones de borde
    rho_bp[fixed_point[0]] = rho_A[0]
    rho_bn[fixed_point[1]] = rho_A[1]
    # COndición de desplazamiento de la densidad de carga debido al viento
    rho_bp[-1,:] = Jplatep/(movp*np.sqrt(np.abs((Eyyi[-1,:]+windy[-1,:]/movp)**2))) # se postula que la dirección de E será solamente vertical
    rho_bn[-1,:] = Jplaten/(movn*np.sqrt(np.abs((Eyyi[-1,:]+windy[-1,:]/movn)**2)))
    rho1 = rho_ini.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
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
    return rho1, rho1p, rho1n
    

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
def algoritmo_V_rho(V, rho1, dx, dy, fd_point, fd_value, max_iter, Vb):
    f_rhs = funcion_f(rho1)
    V_b =  Vb.copy()
    Vi = V.copy()
    for iteration in range(max_iter):
        Vold = Vi.copy()
        Vi = update_v(V, f_rhs, dx,dy,fixed_point=fd_point, fixed_value=fd_value, V_boundary=V_b, Exx=Exx,Eyy=Eyy)
        condicion, diff = convergencia(V, Vold, 0.01)
        if condicion:
            #print(f"Convergencia alcanzada para V en la iteración {iteration}")
            break
        V = Vi
    return Vi

# ALGORITMO PRINCIPAL
# Se calcula primero el potencial incial en base a CSM
## Cálculo de potencial eléctrico inicial en base a la carga calculada por CSM
Vmi,Vmi2 = potencial_electrostático(fixed_point, fixed_value, X, Y, [Rp,Rn], 0, carga=Q)

##### Resolución ecuación de continuidad
Jpp = 30.744*10**(-8) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
Jp = 13*10**(-8)
Jpn = Jp - Jpp
# Condiciones de borde
#rho_ip, rho_in = rhoA(Sx, Jp, np.array([movp,movn]), np.array([Rp,Rn]), np.array([mp,mn]), delta, Jpp)
# Coordenadas de los nodos específicos (índices en la malla)
indice1 = (2, 3)  # Nodo en la posición (x=2, y=3) en términos de índices
indice2 = (7, 8)  # Nodo en la posición (x=7, y=8)

# Obtener las coordenadas físicas de los nodos
pos_central = ((x_coor1+x_coor2)/2, (y_coor1+y_coor2)/2)
nod_central =encuentra_nodos(pos_central[0], pos_central[1])
# Calcular el punto central
Ey = E[nod_central[1],nod_central[0]]
D = np.abs(x_coor2-x_coor1)
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
Ecritp= E_onset(mp, delta, R[0])
Ecritn= E_onset(mn, delta, R[1])
V0p = V_onset(mp, Ecritp, R[0])
V0n = V_onset(mn, Ecritn, R[1])
rho_ip, rho_in = rhoA(Ey, D, epsilon0, V0p, V0n, V, Ecritp, Ecritn, R)
rho_inicial = np.zeros((nodosx, nodosy))
'''
Exxini, Eyyini, Em = calcular_campo_electrico(Vmi, dx, dy) # campo electrostático inicial

rho_inicial[fixed_point[0]] = rho_ip
rho_inicial[fixed_point[1]] = rho_in
rho_inicial[-1,:] = Jpp/(movp*np.sqrt(np.abs((Eyyini[-1,:]+windy[-1,:]/movp)**2)))
rho_inicial[-1,:] += Jpn/(movn*np.sqrt(np.abs((Eyyini[-1,:]-windy[-1,:]/movn)**2)))
#print(Vmi)
Vm = Vmi.copy()
it_global = 1
for n in range(it_global):
    Volder = Vm.copy()
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    if n==0:
        rho_d, rho_p, rho_n = algoritmo_rho_v(Vmi, rho_inicial, dx, dy, windx, windy, max_iter_rho, Jpp, Jpn, [rho_ip, rho_in])
    else:
        rho_d, rho_p, rho_n = algoritmo_rho_v(Vm, rho_d0, dx, dy, windx, windy, max_iter_rho, Jpp, Jpn, [rho_ip, rho_in])
    rho_d0 = rho_p-rho_n # LOS rho_pn son puestos en sus valores absolutos
    Vm = algoritmo_V_rho(Vm, rho_d0, dx, dy, fixed_point, fixed_value, max_iter, Vmi)
    condicion,diff = convergencia(Vm, Volder, 0.01)
    if condicion:
        print(f"Convergencia alcanzada para V en la iteración {n}")
        print(r'Diferencia relativa V y Vold: '+str(diff))
        break

#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
print('Potencial calculado')
##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica

Edefx, Edefy, Edef = calcular_campo_electrico(Vm, dx, dy)
Jplus = rho_p*movp*np.sqrt((Edefx+(windx/movp))**2 + (Edefy+(windy/movp))**2)
Jdiff = rho_n*movn*np.sqrt((Edefx-(windx/movn))**2 + (Edefy-(windy/movn))**2)
Ei = Edef[int((Sy-l)/dy),:] # Magnitud Campo eléctrico a nivel de piso
J = Jplus + Jdiff
Ei = Edef[int((Sy-l)/dy),:] # Magnitud Campo eléctrico a nivel de piso
Ji = J[int((Sy-l)/dy),:] # Densidad de corriente a nivel de piso
Jave = np.mean(Ji)
print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:])/(10**(-9)))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')
#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')


######## GRAFICOS
plt.figure(1)
#mod = np.sqrt(U**2+Ww**2)
plt.quiver(X, Y, U, Ww, E, cmap='plasma', scale_units='xy')
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
#plt.show()


# Crear el gráfico usando SymLogNorm para escala logarítmica simétrica
fig, ax = plt.subplots(num=2)
cmap = plt.get_cmap("coolwarm")  # Mapa de color con tonos fríos y cálidos
# La norma SymLogNorm define la escala logarítmica simétrica alrededor de cero
norm = SymLogNorm(linthresh=0.1, linscale=0.1, vmin=Vmi.min(), vmax=Vmi.max(), base=10)
# Graficar el mapa de colores
c = ax.pcolormesh(X, Y, Vmi, cmap=cmap, norm=norm, shading='auto')
# Agregar una barra de color
cbar = fig.colorbar(c, ax=ax, extend='both', ticks=[Vmi.min(), Vmi.min()/10, Vmi.min()/100, 0, Vmi.max()/100, Vmi.max()/10, Vmi.max()])
#cbar = plt.colorbar()
cbar.set_label(r'Potencial $kV$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Potencial electrostático', fontsize=15)
plt.tight_layout()


plt.figure(3)
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

fig, ax = plt.subplots(num=4)
cmap = plt.get_cmap("coolwarm")  # Mapa de color con tonos fríos y cálidos
# La norma SymLogNorm define la escala logarítmica simétrica alrededor de cero
norm = SymLogNorm(linthresh=0.1, linscale=0.1, vmin=rho_p.min(), vmax=rho_p.max(), base=10)
# Graficar el mapa de colores
c = ax.pcolormesh(X, Y, rho_p, cmap=cmap, norm=norm, shading='auto')
# Agregar una barra de color
cbar = fig.colorbar(c, ax=ax, extend='both', ticks=[rho_d.min(), rho_p.min()/10, rho_p.min()/100, 0, rho_p.max()/100, rho_p.max()/10, rho_p.max()])
#cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title(r'Densidad de carga $\rho_+$',fontsize=15)
plt.tight_layout()

fig, ax = plt.subplots(num=5)
cmap = plt.get_cmap("coolwarm")  # Mapa de color con tonos fríos y cálidos
# La norma SymLogNorm define la escala logarítmica simétrica alrededor de cero
norm = SymLogNorm(linthresh=0.1, linscale=0.1, vmin=rho_n.min(), vmax=rho_n.max(), base=10)
# Graficar el mapa de colores
c = ax.pcolormesh(X, Y, rho_n, cmap=cmap, norm=norm, shading='auto')
# Agregar una barra de color
cbar = fig.colorbar(c, ax=ax, extend='both', ticks=[rho_n.min(), rho_n.min()/10, rho_n.min()/100, 0, rho_n.max()/100, rho_n.max()/10, rho_n.max()])
#cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title(r'Densidad de carga final $\rho_-$',fontsize=15)
plt.tight_layout()


#plt.figure(4)
fig, ax = plt.subplots(num=6)
cmap = plt.get_cmap("coolwarm")  # Mapa de color con tonos fríos y cálidos
# La norma SymLogNorm define la escala logarítmica simétrica alrededor de cero
norm = SymLogNorm(linthresh=0.1, linscale=0.1, vmin=rho_d.min(), vmax=rho_d.max(), base=10)
# Graficar el mapa de colores
c = ax.pcolormesh(X, Y, rho_d, cmap=cmap, norm=norm, shading='auto')
# Agregar una barra de color
cbar = fig.colorbar(c, ax=ax, extend='both', ticks=[rho_d.min(), rho_d.min()/10, rho_d.min()/100, 0, rho_d.max()/100, rho_d.max()/10, rho_d.max()])
#cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga final',fontsize=15)
plt.tight_layout()

########## 
#plt.figure(5)
#plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
#plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto',norm=LogNorm())
fig, ax = plt.subplots(num=7)
cmap = plt.get_cmap("coolwarm")  # Mapa de color con tonos fríos y cálidos
# La norma SymLogNorm define la escala logarítmica simétrica alrededor de cero
norm = SymLogNorm(linthresh=0.1, linscale=0.1, vmin=Vm.min(), vmax=Vm.max(), base=10)
# Graficar el mapa de colores
c = ax.pcolormesh(X, Y, Vm, cmap=cmap, norm=norm, shading='auto')
# Agregar una barra de color
cbar = fig.colorbar(c, ax=ax, extend='both', ticks=[Vm.min(), Vm.min()/10, Vm.min()/100, 0, Vm.max()/100, Vm.max()/10, Vm.max()])
#cbar = plt.colorbar()
cbar.set_label(r'Potencial $kV$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
plt.title('Potencial definitivo',fontsize=15) 
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.tight_layout()
#plt.show()



plt.figure(8)
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


plt.figure(9)
plt.plot(x[10:-10], Ji[10:-10]*(10**9))
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m, $w_x=$'+str(wndx), fontsize=13)
plt.tight_layout()
plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
plt.grid(True)

plt.figure(10)
plt.plot(x[10:-10], Ei[10:-10]/1000)
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+r' m, $w_x=$'+str(wndx), fontsize=13)
plt.tight_layout()
plt.legend([f'$|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV'])
plt.grid(True)

plt.show()
'''