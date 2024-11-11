import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm 
from scipy.interpolate import griddata
plt.close('all')
################################
## Parámetros genericos
R = 6.35 * 10**(-2) # (m) radio equivalente del conductor
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol = 300000 # (V) voltaje de la línea
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
mov = 7*10**(-4) # (m^2/Vs)
m = 1 # (AD) factor de rugosidad
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
windx = 0 # m/s
windy = 0 # m/s
# Sx=3H
# Sy = 2H
Sx = 60 # (m) media longitud del plano de tierra 
Sy = 70 # (m) altura del área de estudio respecto de tierra
l = 0.3 # (m) distancia desde el suelo a la altura de interés de influencia de campo
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

max_iter_rho = 200
max_iter = 1000

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

### Caso unipolar
nodosx = 150
#nodosy = int(1 + (Sy * (nodosx - 1)) / (2 * Sx))
nodosy = 150
windx = np.ones((nodosx, nodosy))*windx
windy = np.ones((nodosx, nodosy))*windy
## Obtención campo eléctrico
E = np.zeros((nodosx,nodosy))
Exx = np.zeros((nodosx,nodosy))
Eyy = np.zeros((nodosx,nodosy))
x = np.linspace(-Sx,Sx,len(E[0,:]))
y = np.linspace(Sy, 0, len(E[:,0]))
dx = np.abs(x[0]-x[1])
dy = np.abs(y[0]-y[1])
for i in range(len(x)):
    for j in range(len(y)):
        Ex=0
        Ey=0
        Exu=0
        Eyu=0
        for z in range(largo):
            N1x= x[i]-w[z]
            N1y= y[j]-h[z]
            N2y = y[j]+h[z]
            Ex += Q[z]*(K/2)*(N1x/(N1x**2+N1y**2)**(3/2)-N1x/(N1x**2+N2y**2)**(3/2))
            Ey += Q[z]*(K/2)*(N1y/(N1x**2+N1y**2)**(3/2) - N2y/(N1x**2+N2y**2)**(3/2))
        Exx[j][i] = Ex
        Eyy[j][i] = Ey
        E[j][i] = np.sqrt(Ex**2 + Ey**2)

X,Y =np.meshgrid(x,y)
U = Exx
Ww = Eyy
print('Campo eléctrico  electrostático libre iones calculado')
plt.figure(1)
mod = np.sqrt(U**2+Ww**2)
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



##### Resolución ecuación de continuidad
Jp = 3.744*10**(-8) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
# Condiciones de borde
rho_i = (Sx*Jp)*10**(-3)/(np.pi*mov*(100*R)*m*(30*delta + 9*np.sqrt(delta/(100*R)))) # radio está en cm
rho = np.zeros((nodosx,nodosy))*10**(-8)
#posx_conductor = int(Sx/dx)
#posy_conductor = int((Sy - coordenada[0][1])/dy)
posx_conductor = int((x_coor + Sx) / dx)  # Para x_conductor = 0
posy_conductor = int((Sy - y_coor) / dy)
fixed_point = (posy_conductor, posx_conductor)

## Resolución para densidad de carga en el espacio

rho1 =np.zeros((nodosx,nodosy))
rho1[fixed_point] = rho_i
#rho1[-1,:] = Jp/(mov*np.sqrt(np.abs((Exx[-1,:]+windx[-1,:]/mov)**2 + (Eyy[-1,:]+windy[-1,:]/mov)**2)))
rho_p = rho1.copy()
Tol = 10**(-4)
con = 0

def gaussian_2d_mesh(row, x0, y0, sigma_x, sigma_y):
    """
    Genera una distribución gaussiana 2D en una malla de tamaño NxN con
    un valor máximo de 1 en el centro.

    Parámetros:
    -----------
    N : int
        Tamaño de la malla (N x N).
    sigma_x : float
        Desviación estándar de la distribución en el eje x.
    sigma_y : float
        Desviación estándar de la distribución en el eje y.

    Retorna:
    --------
    mesh : np.ndarray
        Malla 2D de tamaño NxN con una distribución gaussiana normalizada.
    """
    gaussian_mesh = np.zeros_like(row)
    gaussian_mesh[:,:] = rho_i*np.exp(-((X-x0)**2 / (2 * sigma_x**2) + (Y-y0)**2 / (2 * sigma_y**2)))
    
    # Calcular la distribución gaussiana en la malla
    #gaussian_mesh = rho_i*np.exp(-((X-x0)**2 / (2 * sigma_x**2) + (Y-y0)**2 / (2 * sigma_y**2)))
    
    # Normalizar para que el valor máximo sea 1
    #gaussian_mesh /= np.max(gaussian_mesh)
    
    return gaussian_mesh

def interpolate_specific_value(data, value_to_interpolate=np.nan):
    """
    Interpola valores en una matriz 2D para los nodos con un valor específico.
    
    Parámetros:
    -----------
    data : np.ndarray
        Matriz 2D con valores reales y valores a interpolar.
    value_to_interpolate : float, optional
        Valor en los nodos que se desea interpolar. Por defecto, np.nan.
    
    Retorna:
    --------
    interpolated_data : np.ndarray
        Matriz 2D con los valores interpolados.
    """
    # Crear una copia de los datos para evitar modificar el original
    data_copy = data.copy()
    
    # Crear una máscara para los puntos válidos (excluyendo el valor especificado)
    yii, xii = np.indices(data.shape)
    if np.isnan(value_to_interpolate):
        valid_points = ~np.isnan(data_copy)
    else:
        valid_points = data_copy != value_to_interpolate
    
    # Extraer los puntos válidos y sus coordenadas
    known_coords = np.array([xii[valid_points], yii[valid_points]]).T
    known_values = data_copy[valid_points]
    
    # Interpolar los valores en los nodos con el valor especificado
    interpolated_data = griddata(
        points=known_coords,
        values=known_values,
        xi=(xii, yii),
        method='linear'
    )
    
    # Mantener los valores originales en los nodos que no se desean interpolar
    if np.isnan(value_to_interpolate):
        data_copy[np.isnan(data_copy)] = interpolated_data[np.isnan(data_copy)]
    else:
        data_copy[data_copy == value_to_interpolate] = interpolated_data[data_copy == value_to_interpolate]
    
    return data_copy
# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
# Implementa interpolación lineal  o cúbica
def update_rho_vectorized(rhoini, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_bound):
    #rho_new = np.copy(rho)
    rho = rhoini.copy()
    rho2 = np.zeros_like(rhoini)
    rho2[fixed_point] = (rho_bound[fixed_point]**2)*(10**(30))  # Mantener la condición en el punto central
    rho2[-1, 1:-1] = (rho_bound[-1,1:-1]**2)*(10**(30))  # Mantener la condición en el borde inferior
    # Calcular las derivadas en el interior (diferencias centrales) usando slicing
    d_rho_dx = (rho[1:-1, 2:] - rho[1:-1, :-2]) / (2 * dx)
    d_rho_dy = (rho[2:, 1:-1] - rho[:-2, 1:-1]) / (2 * dy)
    
    Ewx = Ex + wx/mov
    Ewy = Ey + wy/mov
    # Ecuación diferencial discretizada en el interior
    G = -epsilon0 * (Ewx[1:-1, 1:-1] * d_rho_dx + Ewy[1:-1, 1:-1] * d_rho_dy)
    D1 = np.sqrt(np.abs(G))
    rho2[1:-1, 1:-1] = G*(10**(30))
    #rho[1:-1, 1:-1] = np.where(G >= 0, np.sqrt(G),D1)
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
    rho2[1:-1, 0] = Gi*(10**(30))
    #rho[1:-1, 0] = np.where(Gi >= 0, D2, 0)
    #rho[1:-1, 0] = np.where(Gi >= 0, np.sqrt(Gi), 0) # Borde izquierdo
    rho2[1:-1, -1] = Gd*(10**(30))
    #rho[1:-1, -1] = np.where(Gd >= 0, D3, 0)
    #rho[1:-1, -1] = np.where(Gd >= 0, np.sqrt(Gd), 0) # Borde derecho
    rho2[0, 1:-1] = Gs*(10**(30))
    #rho[0, 1:-1] = np.where(Gs >= 0, D4, 0)
    #rho[0, 1:-1] = np.where(Gs >= 0, np.sqrt(Gs), 0) # Borde superior
    rho2[0, 0] =  Geis*(10**(30))
    #rho[0, 0] = np.where(Geis >= 0, D5, 0) # Esquina superior izquierda
    rho2[0, -1] = Geds*(10**(30))
    #rho[0, -1] = np.where(Geds >= 0, D6, 0) # Esquina superior derecha
    rho2[-1, 0] = Geii*(10**(30))
    #rho[-1, 0] = np.where(Geii >= 0, D7, 0)  # Esquina inferior izquierda
    rho2[-1, -1] = Gedi*(10**(30))
    #rho[-1, -1] = np.where(Gedi >= 0, D8, 0)  # Esquina inferior derecha
    
    # Aplicar condiciones de borde
    rho2[fixed_point] = (rho_bound[fixed_point]**2)*(10**(30)) # Mantener la condición en el punto central
    #rho[-1, 1:-1] = rho_bound[-1,1:-1]  # Mantener la condición en el borde inferior
    return rho2

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

## Cálculo de potencial eléctrico inicial en base a la carga calculada por CSM
Vm = np.zeros((nodosx,nodosy))
## Condiciones de borde

fixed_value =Vol
Vm[fixed_point] = fixed_value
# Calcular la distancia entre cada punto de la malla y el punto (w, h)
Mod_coor = np.sqrt((X - w)**2 + (Y - h)**2)
# Calcular el potencial eléctrico en cada nodo de la malla
Vm = Q * (0.5 * K / Mod_coor)
Vm[fixed_point] = fixed_value
Vm[-1,:] = 0 # borde inferior
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
plt.figure(8)
#plt.figure(figsize=(6, 6))
#plt.contourf(X, Y, Vm, levels=200, cmap='plasma')
plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto',norm=LogNorm())
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

# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
def algoritmo_rho_v(V, dx, dy, windx, windy,max_iter_rho):
    conv = 0
    rho_b =np.zeros((nodosx,nodosy))
    rho1 = gaussian_2d_mesh(rho_b, x[posx_conductor], y[posy_conductor], 10, 6)
    # En un comienzo se tiene el campo eléctrico electrostático
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    # Se define a la malla rho_b con las condiciones de borde
    rho_b[fixed_point] = rho_i
    # COndición de desplazamiento de la densidad de carga debido al viento
    rho_b[-1,1:-1] = (Jp/(mov*np.sqrt(np.abs((Exxi[-1,1:-1]+windx[-1,1:-1]/mov)**2 + (Eyyi[-1,1:-1]+windy[-1,1:-1]/mov)**2)))) 
    #rho1 = rho_b.copy() # Parte con una distribución con ceros y las condiciones de borde, luego se actulizan los valores
    #for iteration in range(max_iter_rho):
    for iteration in range(2):
        #rho_1 = rho1.copy()
        rho2 = update_rho_vectorized(rho1, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_b)
        # Criterio de convergencia
        # Define un pequeño valor para evitar divisiones por cero
        '''
        epsilon = 1e-10
        epsilon_array = np.full_like(rho1, epsilon)
        # Calcula la diferencia relativa usando el máximo entre |rho1|, |rho_0| y epsilon
        diferencia_relativa = np.abs(rho1 - rho_0) / np.maximum.reduce([np.abs(rho1), np.abs(rho_0), epsilon_array])
        # Verifica si todas las diferencias relativas son menores al 10%
        condicion = np.all(diferencia_relativa < 0.01)
        diff = np.linalg.norm(diferencia_relativa,ord=np.inf)
        #print(r'Tolerancia $\rho_1$ = '+str(diff))
        if condicion:
            print(f"Convergencia para rho alcanzada en la iteración {iteration}")
            break
        '''
    #rho2 = interpolate_specific_value(rho1, value_to_interpolate=0)
    #print(r'Tolerancia rho_1 = '+str(diff))
    return rho2, rho1
    

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
            u = apply_boundary_conditions_2D(u, V_boundary, Exx=Exx, Eyy=Eyy, dx=dx, dy=dy)
        
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
        epsilon = 1e-10
        epsilon_array = np.full_like(V, epsilon)
        # Calcula la diferencia relativa usando el máximo entre |rho1|, |rho_0| y epsilon
        diferencia_relativa = np.abs(V - Vold) / np.maximum.reduce([np.abs(V), np.abs(Vold), epsilon_array])
        # Verifica si todas las diferencias relativas son menores al 10%
        condicion = np.all(diferencia_relativa < 0.01)
        diff = np.linalg.norm(diferencia_relativa,ord=np.inf)
        #print(r'Tolerancia $V-relativa$ = '+str(diff))
        if condicion:
            #print(f"Convergencia alcanzada para V en la iteración {iteration}")
            break
    return V

# ALGORITMO PRINCIPAL
# Se calcula primero el potencial incial en base a CSM
for n in range(1):
    Volder = Vm.copy()
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    rho_n2, rho_n1 = algoritmo_rho_v(Vm, dx, dy, windx, windy,max_iter_rho)
    #V_n = Vm.copy()
    '''
    Vm = algoritmo_V_rho(Vm, rho_n, dx, dy, fixed_point, fixed_value, max_iter)
    epsilon = 1e-10
    epsilon_array = np.full_like(Vm, epsilon)
    # Calcula la diferencia relativa usando el máximo entre |rho1|, |rho_0| y epsilon
    diferencia_relativa = np.abs(Vm - Volder) / np.maximum.reduce([np.abs(Vm), np.abs(Volder), epsilon_array])
    # Verifica si todas las diferencias relativas son menores al 10%
    condicion = np.all(diferencia_relativa < 0.001)
    diff = np.linalg.norm(diferencia_relativa,ord=np.inf)
    #print(r'Tolerancia $V-relativa$ = '+str(diff))
    print(r'Diferencia relativa V y Vold: '+str(diff))
    if condicion:
        print(f"Convergencia alcanzada para V en la iteración {n}")
        break
    '''
#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
print('Potencial calculado')
rho_n2 = np.asarray(rho_n2, dtype=float)
rho_n2 = np.nan_to_num(rho_n2, nan=0.0, posinf=0.0, neginf=0.0)
rho_n2 = rho_n2.astype(float)
plt.figure(7)
#plt.pcolormesh(X, Y, rho_n2, cmap='viridis', shading='auto',norm=LogNorm())
#plt.pcolormesh(X, Y, rho_n2, cmap='viridis', shading='auto')
maxval = np.max(np.abs([rho_n2.min(), rho_n2.max()]))
plt.imshow(rho_n2, cmap='RdYlGn', interpolation='nearest', vmin=-maxval, vmax=maxval)
plt.colorbar(label="Valor")
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga final',fontsize=15)
# Añadir una barra de colores para mostrar la escala
#cbar = plt.colorbar()
#cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
#plt.show()

plt.figure(10)
plt.pcolormesh(X, Y, rho_n1, cmap='viridis', shading='auto',norm=LogNorm())
#plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto')
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga final',fontsize=15)
# Añadir una barra de colores para mostrar la escala
cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
plt.show()
##########
''' 
plt.figure(3)
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
'''

##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica

'''
Edefx, Edefy, Edef = calcular_campo_electrico(Vm, dx, dy)
J = rho_n*mov*np.sqrt((Edefx+(windx/mov))**2 + (Edefy+(windy/mov))**2)
Ei = Edef[int((Sy-l)/dy),:] # Magnitud Campo eléctrico a nivel de piso
Ji = J[int((Sy-l)/dy),:] # Densidad de corriente a nivel de piso
Jave = np.mean(Ji)



print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:])/(10**(-9)))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')
#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')
plt.figure(4)
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


plt.figure(5)
plt.plot(x[10:-10], Ji[10:-10]*(10**9))
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
plt.grid(True)

plt.figure(6)
plt.plot(x[10:-10], Ei[10:-10]/1000)
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.legend([f'$|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV'])
plt.grid(True)
'''
plt.show()