import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm 
plt.close('all')
################################
## Parámetros genericos
R = 6.35 * 10**(-2) # (m) radio equivalente del conductor
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol = 500000 # (V) voltaje de la línea
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
mov = 7*10**(-5) # (m^2/Vs)
m = 1 # (AD) factor de rugosidad
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
windx = 0 # m/s
windy = 0 # m/s

Sx = 30 # (m) media longitud del plano de tierra 
Sy = 50 # (m) altura del área de estudio respecto de tierra
l = 0.3 # (m) distancia desde el suelo a la altura de interés de influencia de campo
## Definición coordenadas conductores caso bipolar
#coordenada = [(6,7), (-6,7)] # (m) coordenadas cargas en posición de los conductores
x_coor = 0
y_coor = 23
coordenada = [(x_coor,y_coor)]
coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
#Hor = np.linspace(-10,10,100) # (m) Ancho total del área de estudio
#Ver = 1 # (m) Altura donde interesa medir el campo eléctrico
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
nodosx = 100
#nodosy = int(1 + (Sy * (nodosx - 1)) / (2 * Sx))
nodosy = 100
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
plt.tight_layout()
#plt.show()



##### Resolución ecuación de continuidad
Jp = 1.744*10**(-8) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
# Condiciones de borde
rho_i = (Sx*Jp)*10**(-3)/(np.pi*mov*(100*R)*m*(30*delta + 9*np.sqrt(delta/(100*R)))) # radio está en cm
rho = np.zeros((nodosx,nodosy))*10**(-8)
#posx_conductor = int(Sx/dx)
#posy_conductor = int((Sy - coordenada[0][1])/dy)
posx_conductor = int((x_coor + Sx) / dx)  # Para x_conductor = 0
posy_conductor = int((Sy - y_coor) / dy)
fixed_point = (posy_conductor, posx_conductor)

## Resolución para densidad de carga en el espacio
'''
con = 1
c = 0
rho_n = rho.copy() # se copia el primer valor
#while np.all(con > 10**(-7)):
#    print(rho_n)
while True:
    if c<200:
        for i in range(1, len(rho[:,0])-1):
            for j in range(1, len(rho[0,:])-1):
                if i == posy_conductor and j == posx_conductor:
                    rho_n[i][j] = rho_i
                else:
                    #alpha = 0.5*(np.abs(Exx[i][j])/dx + np.abs(Eyy[i][j])/dy)*epsilon0
                    #beta = epsilon0*(np.abs(Eyy[i][j])*rho_n[i-1][j]/dy + np.abs(Exx[i][j])*rho_n[i][j-1]/dx)
                    #alpha = 0.5*(Exx[i][j]/dx + Eyy[i][j]/dy)*epsilon0
                    #beta = epsilon0*(Eyy[i][j]*rho_n[i-1][j]/dy + Exx[i][j]*rho_n[i][j-1]/dx)
                    #rho_n[i][j] = -alpha + np.sqrt(np.abs(alpha**2 + beta)) # Se toma solución positiva dado que se asume que todas los iones son de carga positiva
                    a = rho_n[i][j+1]- rho_n[i][j-1]
                    b = rho_n[i+1][j]- rho_n[i-1][j]
                    rho_n[i][j] = np.sqrt(0.5*epsilon0*(np.abs(Exx[i][j]*(-a))/dx + np.abs(Eyy[i][j]*(-b))/dy))
                    #if alpha**2 + beta < 0:
                    #    print(alpha**2 + beta)
                    #    break
        #print(r'rho_n: '+str(rho_n))
        #print(r'rho: '+str(rho))
        #print(r'iteracion: '+str(c))
        con = np.abs(rho_n - rho)
        #print(r'conv: '+str(con))
        c +=1 
        rho = rho_n.copy()
    else:
        break
'''
rho1 =np.zeros((nodosx,nodosy))
rho1[fixed_point] = rho_i
rho1[-1,:] = Jp/(mov*np.sqrt(np.abs((Exx[-1,:]+windx[-1,:]/mov)**2 + (Eyy[-1,:]+windy[-1,:]/mov)**2)))
rho_p = rho1.copy()
Tol = 10**(-4)
con = 0
max_iter_rho = 500
max_iter = 1
'''
while True:
    if con < max_iter:
        b = rho_p[:-2, 1:-1] - rho_p[2:, 1:-1]
        a = rho_p[1:-1, :-2] - rho_p[1:-1, 2:]
        rho_p[1:-1,1:-1] = np.sqrt(0.5*epsilon0*(np.abs(Exx[1:-1,1:-1]*-(a)/dx + Eyy[1:-1,1:-1]*(-b))/dy))
        if rho_p[fixed_point] is not None:
            rho_p[fixed_point] = rho_i
            rho_p[-1,:] = np.abs(Jp/(mov*E[-1,:]))
        con += 1
        div = np.abs(rho - rho_p)
        rho = rho_p
        #if np.linalg.norm(div) < Tol:
        #    break
    else:
        break
'''
# Función para actualizar rho utilizando operaciones vectorizadas
# utiliza Gaus-Seidel donde utiliza el mismo valor actualizado para hacer la operación en el nodo adyacente
def update_rho_vectorized(rho, Ex, Ey, dx, dy, epsilon0, wx, wy, rho_bound):
    #rho_new = np.copy(rho)
    
    # Calcular las derivadas en el interior (diferencias centrales) usando slicing
    d_rho_dx = (rho[1:-1, 2:] - rho[1:-1, :-2]) / (2 * dx)
    d_rho_dy = (rho[2:, 1:-1] - rho[:-2, 1:-1]) / (2 * dy)
    
    Ewx = Ex + wx/mov
    Ewy = Ey + wy/mov

    # Ecuación diferencial discretizada en el interior
    rho[1:-1, 1:-1] = np.sqrt(np.abs(-epsilon0 * (Ewx[1:-1, 1:-1] * d_rho_dx + Ewy[1:-1, 1:-1] * d_rho_dy)))
    
    # Condiciones en los bordes
    rho[1:-1, 0] = np.sqrt(np.abs(-epsilon0 * (Ewx[1:-1, 0] * (rho[1:-1, 1] - rho[1:-1, 0]) / dx + Ewy[1:-1, 0] * (rho[2:, 0] - rho[:-2, 0]) / (2 * dy))))
    rho[1:-1, -1] = np.sqrt(np.abs(-epsilon0 * (Ewx[1:-1, -1] * (rho[1:-1, -1] - rho[1:-1, -2]) / dx + Ewy[1:-1, -1] * (rho[2:, -1] - rho[:-2, -1]) / (2 * dy))))
    rho[0, 1:-1] = np.sqrt(np.abs(-epsilon0 * (Ewx[0, 1:-1] * (rho[0, 2:] - rho[0, :-2]) / (2 * dx) + Ewy[0, 1:-1] * (rho[1, 1:-1] - rho[0, 1:-1]) / dy)))
    rho[0, 0] = np.sqrt(np.abs(-epsilon0 * (Ewx[0, 0] * (rho[0, 1] - rho[0, 0]) / dx + Ewy[0, 0] * (rho[1, 0] - rho[0, 0]) / dy)))
    rho[0, -1] = np.sqrt(np.abs(-epsilon0 * (Ewx[0, -1] * (rho[0, -1] - rho[0, -2]) / dx + Ewy[0, -1] * (rho[1, -1] - rho[0, -1]) / dy)))
    
    # Aplicar condiciones de borde
    rho[fixed_point] = rho_bound[fixed_point]  # Mantener la condición en el punto central
    rho[-1, :] = rho_bound[-1,:]  # Mantener la condición en el borde inferior
    
    return rho

## Ocupa el negativo del gradiente del potencial calculado
def calcular_campo_electrico(V, dx, dy):
    # Inicializar las matrices para los componentes del campo eléctrico
    Ex = np.zeros_like(V)
    Ey = np.zeros_like(V)
    
    # Calcular las derivadas centrales (excepto en los bordes)
    Ey[1:-1, :] = (V[2:, :] - V[:-2, :]) / (2 * dy)
    Ex[:, 1:-1] = (V[:, 2:] - V[:, :-2]) / (2 * dx)
    
    # Bordes: Aproximaciones hacia adelante y hacia atrás
    Ey[0, :] = (V[1, :] - V[0, :]) / dy  # Borde inferior
    Ey[-1, :] = (V[-1, :] - V[-2, :]) / dy  # Borde superior
    Ex[:, 0] = (V[:, 1] - V[:, 0]) / dx  # Borde izquierdo
    Ex[:, -1] = (V[:, -1] - V[:, -2]) / dx  # Borde derecho
    
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
plt.contourf(X, Y, Vm, levels=200, cmap='plasma')
plt.colorbar(label='Potencial eléctrico (V)')
plt.title('Malla de Potencial Eléctrico inicial libre de iones')
plt.xlabel('x')
plt.ylabel('y')

# Algoritmo iterativo para resolver rho
# Em base a la malla de potencial V que exista,calcula rho donde los campos eléctricos 
# son calculados en base a -\nabla V
def algoritmo_rho_v(V, dx, dy, windx, windy,max_iter_rho):
    conv = 0
    rho_b =np.zeros((nodosx,nodosy))
    Exxi, Eyyi, Em = calcular_campo_electrico(V, dx, dy)
    # Se define a la malla rho_b con las condiciones de borde
    rho_b[fixed_point] = rho_i
    rho_b[-1,:] = Jp/(mov*np.sqrt(np.abs((Exxi[-1,:]+windx[-1,:]/mov)**2 + (Eyyi[-1,:]+windy[-1,:]/mov)**2)))
    rho1 = rho_b.copy()
    for iteration in range(max_iter_rho):
        rho_0 = rho1.copy()
        rho1 = update_rho_vectorized(rho1, Exxi, Eyyi, dx, dy, epsilon0, windx, windy, rho_b)
        # Criterio de convergencia
        diff = np.linalg.norm((rho1 - rho_0)/rho1,ord=np.inf)
        print(r'Tolerancia $\rho_1$ = '+str(diff))
        if diff < Tol:
            print(f"Convergencia alcanzada en la iteración {iteration}")
            break
    print(r'Tolerancia $\rho_1$ = '+str(diff))
    return rho1
    

print('Densidad de carga preliminar alrededor del conductor calculado')
# Mostrar la matriz en una gráfica con un mapa de colores
'''
plt.figure(2)
plt.pcolormesh(X, Y, rho1, cmap='viridis', shading='auto',norm=LogNorm())
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga preliminar',fontsize=15)
# Añadir una barra de colores para mostrar la escala
cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
plt.show()
'''

#### Calculo del potencial con FMG
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

def restrict(u):
    return 0.25 * (u[0:-1:2, 0:-1:2] + u[1::2, 0:-1:2] + u[0:-1:2, 1::2] + u[1::2, 1::2])


def prolong(u):
    # Dimensiones de la grilla fina
    fine = np.zeros((2 * u.shape[0] - 1, 2 * u.shape[1] - 1))
    
    # Asignar valores a las posiciones originales
    fine[::2, ::2] = u  # Posiciones originales (los puntos gruesos)
    
    # Interpolar en las posiciones intermedias entre las filas
    fine[1::2, ::2] = 0.5 * (u[:-1, :] + u[1:, :])
    
    # Interpolar en las posiciones intermedias entre las columnas
    fine[::2, 1::2] = 0.5 * (u[:, :-1] + u[:, 1:])
    
    # Interpolar en los puntos centrales (diagonales)
    fine[1::2, 1::2] = 0.25 * (u[:-1, :-1] + u[:-1, 1:] + u[1:, :-1] + u[1:, 1:])
    
    return fine


# Definir las condiciones de borde usando una matriz 2D
def apply_boundary_conditions_2D(u, V_boundary, Exx=None, Eyy=None, dx=None, dy=None):
    """
    Aplica las condiciones de borde a la malla u utilizando la matriz V_boundary.
    También aplica condiciones de Neumann si se proporcionan Exx, Eyy, dx y dy.
    """
    # Borde izquierdo (Neumann si Exx está definido)
    if Exx is not None and dx is not None:
        u[:, 0] = u[:, 1] + dx * (Exx[:, 0]+windx[:,0]/mov)
    else:
        u[:, 0] = V_boundary[:, 0]  # Dirichlet

    # Borde derecho (Neumann si Exx está definido)
    if Exx is not None and dx is not None:
        u[:, -1] = u[:, -2] - dx * (Exx[:, -1]+windx[:,-1]/mov)
    else:
        u[:, -1] = V_boundary[:, -1]  # Dirichlet

    # Borde superior (Neumann si Eyy está definido)
    if Eyy is not None and dy is not None:
        u[0, :] = u[1, :] + dy * (Eyy[0, :]+windy[0,:]/mov)
    else:
        u[0, :] = V_boundary[0, :]  # Dirichlet

    # Borde inferior (siempre Dirichlet)
    u[-1, :] = V_boundary[-1, :]  

    return u


# Función modificada de v_cycle
def v_cycle(u, f_rhs, dx, dy, n_cycles=1, fixed_point=None, fixed_value=None, V_boundary=None, Exx=None, Eyy=None):
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
        
        # Calcular el residuo
        res = f_rhs - apply_laplacian_asymmetric(u, dx, dy)
        
        # Restricción (proyectar a una grilla más gruesa)
        res_coarse = restrict(res)
        u_coarse = np.zeros_like(res_coarse)
        
        # Resolver en la grilla gruesa
        u_correction_coarse = v_cycle(u_coarse, res_coarse, 2*dx, 2*dy, fixed_point=None, fixed_value=None)
        
        # Prolongación (interpolar a la grilla fina)
        u_correction = prolong(u_correction_coarse)
        
        # Ajustar la corrección al tamaño de 'u' si es necesario
        if u_correction.shape != u.shape:
            diff_x = u.shape[0] - u_correction.shape[0]
            diff_y = u.shape[1] - u_correction.shape[1]

            if diff_x > 0:
                u_correction = np.pad(u_correction, ((0, diff_x), (0, 0)), mode='constant')
            if diff_y > 0:
                u_correction = np.pad(u_correction, ((0, 0), (0, diff_y)), mode='constant')
        
        # Actualizar la solución
        u += u_correction
        
        # Imponer condiciones de borde nuevamente tras corrección
        if V_boundary is not None:
            u = apply_boundary_conditions_2D(u, V_boundary, Exx=Exx, Eyy=Eyy, dx=dx, dy=dy)
        
        # Imponer condición fija tras la corrección
        if fixed_point is not None and fixed_value is not None:
            u[fixed_point] = fixed_value
        
        # Post-suavizado
        u = jacobi_step(u, f_rhs, dx, dy)
        
        # Imponer condiciones de borde después del post-suavizado
        if V_boundary is not None:
            u = apply_boundary_conditions_2D(u, V_boundary, Exx=Exx, Eyy=Eyy, dx=dx, dy=dy)
        
        # Imponer condición fija en el punto después del post-suavizado
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
        V = v_cycle(V, f_rhs, dx,dy,fixed_point=fixed_point, fixed_value=fixed_value, V_boundary=V_b, Exx=Exx,Eyy=Eyy)
        diff = np.linalg.norm((V-Vold)/V,ord=np.inf)
        if  diff < Tol:
            print(r'COnvergencia potencial: '+str(diff))
    return V

# ALGORITMO PRINCIPAL
# Se calcula primero el potencial incial en base a CSM
for n in range(1):
    # Parte estimando rho inicial en base a Vm inicial y rho1 que contiene las condciones de borde para rho
    # Luego, con el nuevo rho_n, calcula el potencial en el espacio en base al potencial anterior Vm
    rho_n = algoritmo_rho_v(Vm, dx, dy, windx, windy,max_iter_rho)
    #V_n = Vm.copy()
    Vm = algoritmo_V_rho(Vm, rho_n, dx, dy, fixed_point, fixed_value, max_iter)

#for iteration in range(max_iter):
#rho_nuevo = densidad_voltaje(Vm, rho1)
#f_rhs = funcion_f(rho_nuevo)
print('Potencial calculado')

plt.figure(7)
plt.pcolormesh(X, Y, rho_n, cmap='viridis', shading='auto',norm=LogNorm())
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga final',fontsize=15)
# Añadir una barra de colores para mostrar la escala
cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
plt.show()
########## 
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


##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica


Edefx, Edefy, Edef = calcular_campo_electrico(Vm, dx, dy)
J = rho_n*mov*np.sqrt((Edefx+(windx/mov))**2 + (Edefy+(windy/mov))**2)
Ei = Edef[int((Sy-l)/dy),:] # Magnitud Campo eléctrico a nivel de piso
Ji = J[int((Sy-l)/dy),:] # Densidad de corriente a nivel de piso
Jave = np.mean(Ji)



print(r'Jp promedio calculado a l=0 m: '+str(np.mean(J[-1,:]))+' nA/m^2, y a l='+str(l)+' m, Jp ='+str(Jave*(10**9))+' nA/m^2')
print(r'Jp promedio propuesto: '+str(Jp*(10**9))+' nA/m^2')
#Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
print('Campo eléctrico y densidad de corriente iónica ya calculados')
plt.figure(4)
#plt.quiver(Xe, Ye, Edefx[1:-1, 1:-1], Edefy[1:-1, 1:-1], Edef[1:-1, 1:-1], cmap='plasma')
plt.quiver(X, Y, -Edefx, Edefy, Edef, cmap='plasma', scale_units='xy')
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
plt.plot(x, Ji*(10**9))
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Densidad de corriente iónica ($nA/m^2$)',fontsize=11)
plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.legend([f'$J_p$ = {str(np.round(Jave*(10**9),3))} $nA/m^2$'])
plt.grid(True)

plt.figure(6)
plt.plot(x, Ei/1000)
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Campo eléctrico (kV/m)',fontsize=11)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.legend([f'$|E|_a$ = {str(np.round(np.mean(Ei/1000),3))} kV'])
plt.grid(True)

plt.show()