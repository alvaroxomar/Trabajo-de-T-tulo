import numpy as np
import matplotlib.pyplot as plt

# Parámetros de la malla
n = 64  # Número de puntos en cada dirección
L = 1.0  # Tamaño del dominio
h = L / (n - 1)  # Tamaño del paso en la malla

# Crear la malla en 2D
x = np.linspace(0, L, n)
y = np.linspace(0, L, n)
X, Y = np.meshgrid(x, y)

# Función fuente (f)
def f(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)

# Inicializar la función u y el lado derecho f en la malla
u = np.zeros((n, n))
f_rhs = f(X, Y) * h**2  # f(x, y) multiplicado por h^2 (parte de la discretización)

def apply_laplacian(u):
    laplacian = np.zeros_like(u)
    laplacian[1:-1, 1:-1] = (
        u[0:-2, 1:-1] + u[2:, 1:-1] + u[1:-1, 0:-2] + u[1:-1, 2:]
        - 4 * u[1:-1, 1:-1]
    )
    return laplacian

# Método Jacobi para resolver el sistema de ecuaciones
def jacobi_step(u, f_rhs):
    unew = np.copy(u)
    unew[1:-1, 1:-1] = 0.25 * (
        u[0:-2, 1:-1] + u[2:, 1:-1] + u[1:-1, 0:-2] + u[1:-1, 2:]
        - f_rhs[1:-1, 1:-1]
    )
    return unew

def restrict(u):
    return 0.25 * (u[0:-1:2, 0:-1:2] + u[1::2, 0:-1:2] + u[0:-1:2, 1::2] + u[1::2, 1::2])

#def prolong(u):
#    v = np.zeros((2 * u.shape[0] - 1, 2 * u.shape[1] - 1))
#    v[::2, ::2] = u
#    v[1::2, ::2] = 0.5 * (u[:-1, :] + u[1:, :])
#    v[::2, 1::2] = 0.5 * (u[:, :-1] + u[:, 1:])
#    v[1::2, 1::2] = 0.25 * (u[:-1, :-1] + u[:-1, 1:] + u[1:, :-1] + u[1:, 1:])
#    return v
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



def v_cycle(u, f_rhs, n_cycles=1):
    if u.shape[0] <= 3:  # Resolver directamente cuando la grilla es muy pequeña
        return jacobi_step(u, f_rhs)

    for _ in range(n_cycles): 
        # Pre-suavizado
        u = jacobi_step(u, f_rhs)
        print(r'dimension u inicial: '+str(u.shape))
        
        # Calcular el residuo
        res = f_rhs - apply_laplacian(u)
        print(r'dimension res: '+str(res.shape))
        
        # Restricción (proyectar a una grilla más gruesa)
        res_coarse = restrict(res)
        print(r'dimension res coarse: '+str(res_coarse.shape))
        u_coarse = np.zeros_like(res_coarse)
        
        # Resolver en la grilla gruesa
        u_correction_coarse = v_cycle(u_coarse, res_coarse)
        
        # Prolongación (interpolar a la grilla fina)
        u_correction = prolong(u_correction_coarse)
        print(r'Dimensión u_prolongado: '+str(u_correction.shape))
        
        # Ajustar la corrección al tamaño de 'u' si es necesario
        if u_correction.shape != u.shape:
            # Ajustar el tamaño de u_correction mediante interpolación adecuada
            diff_x = u.shape[0] - u_correction.shape[0]
            diff_y = u.shape[1] - u_correction.shape[1]

            if diff_x > 0:
                u_correction = np.pad(u_correction, ((0, diff_x), (0, 0)), mode='constant')
            if diff_y > 0:
                u_correction = np.pad(u_correction, ((0, 0), (0, diff_y)), mode='constant')

            print(r'dimension u_correction ajustada: '+str(u_correction.shape))
        
        # Actualizar la solución
        u += u_correction
        
        # Post-suavizado
        u = jacobi_step(u, f_rhs)
    
    return u

# Resolver usando FMG
u = v_cycle(u, f_rhs)

# Visualizar la solución
plt.contourf(X, Y, u, 20, cmap='RdBu_r')
plt.colorbar(label='u(x, y)')
plt.title('Solución de la ecuación de Poisson (2D) usando FMG')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
