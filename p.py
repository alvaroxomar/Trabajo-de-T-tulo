import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
# Parámetros y constantes
nodos = 100  # Tamaño de la malla
epsilon0 = 8.85e-12  # Permitividad eléctrica
dx = 1.0  # Diferencial en x
dy = 1.0  # Diferencial en y
max_iter = 2000  # Número máximo de iteraciones
Tol = 1e-5  # Tolerancia de convergencia

# Inicialización de las matrices para rho, Ex, y Ey
rho = np.zeros((nodos, nodos))  # Matriz de densidad inicial
Ex = np.random.rand(nodos, nodos)  # Ejemplo de campo eléctrico en x
Ey = np.random.rand(nodos, nodos)  # Ejemplo de campo eléctrico en y

# Condiciones iniciales en los bordes
fixed_point = (50, 50)  # Punto central donde se impone una condición de borde
rho[fixed_point] = 1.0  # Condición de borde en el punto central
rho[-1, :] = 0.5  # Condición de borde en el borde inferior

# Función para actualizar rho utilizando operaciones vectorizadas
def update_rho_vectorized(rho, Ex, Ey, dx, dy, epsilon0):
    rho_new = np.copy(rho)
    
    # Calcular las derivadas en el interior (diferencias centrales) usando slicing
    d_rho_dx = (rho[1:-1, 2:] - rho[1:-1, :-2]) / (2 * dx)
    d_rho_dy = (rho[2:, 1:-1] - rho[:-2, 1:-1]) / (2 * dy)
    
    # Ecuación diferencial discretizada en el interior
    rho_new[1:-1, 1:-1] = np.sqrt(np.abs(-epsilon0 * (Ex[1:-1, 1:-1] * d_rho_dx + Ey[1:-1, 1:-1] * d_rho_dy)))
    
    # Condiciones en los bordes
    rho_new[1:-1, 0] = np.sqrt(np.abs(-epsilon0 * (Ex[1:-1, 0] * (rho[1:-1, 1] - rho[1:-1, 0]) / dx + Ey[1:-1, 0] * (rho[2:, 0] - rho[:-2, 0]) / (2 * dy))))
    rho_new[1:-1, -1] = np.sqrt(np.abs(-epsilon0 * (Ex[1:-1, -1] * (rho[1:-1, -1] - rho[1:-1, -2]) / dx + Ey[1:-1, -1] * (rho[2:, -1] - rho[:-2, -1]) / (2 * dy))))
    rho_new[0, 1:-1] = np.sqrt(np.abs(-epsilon0 * (Ex[0, 1:-1] * (rho[0, 2:] - rho[0, :-2]) / (2 * dx) + Ey[0, 1:-1] * (rho[1, 1:-1] - rho[0, 1:-1]) / dy)))
    
    # Aplicar condiciones de borde
    rho_new[fixed_point] = 1.0  # Mantener la condición en el punto central
    rho_new[-1, :] = 0.5  # Mantener la condición en el borde inferior
    
    return rho_new

# Algoritmo iterativo para resolver rho
conv = 0
for iteration in range(max_iter):
    rho_new = update_rho_vectorized(rho, Ex, Ey, dx, dy, epsilon0)
    # Criterio de convergencia
    diff = np.linalg.norm(rho_new - rho)
    #if diff < Tol:
    #    print(f"Convergencia alcanzada en la iteración {iteration}")
    #    break
    
    rho = rho_new

# Mostrar el resultado final
print("Solución final para rho:")

# Ploteo de la solución final con una barra de colores
plt.figure(figsize=(8, 6))
plt.imshow(rho, cmap='viridis', origin='lower')
plt.colorbar(label='Densidad de carga ρ')
plt.title('Distribución de densidad de carga (ρ)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
