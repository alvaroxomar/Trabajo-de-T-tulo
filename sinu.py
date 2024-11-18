import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# Función que define la trayectoria sinusoidal en el rango de x
def s(x, A=8, k=np.pi/10):
    return A * np.cos(k * x) + 10  # Trayectoria sinusoidal

# Función d(x, y) en la malla
def d(x, y):
    return 10 * np.abs(x * y)  # Ejemplo: d(x, y) = 0.1 * x * y (puedes cambiarlo según el problema)

# Resolver la ecuación de difusión sobre la trayectoria utilizando lambda para la EDO
def rk4_trajectory(x_range, y_range, rho0, s, d, grid_shape):
    # Inicializar la malla para rho(x, y)
    rho_grid = np.zeros(grid_shape)  # Matriz para almacenar los valores de rho

    # Inicializar rho en la trayectoria para cada x en el rango
    rho_curr = rho0  # Asignamos el valor inicial de rho

    # Definir la ecuación diferencial usando lambda
    f = lambda rho, x, y: -rho * (10*rho - d(x, y))/(10**(5))  # Ecuación: d rho/ds = -rho * (rho - d(x, y))

    for i in range(1, len(x_range)):  # Comenzamos desde la segunda posición de x
        x_i = x_range[i-1]  # Coordenada x del nodo anterior
        y_i = s(x_i)        # Coordenada y del nodo anterior
        x_ip1 = x_range[i]  # Coordenada x del nodo actual
        y_ip1 = s(x_ip1)    # Coordenada y del nodo actual
        if -5 <= x_i <= 5:
            # Calcular el paso ds como la distancia euclidiana entre dos puntos consecutivos
            ds = np.sqrt((x_ip1 - x_i)**2 + (y_ip1 - y_i)**2)

            # Resolver la ecuación diferencial usando RK4 a lo largo de la trayectoria
            k1 = ds * f(rho_curr, x_i, y_i)
            k2 = ds * f(rho_curr + k1 / 2, x_i, y_i)
            k3 = ds * f(rho_curr + k2 / 2, x_i, y_i)
            k4 = ds * f(rho_curr + k3, x_i, y_i)
            if x_i == -5.07537688:
                print(f'k1 es {k1},k2 es {k2}, k3 es {k3} y k4 es {k4}')
            # Actualizar rho en el punto
            rho_curr += (k1 + 2*k2 + 2*k3 + k4) / 6
            print(f'rho_curr es {rho_curr}')

            # Asignar el valor de rho en la malla en el punto (x, y) usando las coordenadas reales
            if y_range[0] <= y_ip1 < y_range[-1] and x_range[0] <= x_ip1 < x_range[-1]:
                # Convertir las coordenadas reales de x, y a índices de la malla
                paso_x = x_range[1] - x_range[0]  # Paso entre puntos consecutivos en x
                paso_y = y_range[1] - y_range[0]  # Paso entre puntos consecutivos en y

                # Calculando los índices para x y y usando la fórmula propuesta
                x_index = int((x_ip1 - x_range[0]) / paso_x)
                y_index = int((y_ip1 - y_range[0]) / paso_y)

                # Asegurarse de que los índices estén dentro del rango de la malla
                if 0 <= x_index < grid_shape[0] and 0 <= y_index < grid_shape[1]:
                    rho_grid[x_index, y_index] = rho_curr

    return rho_grid

# Parámetros iniciales
x_range = np.linspace(-10, 10, 200)  # Rango de x de -10 a 10
y_range = np.linspace(0, 20, 200)   # Rango de y de 0 a 20
grid_shape = (len(x_range), len(y_range))  # Dimensiones de la malla (x, y)

rho0 = 10  # Valor inicial de rho

# Resolver la ecuación sobre la malla
rho_grid = rk4_trajectory(x_range, y_range, rho0, s, d, grid_shape)

# Graficar los resultados
plt.figure(figsize=(10, 6))
plt.imshow(rho_grid.T, extent=(-10, 10, 0, 20), origin="lower", aspect="auto", cmap="viridis")
plt.colorbar(label="rho(x, y)")
plt.title("Distribución de rho(x, y) en la malla física")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()


