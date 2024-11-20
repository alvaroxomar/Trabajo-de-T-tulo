import numpy as np
import matplotlib.pyplot as plt
plt.close('all')


# Dimensiones de la malla y configuración
n = 10  # Tamaño de la malla (10x10)
cy, cx = 4, 4  # Coordenadas del nodo pivote (cy, cx)
num_iteraciones = 8  # Número de iteraciones

# Inicializar la malla con ceros y asignar el valor fijo en (cy, cx)
matriz = np.zeros((n, n))
matriz[cy, cx] = 10  # Valor inicial fijo en el nodo pivote

# Lista para guardar la matriz en cada iteración
iteraciones = [matriz.copy()]

# Aplicar las fórmulas diferenciadas en cada cuadrante en iteraciones
for _ in range(num_iteraciones):
    nueva_matriz = matriz.copy()
    
    # Primer cuadrante (Inferior derecho): Promedio de arriba e izquierda, excluyendo el borde inferior
    nueva_matriz[cy:-1, cx:] = (matriz[cy-1:-2, cx:] + matriz[cy:-1, cx-1:-1]) / 2

    # Segundo cuadrante (Superior derecho): Suma de los cuadrados de abajo e izquierda, excluyendo el borde inferior
    nueva_matriz[:cy+1, cx:] = matriz[1:cy+2, cx:]**2 + matriz[:cy+1, cx-1:-1]**2

    # Tercer cuadrante (Superior izquierdo): Raíz cuadrada de la suma de abajo y derecha, excluyendo el borde inferior
    nueva_matriz[:cy+1, :cx+1] = np.sqrt(np.abs(matriz[1:cy+2, :cx+1] + matriz[:cy+1, 1:cx+2]))

    # Cuarto cuadrante (Inferior izquierdo): Suma directa de arriba y derecha, excluyendo el borde inferior
    nueva_matriz[cy:-1, :cx+1] = matriz[cy-1:-2, :cx+1] + matriz[cy:-1, 1:cx+2]

    # Mantener el valor fijo en el nodo pivote (cy, cx) en cada iteración
    nueva_matriz[cy, cx] = 10
    matriz = nueva_matriz
    iteraciones.append(matriz.copy())

# Crear gráficos para mostrar la evolución de cada iteración
fig, axes = plt.subplots(1, num_iteraciones + 1, figsize=(15, 4))
for i, ax in enumerate(axes):
    im = ax.imshow(iteraciones[i], cmap='viridis', vmin=0, vmax=10)
    ax.set_title(f"Iteración {i}")
    ax.axis("off")

# Agregar barra de color
fig.colorbar(im, ax=axes, orientation="vertical", fraction=0.02, pad=0.04)
plt.show()


