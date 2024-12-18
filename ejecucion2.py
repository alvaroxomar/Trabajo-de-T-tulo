import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
import numpy as np

# Estado de gráficos
estados_graficos = {"Gráfico 1": False, "Gráfico 2": False, "Gráfico 3": False, "Gráfico 4": False,
                     "Gráfico 5": False, "Gráfico 6": False, "Gráfico 7": False, "Gráfico 8": False}

def alternar_grafico(nombre):
    estados_graficos[nombre] = not estados_graficos[nombre]
    mensaje = "seleccionado" if estados_graficos[nombre] else "deseleccionado"
    messagebox.showinfo(nombre, f"Se ha {mensaje} {nombre}.")

# Botón Bipolar para activar/desactivar las casillas extra
def activar_bipolar():
    # Activar o desactivar las casillas extras
    if boton_bipolar_var.get():  # Si se selecciona el botón Bipolar
        entradas["Voltaje 2(kV)"].config(state="normal")
        entradas["Área 2 (mm²)"].config(state="normal")
        entradas["Posición en x 2 (m)"].config(state="normal")
        entradas["Posición en y 2 (m)"].config(state="normal")
        entradas["factor conductor 2"].config(state="normal")
    else:  # Si el botón Bipolar está desmarcado
        entradas["Voltaje 2(kV)"].config(state="disabled")
        entradas["Área 2 (mm²)"].config(state="disabled")
        entradas["Posición en x 2 (m)"].config(state="disabled")
        entradas["Posición en y 2 (m)"].config(state="disabled")
        entradas["factor conductor 2"].config(state="disabled")

        # Eliminar el contenido actual de las casillas
        entrada_posicion_x2.delete(0, tk.END)  
        entrada_posicion_y2.delete(0, tk.END)  
        
        # Insertar valores vacíos (o None si prefieres)
        entrada_posicion_x2.insert(0, "")  
        entrada_posicion_y2.insert(0, "")  
        
        # Establecer x2 y y2 como None para que no se grafiquen
        #global x2, y2  # Si las variables x2 y y2 son globales
        #x2, y2 = None, None
        
        # Si ya había un gráfico previo, eliminar los puntos asociados a x2, y2
        plt.close('all')  # Cierra el gráfico previo
        # Redibujar solo los puntos x, y (sin x2, y2)
        validar_campos_y_graficar()  # Esto vuelve a graficar los puntos con las nuevas configuraciones



def validar_campos_y_ejecutar():
    # Obtener valor de cantidad
    cantidad = entrada_cantidad.get()
    try:
        cantidad = int(cantidad)
    except ValueError:
        cantidad = 0

    for nombre, entrada in entradas.items():
        if nombre == "Separación (m)" and cantidad <= 1:
            continue  # Ignorar esta casilla
        if nombre == "Rugosidad terreno":
            continue

        valor = entrada.get().strip()  # Obtener valor sin espacios adicionales

        # Validar campo vacío
        if not valor:
            messagebox.showwarning("Advertencia", f"No pueden haber casillas vacías (revisar: {nombre}).")
            return

        # Validar "Modo (str)" para aceptar solo "uniforme" o "gradiente"
        if nombre == "Modo (str)":
            if valor.isdigit():
                messagebox.showwarning(
                    "Advertencia", f"El valor en '{nombre}' no puede ser numérico. Debe ser 'uniforme' o 'gradiente'."
                )
                return
            if valor not in ["uniforme", "gradiente"]:
                messagebox.showwarning(
                    "Advertencia", f"El valor ingresado en '{nombre}' debe ser 'uniforme' o 'gradiente'."
                )
                return
            continue  # Saltar el resto de la validación para este campo

        # Validar que los demás campos sean números válidos
        try:
            float(valor)
        except ValueError:
            messagebox.showwarning("Advertencia", f"El valor ingresado en '{nombre}' debe ser un número.")
            return

    # Validar campos de dimensionamiento
    if not entrada_ancho_total.get() or not entrada_altura_total.get():
        messagebox.showwarning("Advertencia", "Debe ingresar valores numéricos en las casillas de dimensionamiento.")
        return

    try:
        float(entrada_ancho_total.get())
        float(entrada_altura_total.get())
    except ValueError:
        messagebox.showwarning("Advertencia", "Los valores de ancho y altura deben ser números.")
        return

    messagebox.showinfo("Validación exitosa", "Todos los campos han sido completados correctamente.")

def validar_campos_y_graficar():
    try:
        # Validar datos de dimensionamiento
        ancho_mitad = float(entrada_ancho_total.get())
        altura_total = float(entrada_altura_total.get())
        
        if ancho_mitad <= 0 or altura_total <= 0:
            raise ValueError("El ancho y la altura deben ser mayores que 0.")

        # Validar datos de conductores
        cantidad = int(entrada_cantidad.get())
        if cantidad <= 0:
            raise ValueError("La cantidad debe ser mayor que 0.")
        
        x = float(entrada_posicion_x.get())
        y = float(entrada_posicion_y.get())

        # Obtener las posiciones x2 y y2
        x2e = entrada_posicion_x2.get()
        y2e = entrada_posicion_y2.get()
        x2 = float(x2e) if x2e.replace('.', '', 1).isdigit() or x2e.lstrip('-').replace('.', '', 1).isdigit() else None
        y2 = float(y2e) if y2e.replace('.', '', 1).isdigit() or y2e.lstrip('-').replace('.', '', 1).isdigit() else None

        # Validación de rangos
        if y > altura_total:
            messagebox.showwarning("Advertencia", f"La posición 'y' ({y} m) no puede superar la altura total ({altura_total} m).")
            return
        if y2 and y2 > altura_total:  # Solo si x2 y y2 no son None
            messagebox.showwarning("Advertencia", f"La posición 'y2' ({y2} m) no puede superar la altura total ({altura_total} m).")
            return

        if abs(x) > ancho_mitad:
            messagebox.showwarning("Advertencia", f"La posición 'x' ({x} m) no puede superar el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
            return
        if x2 and abs(x2) > ancho_mitad:  # Solo si x2 no es None
            messagebox.showwarning("Advertencia", f"La posición 'x2' ({x2} m) no puede superar el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
            return

        # Validar separación solo si hay más de un conductor
        if cantidad > 1:
            separacion = float(entrada_separacion.get())
            if separacion <= 0:
                raise ValueError("La separación debe ser mayor que 0.")
        else:
            separacion = 0

    except ValueError as e:
        messagebox.showwarning("Advertencia", f"Error en los datos ingresados: {e}")
        return

    # Cerrar figuras previas si existen
    plt.close('all')

    # Graficar los puntos
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # Condicional para graficar solo si x2 y y2 no son None
    if x2 is not None and y2 is not None:
        xs = [x, x2]
        ys = [y, y2]
    else:
        xs = [x]
        ys = [y]

    # Llamar a la función para graficar los puntos
    graficar_puntos(fig, ax, xs, ys, cantidad, separacion, ancho_mitad, altura_total)


def generar_poligono(n_lados, radio, centro_x=0, centro_y=0):
    """
    Genera un conjunto de puntos que forman un polígono regular.

    Parámetros:
        n_lados (int): Número de lados (vértices) del polígono.
        radio (float): Distancia desde el centro hasta cada vértice.
        centro_x (float): Coordenada x del centro del polígono.
        centro_y (float): Coordenada y del centro del polígono.

    Retorna:
        lista_puntos (list): Lista de tuplas (x, y) representando las coordenadas de cada vértice.
    """
    theta0 = 2*math.pi/n_lados
    if n_lados < 2:
        raise ValueError("El polígono debe tener al menos 2 lados.")
    if n_lados >= 2:
        theta_1 = (math.pi/(2*n_lados))*(n_lados-2)
        lista_puntos = []
        for i in range(n_lados):
            theta = theta0 * i  # Ángulo correspondiente a cada vértice
            x = centro_x + radio * math.cos(theta-theta_1)
            y = centro_y + radio * math.sin(theta-theta_1)
            lista_puntos.append((x, y))
    return lista_puntos

def graficar_puntos(fig, ax, x, y, cantidad, separacion, ancho_mitad, altura_total):
    # Configurar límites
    ax.set_xlim(-ancho_mitad, ancho_mitad)
    ax.set_ylim(0, altura_total)
    ax.set_title("Distribución de conductores")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(True)

    n = int(cantidad)
    pos = ["izquierdo", "derecho"]
    colors = ["blue", "green"]

    # Dibujar puntos
    for j in range(len(x)):
        # Verificar si los puntos están dentro de los límites
        if abs(x[j]) > ancho_mitad:
            messagebox.showwarning("Advertencia", f"El punto ({x[j]}, {y[j]}) excede el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
        if n == 1:
            ax.plot(x[j], y[j], 'o', color='blue', label="Conductor")
        elif n >= 2:
            radio = (separacion / 2) * (1 / math.sin(math.pi / n))
            puntos = generar_poligono(n, radio, x[j], y[j])
            for i in range(len(puntos)):
                ax.plot(puntos[i][0], puntos[i][1], 'o', color=colors[j], label=f"Conductor {pos[j]}" if i == 0 else None)

    # Agregar leyenda
    ax.legend(loc="upper right")

    # Mostrar gráfico en la interfaz
    global canvas
    if canvas is not None:
        canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(fig, ventana)
    canvas.get_tk_widget().place(x=600, y=20)


def guardar_datos():
    datos = {nombre: entrada.get() for nombre, entrada in entradas.items()}
    messagebox.showinfo("Datos guardados", f"Datos capturados:\n{datos}")

def limpiar_campos():
    for entrada in entradas.values():
        entrada.delete(0, tk.END)
    entrada_separacion.config(state="disabled")  # Desactivar separación al limpiar

def verificar_separacion(*args):
    try:
        cantidad = int(entrada_cantidad.get())
        if cantidad > 1:
            entrada_separacion.config(state="normal")
        else:
            entrada_separacion.delete(0, tk.END)
            entrada_separacion.config(state="disabled")
    except ValueError:  # Si la entrada no es un número
        entrada_separacion.delete(0, tk.END)
        entrada_separacion.config(state="disabled")

def verificar_modo_viento(*args):
    try:
        modo = str(entrada_modo.get())
        if modo == "gradiente":
            entrada_rugosidad_terreno.config(state="normal")
        else:
            entrada_rugosidad_terreno.delete(0, tk.END)
            entrada_rugosidad_terreno.config(state="disabled")
    except ValueError:  # Si la entrada es un número
        entrada_rugosidad_terreno.delete(0, tk.END)
        entrada_rugosidad_terreno.config(state="disabled")

def validar_natural(texto):
    """Valida que el texto ingresado sea un número natural."""
    if texto == "" or (texto.isdigit() and int(texto) > 0):
        return True
    else:
        messagebox.showwarning("Valor inválido", "La casilla 'Cantidad' solo acepta números naturales (1, 2, 3...).")
        return False
    
# Crear la ventana principal
ventana = tk.Tk()
ventana.title("Ajuste de parámetros")
ventana.geometry("1100x700")

# Contenedores principales
secciones = {
    "Características de los conductores": [
        ("Voltaje (kV)", "400"),
        ("Área (mm²)", "100"), 
        ("Posición en x (m)", "0"), 
        ("Posición en y (m)", "5"),
        ("factor conductor", "1"),
        ("Cantidad", ""), 
        ("Separación (m)", ""),
        ("Voltaje 2(kV)", "-400"), # Agregar nuevo campo (deshabilitado inicialmente)
        ("Área 2 (mm²)", "100"),
        ("Posición en x 2 (m)", "0"), 
        ("Posición en y 2 (m)", "5"),
        ("factor conductor 2", "1"),
    ],
    "Características ambientales": [
        ("Temperatura (°C)", "25"), 
        ("Presión (Pa)", ""), 
        ("Humedad (%)", "50"), 
        ("Viento (m/s)", ""), 
        ("Modo (str)", "uniforme"), 
        ("Rugosidad terreno", "0.2")
    ]
}

entradas = {}
for i, (seccion, campos) in enumerate(secciones.items()):
    marco = ttk.LabelFrame(ventana, text=seccion, padding=10)
    
    if seccion == "Características de los conductores":
        marco.grid(row=i, column=0, sticky="w", padx=10, pady=5)  # Columna 0
    #elif seccion == "Características ambientales":
    #    marco.grid(row=i, column=1, sticky="w", padx=10, pady=5)  # Columna 1
    
    if seccion == "Características de los conductores":  # Para la sección de conductores
        # Organizar las casillas en dos columnas
        for j, (nombre, valor_defecto) in enumerate(campos[:7]):  # Las primeras 7 casillas
            ttk.Label(marco, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)
            entrada = ttk.Entry(marco)
            entrada.insert(0, valor_defecto)
            entrada.grid(row=j, column=1, padx=5, pady=2)
            entradas[nombre] = entrada
            # Deshabilitar las casillas extra inicialmente
            if "2" in nombre:
                entrada.config(state="disabled")
        
        # Las casillas "Área 2", "Posición en x 2" y "Posición en y 2" en la segunda columna
        for j, (nombre, valor_defecto) in enumerate(campos[7:], start=0):  # Las tres casillas adicionales
            ttk.Label(marco, text=nombre).grid(row=j, column=2, sticky="w", padx=5, pady=2)
            entrada = ttk.Entry(marco)
            entrada.insert(0, valor_defecto)
            entrada.grid(row=j, column=3, padx=5, pady=2)
            entradas[nombre] = entrada
            # Deshabilitar las casillas extra inicialmente
            if "2" in nombre:
                entrada.config(state="disabled")
    else:
        # Para las otras secciones, organizarlas como estaban
        for j, (nombre, valor_defecto) in enumerate(campos):
            ttk.Label(marco, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)
            entrada = ttk.Entry(marco)
            entrada.insert(0, valor_defecto)
            entrada.grid(row=j, column=1, padx=5, pady=2)
            entradas[nombre] = entrada

# Crear un contenedor para las secciones de "Características ambientales" y "Dimensionamiento"
marco_contenedor = ttk.Frame(ventana)
marco_contenedor.grid(row=1, column=0, columnspan=2, padx=10, pady=5, sticky="w")  # "sticky" con "w" asegura alineación izquierda

# Subcolumnas para acercar "Dimensionamiento" a "Características ambientales"
# Primero, crear la parte de "Características ambientales"
marco_ambientales = ttk.LabelFrame(marco_contenedor, text="Características ambientales", padding=10)
marco_ambientales.grid(row=0, column=0, padx=10, pady=5, sticky="w")

for j, (nombre, valor_defecto) in enumerate(secciones["Características ambientales"]):
    ttk.Label(marco_ambientales, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)
    entrada = ttk.Entry(marco_ambientales)
    entrada.insert(0, valor_defecto)
    entrada.grid(row=j, column=1, padx=5, pady=2)
    entradas[nombre] = entrada

# Crear la parte de "Dimensionamiento" en la subcolumna
marco_dimensionamiento = ttk.LabelFrame(marco_contenedor, text="Dimensionamiento", padding=10)
marco_dimensionamiento.grid(row=0, column=1, padx=10, pady=5, sticky="w")

ttk.Label(marco_dimensionamiento, text="Ancho (m)").grid(row=0, column=0, sticky="w", padx=5, pady=2)
entrada_ancho_total = ttk.Entry(marco_dimensionamiento)
entrada_ancho_total.insert(0, "5")
entrada_ancho_total.grid(row=0, column=1, padx=5, pady=2)

ttk.Label(marco_dimensionamiento, text="Altura (m)").grid(row=1, column=0, sticky="w", padx=5, pady=2)
entrada_altura_total = ttk.Entry(marco_dimensionamiento)
entrada_altura_total.insert(0, "10")
entrada_altura_total.grid(row=1, column=1, padx=5, pady=2)
# Variables específicas
entrada_cantidad = entradas["Cantidad"]
entrada_separacion = entradas["Separación (m)"]
entrada_posicion_x = entradas["Posición en x (m)"]
entrada_posicion_y = entradas["Posición en y (m)"]

entrada_posicion_x2 =  entradas["Posición en x 2 (m)"]
entrada_posicion_y2 =  entradas["Posición en y 2 (m)"]


# Validación de números naturales
validacion_natural = ventana.register(validar_natural)
entrada_cantidad.config(validate="key", validatecommand=(validacion_natural, "%P"))

entrada_cantidad.bind("<KeyRelease>", verificar_separacion)  # Asociar el evento de escritura
entrada_separacion.config(state="disabled")  # Desactivada inicialmente

# Validación coeficiente rugosidad piso para  viento
entrada_modo = entradas["Modo (str)"]
entrada_rugosidad_terreno = entradas["Rugosidad terreno"]
entrada_modo.bind("<KeyRelease>", verificar_modo_viento)  # Asociar el evento de escritura
entrada_rugosidad_terreno.config(state='disabled') # Desactivada inicialmente

# Botones principales
marco_botones = ttk.Frame(ventana, padding=10)
marco_botones.grid(row=len(secciones) + 1, column=0, columnspan=2, sticky="w", padx=10, pady=5)  # Usar grid

ttk.Button(marco_botones, text="Guardar", command=guardar_datos).grid(row=0, column=0, padx=5)
ttk.Button(marco_botones, text="Limpiar", command=limpiar_campos).grid(row=0, column=1, padx=5)
ttk.Button(marco_botones, text="Graficar", command=validar_campos_y_graficar).grid(row=0, column=2, padx=5)
ttk.Button(marco_botones, text="Ejecutar", command=validar_campos_y_ejecutar).grid(row=0, column=3, padx=5)

# Crear el botón Bipolar
boton_bipolar_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_bipolar = ttk.Checkbutton(ventana, text="Bipolar", variable=boton_bipolar_var, command=activar_bipolar)
boton_bipolar.grid(row=len(secciones) + 2, column=0, padx=5, pady=10)

# Botones de gráficos
marco_graficos = ttk.LabelFrame(ventana, text="Gráficos", padding=10)
marco_graficos.grid(row=len(secciones) + 3, column=0, columnspan=2, padx=10, pady=5)

# Distribuir botones en 2 filas y 4 columnas
num_columnas = 4  # Número de columnas deseadas
for i, grafico in enumerate(estados_graficos):
    fila = i // num_columnas  # División entera para obtener la fila
    columna = i % num_columnas  # Resto de la división para obtener la columna
    ttk.Button(marco_graficos, text=grafico, command=lambda g=grafico: alternar_grafico(g)).grid(row=fila, column=columna, padx=5, pady=5)

# Espacio para el gráfico
canvas = None

# Iniciar la aplicación
ventana.mainloop()










