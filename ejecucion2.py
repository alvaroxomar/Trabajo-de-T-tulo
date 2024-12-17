import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Estado de gráficos
estados_graficos = {"Gráfico 1": False, "Gráfico 2": False}

def alternar_grafico(nombre):
    estados_graficos[nombre] = not estados_graficos[nombre]
    mensaje = "seleccionado" if estados_graficos[nombre] else "deseleccionado"
    messagebox.showinfo(nombre, f"Se ha {mensaje} {nombre}.")

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
        if not entrada.get():  # Campo vacío
            messagebox.showwarning("Advertencia", f"No pueden haber casillas vacías (revisar: {nombre}).")
            return

        try:
            float(entrada.get())
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

        if y > altura_total:
            messagebox.showwarning("Advertencia", f"La posición 'y' ({y} m) no puede superar la altura total ({altura_total} m).")
            return

        if abs(x) > ancho_mitad:
            messagebox.showwarning("Advertencia", f"La posición 'x' ({x} m) no puede superar el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
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
    graficar_puntos(fig, ax, x, y, cantidad, separacion, ancho_mitad, altura_total)

def graficar_puntos(fig, ax, x, y, cantidad, separacion, ancho_mitad, altura_total):
    # Configurar límites
    ax.set_xlim(-ancho_mitad, ancho_mitad)
    ax.set_ylim(0, altura_total)
    ax.set_title("Distribución de conductores")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(True)

    # Dibujar puntos
    for i in range(cantidad):
        punto_x = x + i * separacion
        if abs(punto_x) > ancho_mitad:
            messagebox.showwarning("Advertencia", f"El punto ({punto_x}, {y}) excede el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
            continue
        ax.plot(punto_x, y, 'o', color='blue', label="Conductor" if i == 0 else None)

    # Agregar leyenda
    ax.legend(loc="upper right")

    # Mostrar gráfico en la interfaz
    global canvas  # Para evitar que el gráfico se borre al salir de la función
    if canvas is not None:
        canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(fig, ventana)
    canvas.get_tk_widget().place(x=460, y=20)

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
ventana.geometry("800x550")

# Contenedores principales
secciones = {
    "Características de los conductores": [("Área (mm²)", "100"), ("Cantidad", ""), ("Separación (m)", ""),
                                           ("Posición en x (m)", "0"), ("Posición en y (m)", "0")],
    "Características ambientales": [("Temperatura (°C)", "25"), ("Presión (Pa)", ""), ("Humedad (%)", "50"), ("Viento (m/s)", "")]
}

entradas = {}
for i, (seccion, campos) in enumerate(secciones.items()):
    marco = ttk.LabelFrame(ventana, text=seccion, padding=10)
    marco.pack(fill="x", padx=10, pady=5)
    for j, (nombre, valor_defecto) in enumerate(campos):
        ttk.Label(marco, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)
        entrada = ttk.Entry(marco)
        entrada.insert(0, valor_defecto)
        entrada.grid(row=j, column=1, padx=5, pady=2)
        entradas[nombre] = entrada

# Agregar categoría de dimensionamiento
marco_dimensionamiento = ttk.LabelFrame(ventana, text="Dimensionamiento", padding=10)
marco_dimensionamiento.pack(fill="x", padx=10, pady=5)

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

# Validación de números naturales
validacion_natural = ventana.register(validar_natural)
entrada_cantidad.config(validate="key", validatecommand=(validacion_natural, "%P"))

entrada_cantidad.bind("<KeyRelease>", verificar_separacion)  # Asociar el evento de escritura
entrada_separacion.config(state="disabled")  # Desactivada inicialmente

# Botones principales
marco_botones = ttk.Frame(ventana, padding=10)
marco_botones.pack(fill="x", padx=10, pady=5)

ttk.Button(marco_botones, text="Guardar", command=guardar_datos).grid(row=0, column=0, padx=5)
ttk.Button(marco_botones, text="Limpiar", command=limpiar_campos).grid(row=0, column=1, padx=5)
ttk.Button(marco_botones, text="Graficar", command=validar_campos_y_graficar).grid(row=0, column=2, padx=5)
ttk.Button(marco_botones, text="Ejecutar", command=validar_campos_y_ejecutar).grid(row=0, column=3, padx=5)

# Botones de gráficos
marco_graficos = ttk.LabelFrame(ventana, text="Gráficos", padding=10)
marco_graficos.pack(fill="x", padx=10, pady=5)

for i, grafico in enumerate(estados_graficos):
    ttk.Button(marco_graficos, text=grafico, command=lambda g=grafico: alternar_grafico(g)).grid(row=0, column=i, padx=5)

# Espacio para el gráfico
canvas = None

# Iniciar la aplicación
ventana.mainloop()









