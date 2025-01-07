import tkinter as tk
from tkinter import ttk, messagebox, font
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
import numpy as np
import json
import subprocess

plt.close('all')

# Estado de gráficos
#'Eele', 'Vele', 'Rhof', 'Vf', 'Vdef', 'Ef', 'J1', 'E1', 'SubPl'
estados_graficos = {
    "Eele": False, "Vele": False, "Rhof": False, "Vf": False, "SPd": False, "Rhop": False,
    "Vdef": False, "Ef": False, "J1": False, "E1": False, "SPv": False, "Rhon": False
}
graficos = []

def alternar_grafico(nombre):
    """
    Alterna el estado de un gráfico, actualiza la lista de gráficos seleccionados y muestra un mensaje.
    
    Args:
        nombre (str): Nombre del gráfico.
    """
    # Cambiar el estado del gráfico
    estados_graficos[nombre] = not estados_graficos[nombre]
    mensaje = "seleccionado" if estados_graficos[nombre] else "deseleccionado"
    # Actualizar la lista de gráficos activados
    if estados_graficos[nombre]:
        if nombre not in graficos:
            graficos.append(nombre)  # Agregar a la lista si está seleccionado
    else:
        if nombre in graficos:
            graficos.remove(nombre)  # Eliminar de la lista si está deseleccionado
    
    # Mostrar el mensaje
    messagebox.showinfo(nombre, f"Se ha {mensaje} {nombre}.")
    
    # Imprimir el estado actual de la lista (opcional)
    print(f"Gráficos seleccionados: {graficos}")



# Botón Bipolar para activar/desactivar las casillas extra
def activar_bipolar():
    # Activar o desactivar las casillas extras
    if boton_bipolar_var.get():  # Si se selecciona el botón Bipolar
        entradas["Voltaje 2 (kV)"].config(state="normal")
        entradas["Área 2 (mm²)"].config(state="normal")
        entradas["Posición en x 2 (m)"].config(state="normal")
        entradas["Posición en y 2 (m)"].config(state="normal")
        entradas["factor conductor 2"].config(state="normal")
        entradas["Movilidad iónica 2 (m2/kVs)"].config(state="normal")
        entradas["Recombinación (μm^3/s)"].config(state="normal")
        entradas["Jp (nA/m^2)"].config(state="disabled")
        # Verifica si alguna casilla está activada
        boton_Spd.config(state=tk.NORMAL)
        boton_Rhop.config(state=tk.NORMAL)
        boton_Rhon.config(state=tk.NORMAL)
    else:  # Si el botón Bipolar está desmarcado
        # Eliminar el contenido actual de las casillas
        entradas["Voltaje 2 (kV)"].config(state="disabled")
        entradas["Área 2 (mm²)"].config(state="disabled")
        entradas["Posición en x 2 (m)"].config(state="disabled")
        entradas["Posición en y 2 (m)"].config(state="disabled")
        entradas["factor conductor 2"].config(state="disabled")
        entradas["Movilidad iónica 2 (m2/kVs)"].config(state="disabled")
        entradas["Recombinación (μm^3/s)"].config(state="disabled")
        entradas["Jp (nA/m^2)"].config(state="normal")
        boton_Spd.config(state=tk.DISABLED)
        boton_Rhop.config(state=tk.DISABLED)
        boton_Rhon.config(state=tk.DISABLED)

        
        # Establecer x2 y y2 como None para que no se grafiquen
        #global x2, y2  # Si las variables x2 y y2 son globales
        #x2, y2 = None, None
        
        # Si ya había un gráfico previo, eliminar los puntos asociados a x2, y2
        plt.close('all')  # Cierra el gráfico previo
        # Redibujar solo los puntos x, y (sin x2, y2)
        validar_campos_y_graficar()  # Esto vuelve a graficar los puntos con las nuevas configuraciones

def activar_auto_red():
    """Controla el estado de las casillas según el botón 'Auto red'."""
    if boton_auto_red_var.get():
        # Habilitar todas las casillas temporalmente para permitir la selección
        for entrada in ["Ancho (m)", "Altura (m)", "nodos x", "nodos y"]:
            entradas[entrada].config(state="normal")

        # Vincular eventos para detectar el primer grupo que se edite
        entradas["Ancho (m)"].bind("<FocusIn>", lambda e: habilitar_grupo("dimensiones"))
        entradas["Altura (m)"].bind("<FocusIn>", lambda e: habilitar_grupo("dimensiones"))
        entradas["nodos x"].bind("<FocusIn>", lambda e: habilitar_grupo("nodos"))
        entradas["nodos y"].bind("<FocusIn>", lambda e: habilitar_grupo("nodos"))
    else:
        for entrada in ["Ancho (m)", "Altura (m)", "nodos x", "nodos y"]:
            entradas[entrada].config(state="normal")
        # Restablecer los valores y habilitar todas las casillas
        restablecer_valores()
        # Desvincular eventos
        for entrada in ["Ancho (m)", "Altura (m)", "nodos x", "nodos y"]:
            entradas[entrada].unbind("<FocusIn>")

def habilitar_grupo(grupo):
    """Habilitar un grupo y deshabilitar el opuesto después de la interacción del usuario."""
    if grupo == "dimensiones":
        # Habilitar dimensiones, deshabilitar nodos
        entradas["Ancho (m)"].config(state="normal")
        entradas["Altura (m)"].config(state="normal")
        # Limpiar y deshabilitar nodos
        entradas["nodos x"].delete(0, tk.END)
        entradas["nodos y"].delete(0, tk.END)
        entradas["nodos x"].config(state="disabled")
        entradas["nodos y"].config(state="disabled")
    elif grupo == "nodos":
        # Habilitar nodos, deshabilitar dimensiones
        entradas["nodos x"].config(state="normal")
        entradas["nodos y"].config(state="normal")
        # Limpiar y deshabilitar dimensiones
        entradas["Ancho (m)"].delete(0, tk.END)
        entradas["Altura (m)"].delete(0, tk.END)
        entradas["Ancho (m)"].config(state="disabled")
        entradas["Altura (m)"].config(state="disabled")

    # Desvincular eventos de FocusIn una vez que se seleccionó un grupo
    for entrada in ["Ancho (m)", "Altura (m)", "nodos x", "nodos y"]:
        entradas[entrada].unbind("<FocusIn>")

def restablecer_valores():
    """Rellena las casillas inhabilitadas con valores por defecto y conserva los valores existentes en las habilitadas."""
    valores_por_defecto = {"Ancho (m)": "5","Altura (m)": "10","nodos x": "100","nodos y": "100",}

    for entrada in ["Ancho (m)", "Altura (m)", "nodos x", "nodos y"]:
        # Si la casilla está vacía, colocar el valor por defecto
        if not entradas[entrada].get():
            entradas[entrada].insert(0, valores_por_defecto[entrada])


def ejecutar_script():
    # Guardar todos los parametros en un diccionario
    params = {nombre: entrada.get() for nombre, entrada in entradas.items()}
    params['graficos'] = graficos
    params['guardar'] = estado_guardado
    # Crear un mensaje con los valores
    mensaje = "\n".join([f"{clave}: {valor}" for clave, valor in params.items()])
    # Mostrar un cuadro de diálogo con los valores
    messagebox.showinfo("Parámetros Guardados", mensaje)
    # Convertir los parámetros a JSON
    params_json = json.dumps(params)
    if boton_bipolar_var.get():
        # Ejecutar el script bipolar con los parámetros
        subprocess.run(["python", "script_ejecucion_bip.py", params_json], check=True)
        messagebox.showinfo("Éxito", "El archivo bipolar.py se ejecutó correctamente.")
    else:
        # Ejecutar el script unipolar con los parámetros
        subprocess.run(["python", "script_ejecucion_uni.py", params_json], check=True)
        messagebox.showinfo("Éxito", "El archivo unipolar.py se ejecutó correctamente.")
    

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
        if nombre == "Ancho (m)":
            continue  # Ignorar esta casilla
        if nombre == "Altura (m)":
            continue
        if nombre == "nodos x":
            continue  # Ignorar esta casilla
        if nombre == "nodos y":
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
        # Validar "Modo (str)" para aceptar solo "uniforme" o "gradiente"
        if nombre == "Interior conductor":
            if valor.isdigit():
                messagebox.showwarning(
                    "Advertencia", f"El valor en '{nombre}' no puede ser numérico. Debe ser 'si', 'no' o 'noV'."
                )
                return
            if valor not in ["si", "no", "noV"]:
                messagebox.showwarning(
                    "Advertencia", f"El valor ingresado en '{nombre}' debe ser 'si', 'no' o 'noV'."
                )
                return
            continue  # Saltar el resto de la validación para este campo

        # Validar que los demás campos sean números válidos
        try:
            float(valor)
        except ValueError:
            messagebox.showwarning("Advertencia", f"El valor ingresado en '{nombre}' debe ser un número.")
            return

    campos_dimensionamiento = ["Ancho (m)", "Altura (m)", "nodos x", "nodos y", "Medición (m)"]

    # Verificar si los valores no vacíos son números válidos
    for campo in campos_dimensionamiento:
        valor = entradas[campo].get()
        if valor != "":  # Permitir valores vacíos
            try:
                # Intentar convertir el valor a float
                float(valor)
            except ValueError:
                messagebox.showwarning("Advertencia", f"El valor de {campo} debe ser un número válido o estar vacío.")
                return

    messagebox.showinfo("Validación exitosa", "Todos los campos han sido completados correctamente.")
    ejecutar_script()

# Inicialización global de canvas
canvas = None
# Llamar al método al iniciar la ventana para mostrar un gráfico por defecto
def inicializar_grafico():
    # Valores por defecto
    x = float(entrada_posicion_x.get())
    y = float(entrada_posicion_y.get())

    # Validar datos de dimensionamiento
    ancho_mitad = float(entradas["Ancho (m)"].get())
    altura_total = float(entradas["Altura (m)"].get())
    cantidad = float(entradas["Subconductores"].get())
    separacion = entrada_separacion


    # Crear una figura inicial
    fig, ax = plt.subplots(figsize=(5, 5))
    xs = [x]
    ys = [y]

    # Graficar puntos iniciales
    graficar_puntos(fig, ax, xs, ys, cantidad, separacion, ancho_mitad, altura_total)


def validar_campos_y_graficar():
    # Cerrar todas las figuras de Matplotlib previas si existen
    plt.close('all')  # Asegura que las figuras de Matplotlib se cierren correctamente
    try:
        x = float(entrada_posicion_x.get())
        y = float(entrada_posicion_y.get())

        # Validar datos de dimensionamiento
        ancho_mitad = entradas["Ancho (m)"].get()
        altura_total = entradas["Altura (m)"].get()

        # Si están vacíos, asignar valores de entrada_posicion_x y entrada_posicion_y
        if not ancho_mitad:  # Si ancho_mitad está vacío
            ancho_mitad = y * 2  # Asignar valor de nodos x
        else:
            ancho_mitad = float(ancho_mitad)  # Convertir el valor a float si no está vacío

        if not altura_total:  # Si altura_total está vacío
            altura_total = y * 2  # Asignar valor de nodos y
        else:
            altura_total = float(altura_total)  # Convertir el valor a float si no está vacío

        if ancho_mitad <= 0 or altura_total <= 0:
            raise ValueError("El ancho y la altura deben ser mayores que 0.")

        # Validar datos de conductores
        cantidad = int(entrada_cantidad.get())
        if cantidad <= 0:
            raise ValueError("La cantidad debe ser mayor que 0.")
        
        # Obtener las posiciones x2 y y2
        x2e = entrada_posicion_x2.get()
        y2e = entrada_posicion_y2.get()
        x2 = float(x2e) if x2e.replace('.', '', 1).isdigit() or x2e.lstrip('-').replace('.', '', 1).isdigit() else None
        y2 = float(y2e) if y2e.replace('.', '', 1).isdigit() or y2e.lstrip('-').replace('.', '', 1).isdigit() else None

        # Validación de rangos
        if y > altura_total:
            messagebox.showwarning("Advertencia", f"La posición 'y' ({y} m) no puede superar la altura total ({altura_total} m).")
            return
        if y2 and abs(y2) > altura_total:  # Solo si x2 y y2 no son None
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

    # Graficar los puntos
    fig, ax = plt.subplots(figsize=(5, 5))

    # Condicional para graficar solo si x2 y y2 no son None
    if x2 is not None and y2 is not None and boton_bipolar_var.get():
        xs = [x, x2]
        ys = [y, y2]
    else:
        xs = [x]
        ys = [y]

    # Llamar a la función para graficar los puntos
    graficar_puntos(fig, ax, xs, ys, cantidad, separacion, ancho_mitad, altura_total)

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
        elif n == 1:
            if len(x) == 1:
                ax.plot(x[j], y[j], 'o', color=colors[j], label=f"Conductor")
            else:
                ax.plot(x[j], y[j], 'o', color=colors[j], label=f"Conductor {pos[j]}")
        elif n >= 2:
            radio = ((separacion / 100) / 2) * (1 / math.sin(math.pi / n))
            puntos = generar_poligono(n, radio, x[j], y[j])
            for i in range(len(puntos)):
                ax.plot(puntos[i][0], puntos[i][1], 'o', color=colors[j], label=f"Conductor {pos[j]}" if i == 0 else None)

    # Agregar leyenda
    ax.legend(loc="upper right")

    # Mostrar gráfico en la interfaz
    global canvas
    if canvas is not None:
        canvas.get_tk_widget().destroy()  # Destruir el canvas anterior

    # Crear y mostrar el nuevo canvas
    canvas = FigureCanvasTkAgg(fig, ventana)  # Crear nuevo canvas
    canvas.get_tk_widget().place(x=650, y=20)  # Mostrar el nuevo gráfico



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

estado_guardado = False

def guardar_datos():
    global estado_guardado  # Referencia la variable global
    
    # Alternar entre True y False
    if estado_guardado:
        estado_guardado = False
        messagebox.showinfo("Datos no guardados", "No se guardaran los datos una vez ejecutado el programa")
    else:
        datos = {nombre: entrada.get() for nombre, entrada in entradas.items() if isinstance(entrada, ttk.Entry)}
        messagebox.showinfo("Datos guardados", f"Datos capturados:\n{datos}")
        estado_guardado = True




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
            entrada_separacion.insert(0, "0")
            entrada_separacion.config(state="disabled")
    except ValueError:  # Si la entrada no es un número
        entrada_separacion.delete(0, tk.END)
        entrada_separacion.insert(0, "0")
        entrada_separacion.config(state="disabled")

def verificar_modo_viento(*args):
    try:
        modo = str(entrada_modo.get())
        if modo == "gradiente":
            entrada_rugosidad_terreno.config(state="normal")
        else:
            entrada_rugosidad_terreno.delete(0, tk.END)
            entrada_rugosidad_terreno.insert(0, "0.2")
            entrada_rugosidad_terreno.config(state="disabled")
    except ValueError:  # Si la entrada es un número
        entrada_rugosidad_terreno.delete(0, tk.END)
        entrada_rugosidad_terreno.insert(0, "0.2")
        entrada_rugosidad_terreno.config(state="disabled")


def validar_natural(texto, nombre):
    """Valida que el texto ingresado sea un número natural."""
    if texto == "" or (texto.isdigit() and int(texto) > 0):
        return True
    else:
        messagebox.showwarning("Valor inválido", f"La casilla '{nombre}' solo acepta números naturales (1, 2, 3...).")
        return False
def validar_uno(texto, nombre):
    """Valida que el texto ingresado sea un número real entre 0 y 1."""
    if texto == "":  # Permitir casillas vacías mientras el usuario escribe
        return True
    try:
        valor = float(texto)  # Intentar convertir el texto a flotante
        if 0 <= valor <= 1:  # Verificar si está en el rango permitido
            return True
        else:
            raise ValueError  # Forzar excepción si el valor no está en rango
    except ValueError:
        # Mostrar advertencia solo una vez
        messagebox.showwarning(
            "Valor inválido",
            f"La casilla '{nombre}' solo acepta números reales entre 0 y 1."
        )
        return False
    
def cerrar_programa():
    # Confirmación para evitar cierres accidentales
    if messagebox.askokcancel("Salir", "¿Estás seguro de que quieres salir?"):
        # Liberar recursos, cerrar gráficos
        plt.close('all')
        # Terminar el programa
        if ventana.winfo_exists():  # Verifica si la ventana aún existe
            ventana.destroy()
        plt.close('all')


# Crear la ventana principal
ventana = tk.Tk()
ventana.title("Ajuste de parámetros")
ventana.geometry("1200x700")

# Diccionario con explicaciones de los parámetros
explicaciones = {
    "Voltaje (kV)": "Voltaje aplicado al conductor, medido en kilovoltios (kV).",
    "Área (mm²)": "Sección transversal del conductor, medida en milímetros cuadrados.",
    "Posición en x (m)": "Posición horizontal del conductor 1 en metros.",
    "Posición en y (m)": "Posición vertical del conductor 1 en metros.",
    "Posición en x 2 (m)": "Posición horizontal del conductor 2 en metros.",
    "Posición en y 2 (m)": "Posición vertical del conductor 2 en metros.",
    "factor conductor": "Factor de rugosidad para las características del conductor.",
    "factor conductor 2": "Factor de rugosidad para las características del conductor 2.",
    "Subconductores": "Número de subconductores usados en cada polo.",
    "Separación (cm)": "Distancia entre subconductores, medida en centímetros.",
    "Voltaje 2 (kV)": "Voltaje aplicado al conductor secundario en kilovoltios (kV).",
    "Área 2 (mm²)": "Sección transversal del segundo conductor, medida en milímetros cuadrados.",
    "Recombinación (μm^3/s)": "Tasa de recombinación del material, medida en micrómetros cúbicos por segundo.",
    "Jp (nA/m^2)": "Densidad de corriente de polarización, medida en nanoamperios por metro cuadrado.",
    "Movilidad iónica (m2/kVs)": "Movilidad de iones en el sistema.",
    "Temperatura (°C)": "Temperatura ambiente en grados Celsius.",
    "Altitud (m)": "Altitud por sobre el nivel del mar (m)",
    "Viento x (m/s)": "Componente horizontal de la velocidad del viento.",
    "Viento y (m/s)": "Componente vertical de la velocidad del viento.",
    "Modo (str)": (
        "Modo de distribución del viento:\n"
        "- uniforme: significa que a toda altura existe la misma velocidad horizontal de viento\n"
        "- gradiente: significa que la velocidad del viento se distribuye de manera gradual aumentando proporcionalmente "
        "con la altura de acuerdo con la fórmula wx(y)=wx0(y/y0)^alpha, con wx0 una velocidad de referencia a una altura y0 "
        "y alpha es un coeficiente de rugosidad del suelo (0,1)."
    ),
    "Rugosidad terreno": "Coeficiente que describe la rugosidad del terreno, su valor comprende entre 0 y 1.",
    "Interior conductor": (
        "Indicación sobre la forma en que se impone la condición de borde en el conductor:\n"
        "- si: significa que la densidad de carga y el voltaje se distribuyen en un único nodo ubicado en el centro del conductor equivalente\n"
        "- no: significa que la densidad de carga se distribuye en la periferia del conductor de manera no uniforme, "
        "con el máximo en el nodo sur, el mínimo en el nodo norte e intermedio en los nodos este y oeste\n"
        "- noV: es exclusivo para la distribución del voltaje en el cual este se distribuye también en la periferia del conductor equivalente."
    ),
    "Gráficos disponibles": (
        "Gráficos disponibles:\n"
        "- Eele: Campo electrostático\n"
        "- Vele: Potencial electrostático\n"
        "- Rhof: Densidad de carga espacial total\n"
        "- Vf: Potencial iónico\n"
        "- SPd: Densidades de carga vista 3D, disponible solamente para la configuración bipolar\n"
        "- Rhop: Densidad de carga positiva, disponible solamente para configuración bipolar\n"
        "- Rhon: Densidad de carga negativa, disponible solamente para configuración bipolar\n"
        "- Vdef: Potencial definitivo\n"
        "- Ef: Campo definitivo\n"
        "- J1: Perfil densidad de corriente a nivel de piso\n"
        "- E1: Perfil magnitud campo eléctrico total a nivel de piso\n"
        "- SPv: Potenciales vista 3D\n"
    ),
    "Ancho (m)": "Semi ancho del entorno en el que se encuentran los conductores.",
    "Altura (m)": "Altura total del entorno en el que se encuentran los conductores.",
    "nodos x": "Cantidad total de nodos en dirección horizontal.",
    "nodos y": "Cantidad total de nodos en dirección vertical.",
    "Medición (m)": "Altura en metros en donde se realizará el cálculo del campo eléctrico y la densidad de corriente iónica.",
    "Max iter rho": "Máxima cantidad de iteraciones para calcular la densidad de carga en el espacio.",
    "Max iter V": "Máxima cantidad de iteraciones para calcular el potencial en el espacio.",
    "Max iter Gob": "Máxima cantidad de iteraciones para el algoritmo principal.",
}

# Configurar una fuente más pequeña para el botón
fuente_boton = ("Arial", 6, "bold")  # Tamaño 6 para reducir altura
# Función para mostrar una ventana emergente con la explicación
def mostrar_explicacion(parametro):
    if parametro in explicaciones:
        # Crear una ventana emergente
        ventana_emergente = tk.Toplevel(ventana)
        ventana_emergente.title(f"Explicación: {parametro}")
        ventana_emergente.geometry("500x400")  # Tamaño fijo de la ventana

        # Crear un marco para contener el desplazamiento
        marco = ttk.Frame(ventana_emergente)
        marco.pack(fill="both", expand=True)

        # Crear un Canvas para el contenido desplazable
        canvas = tk.Canvas(marco)
        canvas.pack(side="left", fill="both", expand=True)

        # Barra de desplazamiento vertical
        barra_desplazamiento = ttk.Scrollbar(
            marco, orient="vertical", command=canvas.yview
        )
        barra_desplazamiento.pack(side="right", fill="y")

        # Conectar el Canvas con la barra de desplazamiento
        canvas.configure(yscrollcommand=barra_desplazamiento.set)

        # Crear un marco dentro del Canvas para el contenido
        marco_contenido = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=marco_contenido, anchor="nw")

        # Agregar el texto de explicación al marco interno
        texto = tk.Label(
            marco_contenido,
            text=explicaciones[parametro],
            wraplength=450,  # Ajusta el ancho del texto
            justify="left",
            padx=10,
            pady=10
        )
        texto.pack(fill="both", expand=True)

        # Ajustar el área visible del Canvas al tamaño del contenido
        marco_contenido.update_idletasks()
        canvas.config(scrollregion=canvas.bbox("all"))

        # Crear un marco para el botón de cierre
        marco_boton = ttk.Frame(ventana_emergente)
        marco_boton.pack(fill="x", pady=10)

        # Botón para cerrar la ventana emergente
        boton_cerrar = ttk.Button(
            marco_boton,
            text="Cerrar",
            command=ventana_emergente.destroy
        )
        boton_cerrar.pack(pady=20)  # Centrado automáticamente al usar .pack()


       

# Contenedores principales
secciones = {
    "Características de los conductores": [
        ("Voltaje (kV)", "400"),
        ("Área (mm²)", "100"),
        ("Posición en x (m)", "0"),
        ("Posición en y (m)", "5"),
        ("factor conductor", "1"),
        ("Subconductores", "1"),
        ("Separación (cm)", "0"),
        ("Voltaje 2 (kV)", "-400"),
        ("Área 2 (mm²)", "100"),
        ("Posición en x 2 (m)", "0"),
        ("Posición en y 2 (m)", "5"),
        ("factor conductor 2", "1"),
        ("Recombinación (μm^3/s)", "1.8"),
        ("Jp (nA/m^2)", "1988e-8"),
    ],
    "Características ambientales": [
        ("Movilidad iónica (m2/kVs)", "0.15"),
        ("Movilidad iónica 2 (m2/kVs)", "0.15"),
        ("Temperatura (°C)", "25"),
        ("Altitud (m)", "100"),
        ("Viento x (m/s)", "0"),
        ("Viento y (m/s)", "0"),
        ("Modo (str)", "uniforme"),
        ("Rugosidad terreno", "0.2"),
        ("Interior conductor", "si")
    ]
}

entradas = {}
opciones_areas = ["50", "70", "95", "120", "150", "185", "240", "300", "400", "500", "630"]  # Áreas en mm²
opciones_modo = ["uniforme",  "gradiente"]
opciones_condct = ["si", "no", "noV"]
for i, (seccion, campos) in enumerate(secciones.items()):
    if seccion == "Características de los conductores":
        # Marco contenedor para incluir título y marco decorado
        marco_contenedor = tk.Frame(ventana)
        marco_contenedor.grid(row=i, column=0, sticky="w", padx=10, pady=5)

        # Título destacado
        titulo = ttk.Label(
            marco_contenedor,
            text="Características de los conductores",
            background="lightblue",
            font=("Arial", 10, "bold")
        )
        titulo.grid(row=0, column=0, sticky="w", padx=5, pady=(0, 5))

        # Marco decorado
        marco = tk.Frame(
            marco_contenedor,
            highlightthickness=3,
            highlightbackground="gray",
            highlightcolor="green",
            bg="lightgray"
        )
        marco.grid(row=1, column=0, sticky="w")

        # Organizar las casillas en dos columnas
        for j, (nombre, valor_defecto) in enumerate(campos[:7]):
            ttk.Label(marco, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)

            # Usar Combobox para los campos de área
            if "Área" in nombre:
                entrada = ttk.Combobox(marco, values=opciones_areas, state="readonly")
                entrada.set(valor_defecto)  # Valor por defecto
            else:
                entrada = ttk.Entry(marco)
                entrada.insert(0, valor_defecto)

            entrada.grid(row=j, column=1, padx=5, pady=2)
            entradas[nombre] = entrada
            # Botón de ayuda
            boton_ayuda = ttk.Button(marco, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
            boton_ayuda.grid(row=j, column=2, padx=2, pady=1)
            # Aplicar fuente personalizada al botón
            boton_ayuda.configure(style="Ayuda.TButton")
            # Crear un estilo para el botón
            style = ttk.Style()
            style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido
            # Deshabilitar las casillas extra inicialmente
            if "2" in nombre and nombre != "Jp (nA/m^2)" or nombre == "Recombinación (μm^3/s)":
                entrada.config(state="disabled")

        # Las casillas adicionales en la segunda columna
        for j, (nombre, valor_defecto) in enumerate(campos[7:], start=0):
            ttk.Label(marco, text=nombre).grid(row=j, column=3, sticky="w", padx=5, pady=2)

            # Usar Combobox para los campos de área
            if "Área" in nombre:
                entrada = ttk.Combobox(marco, values=opciones_areas, state="readonly")
                entrada.set(valor_defecto)  # Valor por defecto
            else:
                entrada = ttk.Entry(marco)
                entrada.insert(0, valor_defecto)

            entrada.grid(row=j, column=4, padx=5, pady=2)
            entradas[nombre] = entrada
            # Botón de ayuda
            boton_ayuda = ttk.Button(marco, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
            boton_ayuda.grid(row=j, column=5, padx=2, pady=1)
            # Aplicar fuente personalizada al botón
            boton_ayuda.configure(style="Ayuda.TButton")
            # Crear un estilo para el botón
            style = ttk.Style()
            style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido
            # Deshabilitar las casillas extra inicialmente
            if "2" in nombre and nombre != "Jp (nA/m^2)" or nombre == "Recombinación (μm^3/s)":
                entrada.config(state="disabled")
    else:
        # Para las otras secciones
        marco = ttk.LabelFrame(ventana, text=seccion, padding=10)
        marco.grid(row=i, column=0, sticky="w", padx=10, pady=5)

        for j, (nombre, valor_defecto) in enumerate(campos):
            ttk.Label(marco, text=nombre).grid(row=j, column=0, sticky="w", padx=5, pady=2)
            # Usar Combobox para los campos de modo
            if "Modo" in nombre:
                entrada = ttk.Combobox(marco, values=opciones_modo, state="readonly")
                entrada.set(valor_defecto)  # Valor por defecto
            else:
                entrada = ttk.Entry(marco)
                entrada.insert(0, valor_defecto)
            entrada.grid(row=j, column=1, padx=5, pady=2)
            entradas[nombre] = entrada
            

# Crear un contenedor para las secciones
marco_contenedor = ttk.Frame(ventana)
marco_contenedor.grid(row=1, column=0, columnspan=2, padx=10, pady=5, sticky="w")  # Contenedor general

# Parte de "Características ambientales"
marco_ambientales = tk.Frame(
    marco_contenedor,
    highlightthickness=3,            # Grosor del borde
    highlightbackground="gray",     # Color del borde sin enfoque
    highlightcolor="cyan",          # Color del borde con enfoque
    bg="lightgray"                  # Fondo del marco
)
marco_ambientales.grid(row=0, column=0, padx=10, pady=5, sticky="w")
# Definir la fuente personalizada
fuente_titulo = font.Font(family="Arial", size=9, weight="bold")
ttk.Label(marco_ambientales, text="Características ambientales", background="lightblue",font=fuente_titulo).grid(row=0, column=0, columnspan=2, sticky="w")

for j, (nombre, valor_defecto) in enumerate(secciones["Características ambientales"]):
    ttk.Label(marco_ambientales, text=nombre).grid(row=j + 1, column=0, sticky="w", padx=5, pady=2)
    # Usar Combobox para los campos de modo
    if "Modo" in nombre:
        entrada = ttk.Combobox(marco_ambientales, values=opciones_modo, state="readonly")
        entrada.set(valor_defecto)  # Valor por defecto
    elif "Interior" in nombre:
        entrada = ttk.Combobox(marco_ambientales, values=opciones_condct, state="readonly")
        entrada.set(valor_defecto)  # Valor por defecto
    else:
        entrada = ttk.Entry(marco_ambientales)
        entrada.insert(0, valor_defecto)
    entrada.grid(row=j + 1, column=1, padx=5, pady=2)
    entradas[nombre] = entrada
    boton_ayuda = ttk.Button(marco_ambientales, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
    boton_ayuda.grid(row=j+1, column=2, padx=2, pady=1)
    # Aplicar fuente personalizada al botón
    boton_ayuda.configure(style="Ayuda.TButton")
    # Crear un estilo para el botón
    style = ttk.Style()
    style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido
    # Deshabilitar las casillas extra inicialmente
    if nombre == "Movilidad iónica 2 (m2/kVs)":
        entradas[nombre].config(state="disabled")

# Contenedor para "Dimensionamiento" e "Iteraciones" (subcolumna)
marco_segunda_columna = ttk.Frame(marco_contenedor)
marco_segunda_columna.grid(row=0, column=1, sticky="nw")  # `sticky="nw"` para alinear en la parte superior izquierda

# Parte de "Dimensionamiento"
marco_dimensionamiento = tk.Frame(
    marco_segunda_columna,
    highlightthickness=3,
    highlightbackground="gray",
    highlightcolor="lime",
    bg="lightgray"
)
marco_dimensionamiento.grid(row=0, column=0, padx=5, pady=(0, 2), sticky="w")  # Reduce espacio inferior con pady=(0, 2)
ttk.Label(marco_dimensionamiento, text="Dimensionamiento & Discretización", background="lightgreen", font=fuente_titulo).grid(row=0, column=0, columnspan=2, sticky="w")

campos_dimensionamiento = [
    ("Ancho (m)", "5"),
    ("Altura (m)", "10"),
    ("nodos x", "100"),
    ("nodos y", "100"),
    ("Medición (m)", "1")
]
for i, (nombre, valor_defecto) in enumerate(campos_dimensionamiento):
    ttk.Label(marco_dimensionamiento, text=nombre).grid(row=i + 1, column=0, sticky="w", padx=5, pady=2)
    entrada = ttk.Entry(marco_dimensionamiento)
    entrada.insert(0, valor_defecto)
    entrada.grid(row=i + 1, column=1, padx=5, pady=2)
    entradas[nombre] = entrada
    # Botón de ayuda
    boton_ayuda = ttk.Button(marco_dimensionamiento, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
    boton_ayuda.grid(row=i+1, column=2, padx=2, pady=1)
    # Aplicar fuente personalizada al botón
    boton_ayuda.configure(style="Ayuda.TButton")
    # Crear un estilo para el botón
    style = ttk.Style()
    style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido

# Parte de "Iteraciones"
marco_iteracion = tk.Frame(
    marco_segunda_columna,
    highlightthickness=3,
    highlightbackground="gray",
    highlightcolor="orange",
    bg="lightgray"
)
marco_iteracion.grid(row=1, column=0, padx=5, pady=(2, 0), sticky="w")  # Reduce espacio superior con pady=(2, 0)
ttk.Label(marco_iteracion, text="Iteraciones", background="lightcoral", font=fuente_titulo).grid(row=0, column=0, columnspan=2, sticky="w")

campos_iteraciones = [
    ("Max iter rho", "400"),
    ("Max iter V", "230"),
    ("Max iter Gob", "10")
]
for i, (nombre, valor_defecto) in enumerate(campos_iteraciones):
    ttk.Label(marco_iteracion, text=nombre).grid(row=i + 1, column=0, sticky="w", padx=5, pady=2)
    entrada = ttk.Entry(marco_iteracion)
    entrada.insert(0, valor_defecto)
    entrada.grid(row=i + 1, column=1, padx=5, pady=2)
    entradas[nombre] = entrada
    # Botón de ayuda
    boton_ayuda = ttk.Button(marco_iteracion, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
    boton_ayuda.grid(row=i+1, column=2, padx=2, pady=1)
    # Aplicar fuente personalizada al botón
    boton_ayuda.configure(style="Ayuda.TButton")
    # Crear un estilo para el botón
    style = ttk.Style()
    style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido


# Variables específicas
entrada_cantidad = entradas["Subconductores"]
entrada_separacion = entradas["Separación (cm)"]
entrada_posicion_x = entradas["Posición en x (m)"]
entrada_posicion_y = entradas["Posición en y (m)"]

entrada_posicion_x2 =  entradas["Posición en x 2 (m)"]
entrada_posicion_y2 =  entradas["Posición en y 2 (m)"]


# Registrar la función de validación
validacion_natural = ventana.register(lambda texto, nombre: validar_natural(texto, nombre))
validacion_uno = ventana.register(lambda texto, nombre: validar_uno(texto, nombre))
# Configurar validaciones con nombres específicos
entrada_cantidad.config(validate="key", validatecommand=(validacion_natural, "%P", "Subconductores"))
entradas["Max iter rho"].config(validate="key", validatecommand=(validacion_natural, "%P", "Max iter rho"))
entradas["Max iter V"].config(validate="key", validatecommand=(validacion_natural, "%P", "Max iter V"))
entradas["Max iter Gob"].config(validate="key", validatecommand=(validacion_natural, "%P", "Max iter Gob"))
entradas["factor conductor"].config(validate="key", validatecommand=(validacion_uno, "%P", "factor conductor"))
entradas["factor conductor 2"].config(validate="key", validatecommand=(validacion_uno, "%P", "factor conductor 2"))

entrada_cantidad.bind("<KeyRelease>", verificar_separacion)  # Asociar el evento de escritura
entrada_separacion.config(state="disabled")  # Desactivada inicialmente

# Validación coeficiente rugosidad piso para  viento
entrada_modo = entradas["Modo (str)"]
entrada_rugosidad_terreno = entradas["Rugosidad terreno"]
entrada_modo.bind("<<ComboboxSelected>>", verificar_modo_viento) # Asociar evento de selección al Combobox
entrada_rugosidad_terreno.config(state='disabled') # Desactivada inicialmente

# Botones principales
marco_botones = ttk.Frame(ventana, padding=10)
marco_botones.grid(row=len(secciones) + 1, column=0, columnspan=3, sticky="w", padx=10, pady=5)  # Usar grid

# Coloca los botones principales en las primeras columnas
ttk.Button(marco_botones, text="Guardar", command=guardar_datos).grid(row=0, column=0, padx=5)
ttk.Button(marco_botones, text="Limpiar", command=limpiar_campos).grid(row=0, column=1, padx=5)
ttk.Button(marco_botones, text="Graficar", command=validar_campos_y_graficar).grid(row=0, column=2, padx=5)
ttk.Button(marco_botones, text="Ejecutar", command=validar_campos_y_ejecutar).grid(row=0, column=3, padx=5)
ttk.Button(marco_botones, text="Salir", command=cerrar_programa).grid(row=0, column=4, padx=5)
ventana.protocol("WM_DELETE_WINDOW", cerrar_programa)  # Para el botón de cerrar ventana


# Crear el botón Bipolar inmediatamente a la derecha
boton_bipolar_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_bipolar = ttk.Checkbutton(marco_botones, text="Bipolar", variable=boton_bipolar_var, command=activar_bipolar)
boton_bipolar.grid(row=0, column=5, padx=5, pady=10)  # Ubicado justo después de los botones principales

# Añadir un nuevo botón "Auto red" inmediatamente a la derecha del Bipolar
boton_auto_red_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_auto_red = ttk.Checkbutton(marco_botones, text="Auto red", variable=boton_auto_red_var, command=activar_auto_red)
boton_auto_red.grid(row=0, column=6, padx=5, pady=10)  # Ubicado justo al lado del Bipolar


# Botones de gráficos
marco_contenedor = tk.Frame(ventana)
marco_contenedor.grid(row=len(secciones) + 3, column=0, columnspan=2, padx=10, pady=5)
'''
# Título destacado
titulo = ttk.Label(
    marco_contenedor,
    text="Gráficos",
    background="lightblue",
    font=("Arial", 9, "bold")
)
titulo.grid(row=0, column=0, sticky="w", padx=5, pady=(0, 5))
'''
# Marco decorado
marco_graficos = tk.Frame(
    marco_contenedor,
    highlightthickness=3,           # Grosor del borde
    highlightbackground="gray",    # Color del borde sin enfoque
    highlightcolor="green",        # Color del borde con enfoque
    bg="lightgray"                 # Fondo del marco
)
marco_graficos.grid(row=1, column=0, sticky="w", padx=10, pady=5)
ttk.Label(marco_graficos, text="Gráficos", background="lightblue", font=fuente_titulo).grid(row=0, column=0, sticky="w")
# Botón de ayuda
boton_ayuda = ttk.Button(
    marco_graficos,
    text="?",
    width=2,
    command=lambda: mostrar_explicacion("Gráficos disponibles"))  # Usar lambda para diferir la ejecución)
boton_ayuda.grid(row=0, column=1, padx=2, pady=1)
# Aplicar fuente personalizada al botón
boton_ayuda.configure(style="Ayuda.TButton")
# Crear un estilo para el botón
style = ttk.Style()
style.configure("Ayuda.TButton", font=fuente_boton, padding=(2, 2))  # Padding interno reducido
boton_Spd = None  # Variable para almacenar el botón de SPd
boton_Rhop = None  # Variable para almacenar el botón de SPd
boton_Rhon = None  # Variable para almacenar el botón de SPd
# Distribuir botones en 2 filas y 4 columnas
num_columnas = 6  # Número de columnas deseadas
for i, grafico in enumerate(estados_graficos):
    fila = i // num_columnas  # División entera para obtener la fila
    columna = i % num_columnas  # Resto de la división para obtener la columna
    boton = ttk.Button(marco_graficos, text=grafico, command=lambda g=grafico: alternar_grafico(g))
    boton.grid(row=fila+1, column=columna, padx=5, pady=5)
    # Guardar referencia al botón "SPd"
    if grafico == "SPd":
        boton_Spd = boton
        boton_Spd.config(state=tk.DISABLED)  # Iniciar deshabilitado
    if grafico == "Rhop":
        boton_Rhop = boton
        boton_Rhop.config(state=tk.DISABLED)  # Iniciar deshabilitado
    if grafico == "Rhon":
        boton_Rhon = boton
        boton_Rhon.config(state=tk.DISABLED)  # Iniciar deshabilitado

# Iniciar la aplicación
inicializar_grafico()
ventana.mainloop()










