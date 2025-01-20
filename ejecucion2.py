import tkinter as tk
from tkinter import ttk, messagebox, font
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math
import numpy as np
import json
import subprocess
import math as ma
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import math

plt.close('all')

# Estado de gráficos
#'Eele', 'Vele', 'Rhof', 'Vf', 'Vdef', 'Ef', 'J1', 'E1', 'SubPl'
estados_graficos = {
    "Eele": False, "Vele": False, "Rhof": False, "Vf": False, "SPd": False, "Rhop": False,
    "Vdef": False, "Ef": False, "J1": False, "E1": False, "SPv": False, "Rhon": False
}
graficos = []
estado_guardado = False
estado_area = False
todos_seleccionados = False
def alternar_grafico(nombre, mostrar_mensaje=True):
    """
    Alterna el estado de un gráfico, actualiza la lista de gráficos seleccionados y muestra un mensaje opcional.
    
    Args:
        nombre (str): Nombre del gráfico.
        mostrar_mensaje (bool): Indica si se debe mostrar un mensaje de confirmación.
    """
    # Cambiar el estado del gráfico
    estados_graficos[nombre] = not estados_graficos[nombre]
    mensaje = "seleccionado" if estados_graficos[nombre] else "deseleccionado"
    
    # Actualizar la lista de gráficos activados
    if estados_graficos[nombre]:
        if nombre not in graficos:
            graficos.append(nombre)
    else:
        if nombre in graficos:
            graficos.remove(nombre)
    
    # Mostrar el mensaje si está habilitado
    if mostrar_mensaje:
        messagebox.showinfo(nombre, f"Se ha {mensaje} {nombre}.")
    
    # Imprimir el estado actual de la lista (opcional)
    print(f"Gráficos seleccionados: {graficos}")

def Hay_corona(Vol, Ev, ev):
    '''
    Vol: en kV
    Sep: en cm
    r: en cm
    '''
    if ev <= np.abs(Vol):
        print(f"{Vol} kV >= {ev} kV, hay efecto corona")
        return True
    else:
        print(f"{Vol} kV < {ev} kV, no hay efecto corona")
        return False
    
def grad_sup(g0, m, delta, r):
    Ev = g0*delta*m*(1+0.301/np.sqrt(delta*r))
    return  Ev

def Vol_crit(Ev, Sep, r):
    print(f"los parametros para Vol_crit son", Ev, Sep, r)
    ev = Ev*r*np.log(Sep/r)
    return ev

def pressure_height(altura,tempe):
    P0  = 101325 # Pa
    M = 0.029 # Kg/mol
    g =  9.8 # m/s^2
    R0 = 8.314 # J/mol K
    P = P0*np.e**(-g*M*altura/(R0*tempe))
    return P/1000 # en kPa

def seleccionar_todos():
    """
    Alterna la selección de todos los gráficos.
    """
    global todos_seleccionados
    
    for grafico in estados_graficos:
        if grafico in ["SPd", "Rhop", "Rhon"] and not boton_bipolar_var.get():
            # No modificar gráficos deshabilitados
            continue
        alternar_grafico(grafico, mostrar_mensaje=False)
    
    # Mostrar mensaje de confirmación
    mensaje_completo = "seleccionados" if not todos_seleccionados else "deseleccionados"
    messagebox.showinfo("Selección Completa", f"Todos los gráficos han sido {mensaje_completo}.")
    
    # Alternar el estado global
    todos_seleccionados = not todos_seleccionados
####################################
# FUNCIONES DE VERIFICACION DE ESPACIO
######################################
def obtener_parametro(param, tipo):
    """
    Verifica si un parámetro está vacío y lo asigna como None o lo convierte al tipo especificado.

    :param param: Valor del parámetro a verificar.
    :param tipo: Tipo al que se debe convertir (float o int).
    :return: Valor convertido o None si está vacío.
    """
    if param == "" or param is None:
        return None
    return tipo(param)

def convierte_a_radio(area, conversion=1, es_mcm=False):
    """
    Calcula el radio de un conductor a partir del área.
    Parámetros:
        area: Área del conductor.
        conversion: Factor de conversión (por ejemplo, de MCM a mm²).
        es_mcm: Indica si el área está en MCM (True) o no (False).
    Retorna:
        Radio en las mismas unidades que el área (por defecto en mm).
    """
    factor = conversion if es_mcm else 1
    return np.sqrt(area * factor / np.pi)

def calcula_lados(numero, separacion):
    """
    Calcula las distancias entre los subconductores según su disposición.
    Parámetros:
        numero: Número de subconductores.
        separacion: Distancia entre subconductores vecinos.
    Retorna:
        Lista con las longitudes de los lados del polígono.
    """
    lados = [separacion] * numero
    if numero == 4:
        lados.extend([separacion * np.sqrt(2)] * 2)
    return lados

def distancia_equivalente(lados):
    """
    Calcula la distancia geométrica media (D_eq) a partir de una lista de lados.
    Parámetros:
        lados: Lista de distancias entre subconductores.
    Retorna:
        Distancia geométrica media (D_eq).
    """
    return (ma.prod(lados))**(1 / len(lados))

def distancia_equivalente_2(sepa, numero):
    """
    Calcula la distancia geométrica media (D_eq) a partir de una lista de lados.
    Parámetros:
        sepa: Separación subconductores vecinos
    Retorna:
        Distancia geométrica media (D_eq).
    """
    return sepa/(2*np.sin(np.pi/numero))

def radio_eq(radio, numero, distancia):
    """
    Calcula el radio equivalente de un subconductor fasciculado.
    Parámetros:
        radio: Radio del subconductor individual.
        numero: Número de subconductores.
        distancia: Distancia equivalente entre subconductores.
    Retorna:
        Radio equivalente del conductor.
    """
    return (radio*numero*distancia**(numero-1))**(1/numero)

def calculo_radio_eq(numero, area_sub, sepa, conversion=1, es_mcm=False, es_cm=False):
    """
    Calcula el radio equivalente de un conductor fasciculado.
    Parámetros:
        numero: Número de subconductores.
        area_sub: Área de cada subconductor.
        sep: Separación entre subconductores vecinos.
        conversion: Factor de conversión del área si está en MCM.
        es_mcm: Indica si el área está en MCM.
        es_cm: Indica si la separación está en cm.
    Retorna:
        Tuple con el radio equivalente y el radio individual.
    """
    #lados = calcula_lados(numero, sepa)
    #distancia = distancia_equivalente(lados)
    distancia_e = distancia_equivalente_2(sepa, numero)
    radio = convierte_a_radio(area_sub, conversion, es_mcm)
    if es_cm:
        radio /= 10  # Convertir de mm a cm si es necesario
    radio_equi = radio_eq(radio, numero, distancia_e)
    if radio_equi==0:
        radio_equi = radio
    return radio_equi, radio

def espacio(anchox,anchoy, nodox, nodoy):
    '''
    Función que entrega las dimensiones espaciales de la malla en base al ancho de separacion de los nodos y en base a
    la cantidad de nodos que se desean en ambas dimensiones
    '''
    sx = np.abs((anchox/2)*(nodox-1))
    sy = np.abs((anchoy)*(nodoy-1))
    return sx, sy

def malla(ancho, Sx, Sy, nodox = None, nodoy = None):
    '''
    Discretiza la malla y define a los  arreglos de distancia que definen la malla
    Si nodox y nodoy no son dados, asume que se calculan a partir de 'ancho'
    Si son dados, solo discretiza el largo x e y
    '''
    if nodox is None and nodoy is None:
        x = np.arange(-Sx, Sx+ancho, ancho)
        y = np.arange(Sy, 0-ancho, -ancho)
        nodox = len(x)
        nodoy = len(y)
    else:
        x = np.linspace(-Sx, Sx, nodox)
        y = np.linspace(Sy, 0, nodoy)
    return x, y, nodox, nodoy

def verificar_discretizacion(min, sx, sy, nodox, nodoy, x_coor, y_coor):
    print(f"tipo sx,sy son {type(sx)}, {type(sy)}")
    if (sx is None and sy is None) and (nodox is not None and nodoy is not None):
        Sx, Sy = espacio(min, min, nodox, nodoy)
        x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=nodox, nodoy=nodoy)
        # Calcular las diferencias y nodos adicionales necesarios
        distx = np.abs(x_coor) - Sx
        disty = np.abs(y_coor) - Sy
        nodos_sx = max(int(distx / min) + 1, 0)
        nodos_sy = max(int(disty / min) + 1, 0)
        # Inicializar un mensaje vacío
        mensajes_error = []
        # Verificar condiciones de dimensionamiento
        if distx > 0:
            mensajes_error.append(
                f"Error de dimensionamiento en X. Si quieres usar {nodox} nodos en el eje horizontal entonces aumenta el radio equivalente o selecciona más nodos: {nodos_sx + nodosx}, "
                f"o bien reduce la distancia x_coor. Diferencia x-sx: {distx:.2f} metros.")
        if disty > 0:
            mensajes_error.append(
                f"Error de dimensionamiento en Y. Si quieres usar {nodoy} nodos en el eje vertical entonces aumenta el radio equivalente o selecciona más nodos: {2 * (nodos_sy + nodosy)}, "
                f"o bien reduce la distancia y_coor. Diferencia y-sy: {disty:.2f} metros.")
        # Mostrar mensajes de error si hay problemas
        if mensajes_error:
            messagebox.showinfo(
                "Verificación dimensionamiento", "\n".join(mensajes_error))
            return False
        # Si no hay problemas, retornar True
        return True
    return True

# Botón Bipolar para activar/desactivar las casillas extra
def activar_bipolar():
    # Activar o desactivar las casillas extras
    if boton_bipolar_var.get():  # Si se selecciona el botón Bipolar
        menu_tipo2.configure(state="normal")
        menu_nombre2.configure(state="normal")
        entradas["Voltaje 2 (kV)"].config(state="normal")
        #entradas["Área 2 (mm²)"].config(state="normal")
        entradas["Posición en x 2 (m)"].config(state="normal")
        entradas["Posición en y 2 (m)"].config(state="normal")
        entradas["Estado sup cond 2"].config(state="normal")
        entradas["Movilidad iónica 2 (m2/kVs)"].config(state="normal")
        entradas["Recombinación (μm^3/s)"].config(state="normal")
        entradas["Jp (nA/m^2)"].config(state="disabled")
        # Verifica si alguna casilla está activada
        boton_Spd.config(state=tk.NORMAL)
        boton_Rhop.config(state=tk.NORMAL)
        boton_Rhon.config(state=tk.NORMAL)
    else:  # Si el botón Bipolar está desmarcado
        menu_tipo2.configure(state="disabled")
        menu_nombre2.configure(state="disabled")
        # Eliminar el contenido actual de las casillas
        entradas["Voltaje 2 (kV)"].config(state="disabled")
        #entradas["Área 2 (mm²)"].config(state="disabled")
        entradas["Posición en x 2 (m)"].config(state="disabled")
        entradas["Posición en y 2 (m)"].config(state="disabled")
        entradas["Estado sup cond 2"].config(state="disabled")
        entradas["Movilidad iónica 2 (m2/kVs)"].config(state="disabled")
        entradas["Recombinación (μm^3/s)"].config(state="disabled")
        entradas["Jp (nA/m^2)"].config(state="normal")

        boton_Spd.config(state=tk.DISABLED)
        boton_Rhop.config(state=tk.DISABLED)
        boton_Rhon.config(state=tk.DISABLED)
        
        # Si ya había un gráfico previo, eliminar los puntos asociados a x2, y2
        plt.close('all')  # Cierra el gráfico previo
        # Redibujar solo los puntos x, y (sin x2, y2)
        validar_campos_y_graficar()  # Esto vuelve a graficar los puntos con las nuevas configuraciones

def activar_desarrollador():
    if boton_desarrollador_var.get():
        entradas["Max iter rho"].config(state="normal")
        entradas["Max iter V"].config(state="normal")
        entradas["Max iter Gob"].config(state="normal")
        entradas["Interior conductor"].config(state="normal")
    else:
        entradas["Max iter rho"].config(state="disabled")
        entradas["Max iter V"].config(state="disabled")
        entradas["Max iter Gob"].config(state="disabled")
        entradas["Interior conductor"].config(state="disabled")

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


def obtener_area_escogida(entrada_area):
    # Obtener la clave y color seleccionados para el área
    llave = str(entrada_area[0].get())
    color = str(entrada_area[1].get())
    lista_opciones = opciones_subcon[llave]
    # Buscar el valor correspondiente al color seleccionado
    for color_opcion, valor1, valor2 in lista_opciones:
        if color_opcion == color:
            val = [valor1, valor2]
            return val # area y diametro mm^2 y mm
    return None  # En caso de que no se encuentre el color

import json
from tkinter import messagebox
import subprocess

def procesar_entradas():
    params = {}
    for nombre, entrada in entradas.items():
        if isinstance(entrada, list):
            # Si la entrada es una lista, obtener los valores serializables
            params[nombre] = [elem.get() if hasattr(elem, 'get') else elem for elem in entrada]
        elif hasattr(entrada, 'get'):  
            # Si es un widget con 'get' (Entry, Combobox, etc.), usar su valor
            params[nombre] = entrada.get()
        else:
            # Si es otro tipo (cadena, entero, etc.), asignarlo directamente
            params[nombre] = entrada
    # Agregar otros parámetros generales
    params['graficos'] = graficos
    params['guardar'] = estado_guardado
    return params

def ejecutar_script():
    # Procesar el área principal
    # Procesar todas las entradas y convertirlas en un diccionario
    params = procesar_entradas()

    # Mostrar los parámetros en un cuadro de diálogo
    mensaje = "\n".join([f"{clave}: {valor}" for clave, valor in params.items()])
    messagebox.showinfo("Parámetros Guardados", mensaje)

    try:
        # Convertir parámetros a JSON
        params_json = json.dumps(params)
        
        # Ejecutar el script correspondiente
        script = "script_ejecucion_bip.py" if boton_bipolar_var.get() else "script_ejecucion_uni.py"
        subprocess.run(["python", script, params_json], check=True)
        
        # Notificar éxito
        messagebox.showinfo("Éxito", f"El archivo {script} se ejecutó correctamente.")
    except Exception as e:
        # Notificar error
        messagebox.showerror("Error", f"Ocurrió un error al ejecutar el script: {str(e)}")


def process_area(area, cantidad, sep):
    Req, R = calculo_radio_eq(cantidad, area, sep, conversion=0.5067, es_mcm=False, es_cm=True) # están en cm
    print(f"radio eq {Req} cm y radio subconductor {R} cm")
    Req /= 100 # en m
    R /= 100
    print(f'{Req},{R} en metros')
    return Req, R # devuelve en metros

def verificar_selecciones():
    if not tipo_conductor.get():
        return "El tipo de conductor no está seleccionado."
    if not nombre_conductor.get():
        return "El nombre del conductor no está seleccionado."
    if boton_bipolar_var.get():
        if not tipo_conductor2.get():
            return "El tipo de conductor 2 no está seleccionado."
        if not nombre_conductor2.get():
            return "El nombre del conductor 2 no está seleccionado."
    return None
def mostrar_error(mensaje):
    messagebox.showwarning("Error", mensaje)

def extrae_radios(cantidad, sepa):
    # devuelve [R1, 10**10] si unipolar
    # devuelve [R1, R2] si bipolar
    radios = []
    area1 = entradas["Área (mm²)"][2][0]
    Req, R = process_area(area1, cantidad, sepa) # en metros
    radios.append(Req)
    if boton_bipolar_var.get():
        area2 = entradas["Área 2 (mm²)"][2][0]
        Req2, R2 = process_area(area2, cantidad, sepa)
        radios.append(Req2)
    else:
        radios.append(10**10)
    return radios

def validar_campos_y_ejecutar():
    global estado_area
    # Obtener valor de cantidad
    cantidad = entrada_cantidad.get()
    try:
        cantidad = int(cantidad)
        print(f"cantidad de subconductores es {cantidad}")
    except ValueError:
        cantidad = 0
        mostrar_error("La cantidad de subconductores debe ser un número entero.")
        return
    sx = obtener_parametro(entradas["Ancho (m)"].get(), float)  # (m) Media longitud del plano de tierra
    sy = obtener_parametro(entradas["Altura (m)"].get(), float)  # (m) Altura del área de estudio respecto de tierra
    nodox = obtener_parametro(entradas["nodos x"].get(), int)  # Número de nodos en x
    nodoy = obtener_parametro(entradas["nodos y"].get(), int)  # Número de nodos en y
    try:
        sep = float(entradas["Separación (cm)"].get())
        print(f"la separación es {sep}")
    except:
        messagebox.showwarning("Advertencia", "La separación debe ser un número.")
        return
    radios =  []
    for nombre, entrada in entradas.items():
        #print(f"la entrada es {entrada}")
        if nombre == "Separación (m)" and cantidad <= 1:
           continue  # Ignorar esta casilla
        elif nombre == "Subconductores":
            continue
        elif nombre == "Rugosidad terreno":
            continue
        elif nombre == "Ancho (m)":
            #sx = entradas[nombre].get()
            continue  # Ignorar esta casilla
        elif nombre == "Altura (m)":
            #sy = entradas[nombre].get()
            continue
        elif nombre == "nodos x":
            #nodox = entradas[nombre].get()
            continue  # Ignorar esta casilla
        elif nombre == "nodos y":
            #nodoy = entradas[nombre].get()
            continue
        elif nombre == "Posición en x (m)":
            x1 = float(entradas[nombre].get().strip())  # Obtener el valor de la casilla "Posición en x (m)"
            x2 = float(entradas["Posición en x 2 (m)"].get().strip())  # Obtener el valor de la casilla "Posición en x 2 (m)"
            if x1 - x2 == 0 and boton_bipolar_var.get():
                messagebox.showwarning(
                    "Advertencia", "En modo Bipolar, los conductores no pueden estar en la misma posición"
                )
                return
        elif nombre == "Onset corona (kV)":
            continue
        elif nombre == "Corona gradient (kV/cm)":
            continue
        elif nombre == "Posición en y (m)":
            y1 = float(entradas[nombre].get())
        
        elif nombre == "Área (mm²)":
            if estado_area:
                radios = extrae_radios(cantidad, sep)
            else:
                verificar_estado_area()
                radios = extrae_radios(cantidad, sep)
            estado_area = False
            continue
        elif nombre == "Área 2 (mm²)":
            continue
        else:
            valor = entradas[nombre].get().strip()  # Obtener valor sin espacios adicionales
        print(f"el valor de {nombre} es {valor}")
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
    print(f"los radios son {radios}")
    print(f"el sx y sy son {sx} y {sy}")
    print(f"el nodox y nodoy son {nodox} y {nodoy}")
    verificado = verificar_discretizacion(np.min(radios), sx, sy, nodox, nodoy, x1, y1)
    if not verificado:
        return
    messagebox.showinfo("Validación exitosa", "Todos los campos han sido completados correctamente.")
    '''
    Vol1 = entradas["Voltaje (kV)"].get() # kV
    m1 = entradas["Estado sup cond 1"].get()
    Vol2 = entradas["Voltaje 2 (kV)"].get() # kV
    m2 = entradas["Estado sup cond 2"].get()
    delta = param_corona()
    Von1, Egrad1 = calcula_corona(sep, radios[0], delta, m1, g0) # Se crean los elementos Ev, ev en el diccionarios
    Von2, Egrad2 = calcula_corona(sep, radios[1], delta, m2, g0) 
    corona1 = Hay_corona(Vol1, Egrad1, Von1)
    corona2 = Hay_corona(Vol2, Egrad2, Von2)
    if corona1 or corona2:
        messagebox.showinfo("Advertencia", f"Se verifica que si hay efecto corona en los conductores.\n"
                            f"Voltaje aplicado {Vol} kV >= Voltaje crítico {Von} kV y gradiente superficial {Egrad} kV/cm")
    else:
        messagebox.showinfo("Advertencia", "Se verifica que no hay efecto corona en los conductores.\n"
                            f"Voltaje aplicado {Vol} kV < Voltaje crítico {Von} kV y gradiente superficial {Egrad} kV/cm")
    '''
    y1 = float(entradas["Posición en y (m)"].get())
    if boton_bipolar_var.get():
        y2 = float(entradas["Posición en y 2 (m)"].get())
    else:
        y2 = 0
    alturas = np.array([y1, y2])
    verificar_corona(alturas*100, np.array(radios)*100) # en centimetros el radio
    
    ejecutar_script()


def verificar_corona(alturas, radios):
    # Obtener valores desde el diccionario `entradas`
    Vol1 = float(entradas["Voltaje (kV)"].get())  # kV
    m1 = float(entradas["Estado sup cond 1"].get())
    print(f"el radio1 es {radios[0]}")
    print(f"el radio2 es {radios[1]}")
    delta = param_corona()
    
    # Calcular para el primer conductor
    Von1, Egrad1 = calcula_corona(alturas[0], radios[0], delta, m1, g0)
    corona1 = Hay_corona(Vol1, Egrad1, Von1)
    
    # Detalles del efecto corona
    detalles = f"Conductor 1:\n"
    if corona1:
        detalles += (f"Voltaje aplicado: {Vol1} kV >= Voltaje crítico: {Von1:.2f} kV\n"
                     f"Gradiente superficial: {Egrad1:.2f} kV/cm\n\n")
    else:
        detalles += (f"Voltaje aplicado: {Vol1} kV < Voltaje crítico: {Von1:.2f} kV\n"
                     f"Gradiente superficial: {Egrad1:.2f} kV/cm\n\n")
    
    # Verificar si se deben calcular los valores para el segundo conductor
    if boton_bipolar_var.get():
        Vol2 = float(entradas["Voltaje 2 (kV)"].get())  # kV
        m2 = float(entradas["Estado sup cond 2"].get())
        
        # Calcular para el segundo conductor
        Von2, Egrad2 = calcula_corona(alturas[1], radios[1], delta, m2, g0)
        corona2 = Hay_corona(Vol2, Egrad2, Von2)
        
        detalles += f"Conductor 2:\n"
        if corona2:
            detalles += (f"Voltaje aplicado: {Vol2} kV >= Voltaje crítico: {Von2:.2f} kV\n"
                         f"Gradiente superficial: {Egrad2:.2f} kV/cm\n\n")
        else:
            detalles += (f"Voltaje aplicado: {Vol2} kV < Voltaje crítico: {Von2:.2f} kV\n"
                         f"Gradiente superficial: {Egrad2:.2f} kV/cm\n\n")
    else:
        detalles += "El cálculo para el conductor 2 no fue realizado.\n"
    
    # Mostrar el mensaje de advertencia
    if corona1 or (boton_bipolar_var.get() and corona2):
        messagebox.showinfo("Advertencia", f"Se verifica que sí hay efecto corona en los conductores.\n\n{detalles}")
    else:
        messagebox.showinfo("Advertencia", f"Se verifica que no hay efecto corona en los conductores.\n\n{detalles}")


P0 = 101.3  # kPa
T0 = 298.15  # Kelvin
g0 = 29.8 # kV/cm
def param_corona():
    altura = float(entradas["Altitud (m)"].get()) # m
    Tem = float(entradas["Temperatura (°C)"].get()) + 273.15 # Kelvin
    Pres = pressure_height(altura, Tem)  # kPa
    delta = Pres * T0 / (P0 * Tem)  # Densidad del aire
    return delta
    


def calcula_corona(Sep, r, delta, m, g0):
    '''
    Vol en kV
    Sep en cm
    r en cm
    g0 en kV/cm
    '''
    Ev= grad_sup(g0, m, delta, r)
    ev = Vol_crit(Ev, Sep, r)
    entradas["Onset corona (kV)"] = ev
    entradas["Corona gradient (kV/cm)"] = Ev
    return ev, Ev

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
    v1 = float(entradas["Voltaje (kV)"].get())
    v2 = float(entradas["Voltaje 2 (kV)"].get())
    volts = [v1, v2] if boton_bipolar_var.get() else [v1]
    # Graficar puntos iniciales
    graficar_puntos(fig, ax, xs, ys, cantidad, separacion, ancho_mitad, altura_total, volts)


def validar_campos_y_graficar():
    # Cerrar todas las figuras de Matplotlib previas si existen
    plt.close('all')  # Asegura que las figuras de Matplotlib se cierren correctamente
    v1 = float(entradas["Voltaje (kV)"].get())
    v2 = float(entradas["Voltaje 2 (kV)"].get())
    volts = [v1, v2] if boton_bipolar_var.get() else [v1]
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
    graficar_puntos(fig, ax, xs, ys, cantidad, separacion, ancho_mitad, altura_total, volts)

def graficar_puntos(fig, ax, x, y, cantidad, separacion, ancho_mitad, altura_total, volts):
    # Configurar límites
    ax.set_xlim(-ancho_mitad, ancho_mitad)
    ax.set_ylim(0, altura_total)
    ax.set_title("Distribución de conductores")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(True)

    # Normalizar los valores absolutos de volts para usar un colormap, con rango fijo [0, 1500]
    vmin, vmax = 0, 1200
    norm = Normalize(vmin=vmin, vmax=vmax)
    colormaps = [plt.cm.Blues, plt.cm.Greens]  # Mapas de colores para los dos puntos

    n = int(cantidad)
    pos = ["izquierdo", "derecho"]

    # Dibujar puntos
    for j in range(len(x)):
        # Calcular color según el voltaje y el colormap correspondiente
        color = colormaps[j % len(colormaps)](norm(abs(volts[j])))

        # Verificar si los puntos están dentro de los límites
        if abs(x[j]) > ancho_mitad:
            messagebox.showwarning("Advertencia", f"El punto ({x[j]}, {y[j]}) excede el rango horizontal ({-ancho_mitad}, {ancho_mitad}).")
        elif n == 1:
            if len(x) == 1:
                ax.plot(x[j], y[j], 'o', color=color, label=f"Conductor")
            else:
                ax.plot(x[j], y[j], 'o', color=color, label=f"Conductor {pos[j]}")
        elif n >= 2:
            radio = ((separacion / 100) / 2) * (1 / math.sin(math.pi / n))
            puntos = generar_poligono(n, radio, x[j], y[j])
            for i in range(len(puntos)):
                ax.plot(puntos[i][0], puntos[i][1], 'o', color=color, label=f"Conductor {pos[j]}" if i == 0 else None)

    # Agregar leyenda
    ax.legend(loc="upper right")

    # Agregar barras de colores personalizadas
    for j in range(len(volts)):
        cmap = colormaps[j % len(colormaps)]  # Seleccionar el colormap correspondiente
        ticks = [0, 400, 800, 1200] if volts[j] >= 0 else [0, -400, -800, -1200]
        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, location="right", pad=0.1)
        cbar.set_label(f"Voltaje {pos[j]} (V)", rotation=270, labelpad=15)
        cbar.set_ticks([abs(t) for t in ticks])  # Posiciones basadas en valores absolutos
        cbar.ax.set_yticklabels([str(t) for t in ticks])  # Etiquetas reflejan el signo del voltaje
    # Ajustar el layout para evitar solapamientos
    fig.tight_layout()
    # Cambiar el tamaño de la figura
    fig.set_size_inches(6, 5)  # Cambiar el tamaño de la figura a 12x8 pulgadas (puedes ajustar estos valores)
    # Ajustar tamaño de la ventana en Tkinter
    ventana.geometry("1300x1100")  # Cambiar el tamaño de la ventana a 1300x900 píxeles
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

def verificar_estado_area():
    # Actualizar el estado de la variable global estado_area
    # También actuliza los estados de area1 y area2
    global estado_area
    #namei == "Área (mm²)"
    error = verificar_selecciones()
    if error:
        mostrar_error(error) # mostrará tanto los errores para el conductor 1  y el 2
        return error
    tipo = tipo_conductor.get()
    nombre_con = nombre_conductor.get()
    # Obtener los valores del diccionario
    datos_conductor = next(
        (x[1:] for x in opciones_subcon[tipo] if x[0] == nombre_con), None
    )
    if datos_conductor:
        # Crear o actualizar el área en el diccionario `entradas`
        entradas["Área (mm²)"] = [tipo, nombre_con, list(datos_conductor)]
    else:
        mostrar_error("Error al buscar los datos del conductor seleccionado.")
        return
    if boton_bipolar_var.get():
        name = "Área 2 (mm²)"
        tipo2 = tipo_conductor2.get()
        nombre_con2 = nombre_conductor2.get()
        # Obtener los valores del diccionario
        datos_conductor2 = next(
            (x[1:] for x in opciones_subcon[tipo2] if x[0] == nombre_con2), None
        )
        if datos_conductor2:
            # Crear o actualizar el área en el diccionario `entradas`
            entradas["Área 2 (mm²)"] = [tipo2, nombre_con2, list(datos_conductor2)]
            area2 = list(datos_conductor2)[0]
            print(f"Entradas actualizadas: {entradas}")
        else:
            mostrar_error("Error al buscar los datos del conductor 2 seleccionado.")
            return
    estado_area = True


def guardar_datos():
    global estado_guardado  # Referencia la variable global
    advertencia = verificar_estado_area()  # Verificar el estado del área antes de guardar
    # Alternar entre True y False
    if advertencia:
        return
    else:
        if estado_guardado:
            estado_guardado = False
            messagebox.showinfo("Datos no guardados", "No se guardarán los datos una vez ejecutado el programa")
        else:
            # Crear un diccionario con los valores capturados
            datos = {}
            for nombre, entrada in entradas.items():
                if isinstance(entrada, ttk.Entry):  # Si es un Entry, capturamos el valor con .get()
                    datos[nombre] = entrada.get()
                elif isinstance(entrada, (tuple, list)):  # Si es una tupla o lista, verificamos sus elementos
                    valores = []
                    for item in entrada:
                        if hasattr(item, "get"):  # Si tiene el método get, lo usamos para capturar el valor
                            valores.append(item.get())
                        else:
                            valores.append(item)  # De lo contrario, añadimos el valor directamente
                    datos[nombre] = type(entrada)(valores)  # Usamos el tipo original (tuple o list)
                else:
                    datos[nombre] = entrada  # Capturamos valores directos para otros tipos
            
            # Mostrar los datos capturados
            messagebox.showinfo("Datos guardados", f"Datos capturados:\n{datos}")
            estado_guardado = True


def limpiar_campos():
    for entrada in entradas.values():
        # Verificar si es una lista o tupla
        if isinstance(entrada, (list, tuple)):
            for sub_entrada in entrada:
                sub_entrada.delete(0, tk.END)
        else:
            # Si no es una lista o tupla, simplemente limpiar el campo
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
            plt.close('all')
            ventana.destroy()


# Crear la ventana principal
ventana = tk.Tk()
ventana.title("Ajuste de parámetros")
ventana.geometry("1300x700") # ajuste tamaño ventana completa

# Diccionario con explicaciones de los parámetros
explicaciones = {
    "Voltaje (kV)": "Voltaje aplicado al conductor, medido en kilovoltios (kV).",
    "Área (mm²)": (
        "Sección transversal del conductor, medida en milímetros cuadrados:\n"
        "- ACSR: Aluminum Conductor Steel Reinforced:\n"
        "    Kiwi: área: 1170 mm^2, diámetro: 44.12 mm \n"
        "    Chukar: área: 976.7 mm^2, diámetro: 40.69 mm \n"
        "    Falcon: área: 908.7 mm^2, diámetro: 39.24 mm \n"
        "    Lapwing: área: 859.8 mm^2, diámetro: 38.15 mm \n"
        "    Plover: área: 818.7 mm^2, diámetro: 37.21 mm \n"
        "    Pheasant: área: 726.8 mm^2, diámetro: 35.10 mm \n"
        "    Bunting: área: 647.7 mm^2, diámetro: 33.07 mm \n"
        "    Curlew: área: 593.6 mm^2, diámetro: 31.65 mm \n"
        "    Rail: área: 517.3 mm^2, diámetro: 29.59 mm \n"
        "    Condor: área: 454.5 mm^2, diámetro: 27.76 mm \n"
        "    Crow: área: 408.5 mm^2, diámetro: 26.28 mm \n"
        "    Gull: área: 361.0 mm^2, diámetro: 25.38 mm \n"
        "    Kingbird: área: 341.0 mm^2, diámetro: 23.88 mm \n"
        "    Osprey: área: 298.2 mm^2, diámetro: 22.33 mm \n"
        "    Pelican: área: 255.1 mm^2, diámetro: 20.68 mm \n"
        "    Oriole: área: 210.3 mm^2, diámetro: 18.82 mm \n"
        "    Ostrich: área: 176.7 mm^2, diámetro: 17.27 mm \n"
        "    Penguin (4/0): área: 125.1 mm^2, diámetro: 14.30 mm \n"
        "    Pigeon (3/0): área: 99.2 mm^2, diámetro: 12.75 mm \n"
        "    Raven (1/0): área: 62.4 mm^2, diámetro: 10.11 mm \n"
        "    Sparrow: área: 39.2 mm^2, diámetro: 8.03 mm \n"
        "    Turkey: área: 15.5 mm^2, diámetro: 5.03 mm \n"
        "\n"
        "- AAC: All Aluminum Conductor:\n"
        "    Cowslip: área: 1013 mm^2, diámetro: 41.41 mm \n"
        "    Jessamine: área: 887 mm^2, diámetro: 38.72 mm \n"
        "    Carnation: área: 725 mm^2, diámetro: 35.02 mm \n"
        "    Narcissus: área: 645 mm^2, diámetro: 32.94 mm \n"
        "    Marigold: área: 564 mm^2, diámetro: 30.88 mm \n"
        "    Goldenrod: área: 483 mm^2, diámetro: 28.60 mm \n"
        "    Lilac: área: 403 mm^2, diámetro: 26.11 mm \n"
        "    Flag: área: 355 mm^2, diámetro: 24.48 mm \n"
        "    Hyacinth: área: 253 mm^2, diámetro: 20.66 mm \n"
        "    Canna: área: 202 mm^2, diámetro: 18.38 mm \n"
        "    Laurel: área: 135 mm^2, diámetro: 15.05 mm \n"
        "    Oxlip: área: 107 mm^2, diámetro: 13.25 mm \n"
        "    Poppy: área: 53.7 mm^2, diámetro: 9.36 mm \n"
        "\n"
        "- ACAR: Aluminum Conductor Alloy Reinforced:\n"
        "    ACAR1: área: 1013 mm^2, diámetro: 41.41 mm \n"
        "    ACAR2: área: 963 mm^2, diámetro: 40.37 mm \n"
        "    ACAR3: área: 887 mm^2, diámetro: 38.72 mm \n"
        "    ACAR4: área: 811 mm^2, diámetro: 37.04 mm \n"
        "    ACAR5: área: 760 mm^2, diámetro: 35.85 mm \n"
        "    ACAR6: área: 659 mm^2, diámetro: 33.37 mm \n"
        "    ACAR7: área: 557 mm^2, diámetro: 30.65 mm \n"
        "    ACAR8: área: 481 mm^2, diámetro: 28.48 mm \n"
        "    ACAR9: área: 380 mm^2, diámetro: 25.32 mm \n"
        "    ACAR10: área: 279 mm^2, diámetro: 21.67 mm \n"
        "    ACAR11: área: 177 mm^2, diámetro: 17.24 mm \n"
        "    ACAR12: área: 125 mm^2, diámetro: 14.31 mm \n"
        "    ACAR13: área: 62.5 mm^2, diámetro: 10.11 mm \n"
    ),
    "Posición en x (m)": "Posición horizontal del conductor 1 en metros.",
    "Posición en y (m)": "Posición vertical del conductor 1 en metros.",
    "Posición en x 2 (m)": "Posición horizontal del conductor 2 en metros.",
    "Posición en y 2 (m)": "Posición vertical del conductor 2 en metros.",
    "Estado sup cond 1": "Estado de superficie para el conductor 1. Factor 'm'.",
    "Estado sup cond 2": "Estado de superficie para el conductor 2. Factor 'm'.",
    "Subconductores": "Número de subconductores usados en cada polo.",
    "Separación (cm)": "Distancia entre subconductores, medida en centímetros.",
    "Voltaje 2 (kV)": "Voltaje aplicado al conductor secundario en kilovoltios (kV).",
    "Área 2 (mm²)": (
        "Sección transversal del segundo conductor, medida en milímetros cuadrados:\n"
        "- ACSR: Aluminum Conductor Steel Reinforced:\n"
        "    Kiwi: área: 1170 mm^2, diámetro: 44.12 mm \n"
        "    Chukar: área: 976.7 mm^2, diámetro: 40.69 mm \n"
        "    Falcon: área: 908.7 mm^2, diámetro: 39.24 mm \n"
        "    Lapwing: área: 859.8 mm^2, diámetro: 38.15 mm \n"
        "    Plover: área: 818.7 mm^2, diámetro: 37.21 mm \n"
        "    Pheasant: área: 726.8 mm^2, diámetro: 35.10 mm \n"
        "    Bunting: área: 647.7 mm^2, diámetro: 33.07 mm \n"
        "    Curlew: área: 593.6 mm^2, diámetro: 31.65 mm \n"
        "    Rail: área: 517.3 mm^2, diámetro: 29.59 mm \n"
        "    Condor: área: 454.5 mm^2, diámetro: 27.76 mm \n"
        "    Crow: área: 408.5 mm^2, diámetro: 26.28 mm \n"
        "    Gull: área: 361.0 mm^2, diámetro: 25.38 mm \n"
        "    Kingbird: área: 341.0 mm^2, diámetro: 23.88 mm \n"
        "    Osprey: área: 298.2 mm^2, diámetro: 22.33 mm \n"
        "    Pelican: área: 255.1 mm^2, diámetro: 20.68 mm \n"
        "    Oriole: área: 210.3 mm^2, diámetro: 18.82 mm \n"
        "    Ostrich: área: 176.7 mm^2, diámetro: 17.27 mm \n"
        "    Penguin (4/0): área: 125.1 mm^2, diámetro: 14.30 mm \n"
        "    Pigeon (3/0): área: 99.2 mm^2, diámetro: 12.75 mm \n"
        "    Raven (1/0): área: 62.4 mm^2, diámetro: 10.11 mm \n"
        "    Sparrow: área: 39.2 mm^2, diámetro: 8.03 mm \n"
        "    Turkey: área: 15.5 mm^2, diámetro: 5.03 mm \n"
        "\n"
        "- AAC: All Aluminum Conductor:\n"
        "    Cowslip: área: 1013 mm^2, diámetro: 41.41 mm \n"
        "    Jessamine: área: 887 mm^2, diámetro: 38.72 mm \n"
        "    Carnation: área: 725 mm^2, diámetro: 35.02 mm \n"
        "    Narcissus: área: 645 mm^2, diámetro: 32.94 mm \n"
        "    Marigold: área: 564 mm^2, diámetro: 30.88 mm \n"
        "    Goldenrod: área: 483 mm^2, diámetro: 28.60 mm \n"
        "    Lilac: área: 403 mm^2, diámetro: 26.11 mm \n"
        "    Flag: área: 355 mm^2, diámetro: 24.48 mm \n"
        "    Hyacinth: área: 253 mm^2, diámetro: 20.66 mm \n"
        "    Canna: área: 202 mm^2, diámetro: 18.38 mm \n"
        "    Laurel: área: 135 mm^2, diámetro: 15.05 mm \n"
        "    Oxlip: área: 107 mm^2, diámetro: 13.25 mm \n"
        "    Poppy: área: 53.7 mm^2, diámetro: 9.36 mm \n"
        "\n"
        "- ACAR: Aluminum Conductor Alloy Reinforced:\n"
        "    ACAR1: área: 1013 mm^2, diámetro: 41.41 mm \n"
        "    ACAR2: área: 963 mm^2, diámetro: 40.37 mm \n"
        "    ACAR3: área: 887 mm^2, diámetro: 38.72 mm \n"
        "    ACAR4: área: 811 mm^2, diámetro: 37.04 mm \n"
        "    ACAR5: área: 760 mm^2, diámetro: 35.85 mm \n"
        "    ACAR6: área: 659 mm^2, diámetro: 33.37 mm \n"
        "    ACAR7: área: 557 mm^2, diámetro: 30.65 mm \n"
        "    ACAR8: área: 481 mm^2, diámetro: 28.48 mm \n"
        "    ACAR9: área: 380 mm^2, diámetro: 25.32 mm \n"
        "    ACAR10: área: 279 mm^2, diámetro: 21.67 mm \n"
        "    ACAR11: área: 177 mm^2, diámetro: 17.24 mm \n"
        "    ACAR12: área: 125 mm^2, diámetro: 14.31 mm \n"
        "    ACAR13: área: 62.5 mm^2, diámetro: 10.11 mm \n"
    ),
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
        "Seleccionar todo: Selecciona a todos los gráficos  disponibles según sean el caso unipolar o bipolar\n"
        "Gráficos disponibles:\n"
        "- Eele: Campo electrostático (V/m)\n"
        "- Vele: Potencial electrostático (V)\n"
        "- Rhof: Densidad de carga espacial total (C/m^3)\n"
        "- Vf: Potencial iónico (V)\n"
        "- SPd: Densidades de carga vista 3D, disponible solamente para la configuración bipolar (C/m^3)\n"
        "- Rhop: Densidad de carga positiva, disponible solamente para configuración bipolar (C/m^3)\n"
        "- Rhon: Densidad de carga negativa, disponible solamente para configuración bipolar (C/m^3)\n"
        "- Vdef: Potencial definitivo (V)\n"
        "- Ef: Campo definitivo (V/m)\n"
        "- J1: Perfil densidad de corriente a nivel de piso (nA/m^2)\n"
        "- E1: Perfil magnitud campo eléctrico total a nivel de piso (V/m)\n"
        "- SPv: Potenciales vista 3D (V)\n"
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
        ("Estado sup cond 1", "1"),
        ("Subconductores", "1"),
        ("Separación (cm)", "0"),
        ("Voltaje 2 (kV)", "-400"),
        ("Área 2 (mm²)", "100"),
        ("Posición en x 2 (m)", "0"),
        ("Posición en y 2 (m)", "5"),
        ("Estado sup cond 2", "1"),
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
        ("Interior conductor", "no")
    ]
}

# Función unificada para actualizar opciones
def actualizar_nombres(tipo_var, nombre_var, menu):
    tipo = tipo_var.get()
    nombres = [x[0] for x in opciones_subcon.get(tipo, [])]
    nombre_var.set("")  # Resetear la selección
    menu["values"] = nombres


entradas = {}
# Configuración de las opciones por tipo de conductor

opciones_subcon = {
    "ACSR": [
        ("Kiwi", 1170, 44.12),
        ("Chukar", 976.7, 40.69),
        ("Falcon", 908.7, 39.24),
        ("Lapwing", 859.8, 38.15),
        ("Plover", 818.7, 37.21),
        ("Pheasant", 726.8, 35.10),
        ("Bunting", 647.7, 33.07),
        ("Curlew", 593.6, 31.65),
        ("Rail", 517.3, 29.59),
        ("Condor", 454.5, 27.76),
        ("Crow", 408.5, 26.28),
        ("Gull", 361.0, 25.38),
        ("Kingbird", 341.0, 23.88),
        ("Osprey", 298.2, 22.33),
        ("Pelican", 255.1, 20.68),
        ("Oriole", 210.3, 18.82),
        ("Ostrich", 176.7, 17.27),
        ("Penguin (4/0)", 125.1, 14.30),
        ("Pigeon (3/0)", 99.2, 12.75),
        ("Raven (1/0)", 62.4, 10.11),
        ("Sparrow", 39.2, 8.03),
        ("Turkey", 15.5, 5.03)
    ],
    "AAC": [
        ("Cowslip", 1013, 41.41),
        ("Jessamine", 887, 38.72),
        ("Carnation", 725, 35.02),
        ("Narcissus", 645, 32.94),
        ("Marigold", 564, 30.88),
        ("Goldenrod", 483, 28.60),
        ("Lilac", 403, 26.11),
        ("Flag", 355, 24.48),
        ("Hyacinth", 253, 20.66),
        ("Canna", 202, 18.38),
        ("Laurel", 135, 15.05),
        ("Oxlip", 107, 13.25),
        ("Poppy", 53.7, 9.36)
    ],
    "ACAR": [
        ("ACAR1", 1013, 41.41),
        ("ACAR2", 963, 40.37),
        ("ACAR3", 887, 38.72),
        ("ACAR4", 811, 37.04),
        ("ACAR5", 760, 35.85),
        ("ACAR6", 659, 33.37),
        ("ACAR7", 557, 30.65),
        ("ACAR8", 481, 28.48),
        ("ACAR9", 380, 25.32),
        ("ACAR10", 279, 21.67),
        ("ACAR11", 177, 17.24),
        ("ACAR12", 125, 14.31),
        ("ACAR13", 62.5, 10.11)
    ]
}


#opciones_areas = ["50", "70", "95", "120", "150", "185", "240", "300", "400", "500", "630"]  # Áreas en mm²
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
            # Extraer las categorías del diccionario
            categorias = list(opciones_subcon.keys())
            # Código principal
            if "Área" in nombre:
                tipo_conductor = tk.StringVar()
                menu_tipo = ttk.Combobox(marco, textvariable=tipo_conductor, state="readonly", width=7)
                menu_tipo["values"] = list(opciones_subcon.keys())
                menu_tipo.grid(row=j, column=1, padx=5, pady=2)
                menu_tipo.bind("<<ComboboxSelected>>", lambda e: actualizar_nombres(tipo_conductor, nombre_conductor, menu_nombre))

                nombre_conductor = tk.StringVar()
                menu_nombre = ttk.Combobox(marco,textvariable=nombre_conductor, state="readonly", width=7)
                menu_nombre.grid(row=j, column=2, padx=5, pady=2)
                entradas[nombre] = []
            # Crear Entry para otros campos
            else:
                entrada = ttk.Entry(marco, width=9)
                entrada.insert(0, valor_defecto)
                entrada.grid(row=j, column=1, padx=5, pady=2)
                entradas[nombre] = entrada
            # Botón de ayuda
            boton_ayuda = ttk.Button(marco, text="?", width=2, command=lambda p=nombre: mostrar_explicacion(p))
            if j==1:
                boton_ayuda.grid(row=j, column=3, padx=2, pady=1)
            else:
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
            ttk.Label(marco, text=nombre).grid(row=j, column=4, sticky="w", padx=5, pady=2)
            categorias2 = list(opciones_subcon.keys())
            # Crear Combobox para "Área"
            if "Área" in nombre:
                ################
                tipo_conductor2 = tk.StringVar()
                menu_tipo2 = ttk.Combobox(marco, textvariable=tipo_conductor2, state="disabled", width=7)
                menu_tipo2["values"] = list(opciones_subcon.keys())
                menu_tipo2.grid(row=j, column=5, padx=5, pady=2)
                menu_tipo2.bind("<<ComboboxSelected>>", lambda e: actualizar_nombres(tipo_conductor2, nombre_conductor2, menu_nombre2))

                nombre_conductor2 = tk.StringVar()
                menu_nombre2 = ttk.Combobox(marco, textvariable=nombre_conductor2, state="disabled", width=7)
                menu_nombre2.grid(row=j, column=6, padx=5, pady=2)
                entradas[nombre] = []
                ################
            # Crear Entry para otros campos
            else:
                entrada = ttk.Entry(marco, width=9)
                entrada.insert(0, valor_defecto)
                entrada.grid(row=j, column=5, padx=5, pady=2)
                entradas[nombre] = entrada
            # Botón de ayuda
            boton_ayuda = ttk.Button(marco, text="?", width=2,command=lambda p=nombre: mostrar_explicacion(p))
            if j==1:
                boton_ayuda.grid(row=j, column=7, padx=2, pady=1)
            else:
                boton_ayuda.grid(row=j, column=6, padx=2, pady=1)
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
        entrada = ttk.Combobox(marco_ambientales, values=opciones_condct, state="disabled")
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
    entrada.config(state="disabled")
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
entradas["Estado sup cond 1"].config(validate="key", validatecommand=(validacion_uno, "%P", "Estado sup cond 1"))
entradas["Estado sup cond 2"].config(validate="key", validatecommand=(validacion_uno, "%P", "Estado sup cond 2"))

entrada_cantidad.bind("<KeyRelease>", verificar_separacion)  # Asociar el evento de escritura
entrada_separacion.config(state="disabled")  # Desactivada inicialmente

# Validación coeficiente rugosidad piso para  viento
entrada_modo = entradas["Modo (str)"]
entrada_rugosidad_terreno = entradas["Rugosidad terreno"]
entrada_modo.bind("<<ComboboxSelected>>", verificar_modo_viento) # Asociar evento de selección al Combobox
entrada_rugosidad_terreno.config(state='disabled') # Desactivada inicialmente


# Crear un estilo para los botones normales y seleccionados
style = ttk.Style()
# Estilo para el botón normal (sin sombra)
style.configure("BotonNormal.TButton",
                relief="flat",  # Sin relieve
                padding=5,
                borderwidth=1,  # Borde fino
                highlightthickness=0,  # Sin contorno
                width=15,  # Establecer un tamaño fijo de ancho
                height=2,  # Establecer un tamaño fijo de altura
                )

# Estilo para el botón seleccionado (con sombra)
style.configure("BotonSeleccionado.TButton",
                relief="flat",  # Sin relieve
                padding=5,
                borderwidth=3,  # Borde más grueso para simular la sombra
                highlightthickness=10,  # Contorno de color gris simulando sombra
                highlightbackground="gray",  # Color de la "sombra"
                highlightcolor="gray",  # Color de la "sombra" cuando está enfocado
                width=15,  # Mantener el tamaño fijo
                height=2,  # Mantener el tamaño fijo
                )

estado_botones = {}

def alternar_estado_boton(boton):
    """
    Alterna el estado de un botón entre seleccionado y no seleccionado.
    
    Args:
        boton (ttk.Button): El botón al que se aplicará el cambio.
    """
    if estado_botones.get(boton, False):
        boton.configure(style="BotonNormal.TButton")  # Cambia al estilo normal
        estado_botones[boton] = False
    else:
        boton.configure(style="BotonSeleccionado.TButton")  # Cambia al estilo seleccionado
        estado_botones[boton] = True

# Crear los botones
marco_botones = ttk.Frame(ventana, padding=10)
marco_botones.grid(row=len(secciones) + 1, column=0, columnspan=3, sticky="w", padx=10, pady=5)

# Definir ancho uniforme
ancho_botones = 7  # Ajusta este valor según tus necesidades

# Crear los botones y definir su ancho
boton_guardar = ttk.Button(marco_botones, text="Guardar", command=lambda: [guardar_datos(), alternar_estado_boton(boton_guardar)], width=ancho_botones)
boton_guardar.grid(row=0, column=0, padx=5)
boton_limpiar = ttk.Button(marco_botones, text="Limpiar", command=limpiar_campos, width=ancho_botones)
boton_limpiar.grid(row=0, column=1, padx=5)
boton_graficar = ttk.Button(marco_botones, text="Graficar", command=validar_campos_y_graficar, width=ancho_botones)
boton_graficar.grid(row=0, column=2, padx=5)
boton_ejecutar = ttk.Button(marco_botones, text="Ejecutar", command=validar_campos_y_ejecutar, width=ancho_botones)
boton_ejecutar.grid(row=0, column=3, padx=5)
boton_salir = ttk.Button(marco_botones, text="Salir", command=cerrar_programa, width=ancho_botones)
boton_salir.grid(row=0, column=4, padx=5)

# Protocolo de cierre
ventana.protocol("WM_DELETE_WINDOW", cerrar_programa)



# Crear el botón Bipolar inmediatamente a la derecha
boton_bipolar_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_bipolar = ttk.Checkbutton(marco_botones, text="Bipolar", variable=boton_bipolar_var, command=activar_bipolar)
boton_bipolar.grid(row=0, column=5, padx=5, pady=10)  # Ubicado justo después de los botones principales

# Añadir un nuevo botón "Auto red" inmediatamente a la derecha del Bipolar
boton_auto_red_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_auto_red = ttk.Checkbutton(marco_botones, text="Auto red", variable=boton_auto_red_var, command=activar_auto_red)
boton_auto_red.grid(row=0, column=6, padx=5, pady=10)  # Ubicado justo al lado del Bipolar

# Añadir un nuevo botón "Desarrollador" inmediatamente a la derecha de autored
boton_desarrollador_var = tk.BooleanVar()  # Variable para controlar el estado del botón
boton_desarrollador = ttk.Checkbutton(marco_botones, text="Desarrollador", variable=boton_desarrollador_var, command=activar_desarrollador)
boton_desarrollador.grid(row=0, column=7, padx=5, pady=10)  # Ubicado justo al lado del Bipolar


# Botones de gráficos
marco_contenedor = tk.Frame(ventana)
marco_contenedor.grid(row=len(secciones) + 3, column=0, columnspan=2, padx=10, pady=5)

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



# Crear el botón "Seleccionar todos"
boton_seleccionar_todos = ttk.Button(marco_graficos, text="Seleccionar Todos", command=lambda: [seleccionar_todos(), alternar_estado_boton(boton_seleccionar_todos)])

# Ubicar el botón en la 3ra columna (índice 2)
boton_seleccionar_todos.grid(row=0, column=2, padx=10, pady=5, sticky="ew")

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










