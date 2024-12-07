import tkinter as tk
from tkinter import messagebox
import json
import subprocess
import matplotlib.pyplot as plt

plt.close('all')
def interpretar_valor(valor):
    """
    Interpreta un valor ingresado como float, int, bool o None.
    """
    if not valor:
        return None  # Si está vacío, retorna None
    if valor.lower() == 'true':
        return True
    if valor.lower() == 'false':
        return False
    try:
        if '.' in valor:
            return float(valor)  # Intenta convertir a float si hay un punto
        return int(valor)  # Intenta convertir a int
    except ValueError:
        raise ValueError(f"El valor ingresado '{valor}' no es válido. Use números, 'True', 'False' o deje vacío.")
def ejecutar_script():
    try:
        # Leer los valores de los parámetros
        params = {}

        # Parámetros con valores por defecto
        params["m"] = float(entry_m.get()) if entry_m.get() else 1
        params["Pr"] = float(entry_Pr.get()) if entry_Pr.get() else 101
        params["Tr"] = float(entry_Tr.get()) if entry_Tr.get() else 303
        params["l"] = float(entry_l.get()) if entry_l.get() else 1
        params["wndx"] = float(entry_wndx.get()) if entry_wndx.get() else 0
        params["wndy"] = float(entry_wndy.get()) if entry_wndy.get() else 0
        params["modo"] = entry_modo.get() if entry_modo.get() else 'uniforme'
        params["coef_drag"] = float(entry_coef_drag.get()) if entry_coef_drag.get() else 0.3
        params["max_iter_rho"] = int(entry_max_iter_rho.get()) if entry_max_iter_rho.get() else 400
        params["max_iter"] = int(entry_max_iter.get()) if entry_max_iter.get() else 250
        params["it_global"] = int(entry_it_global.get()) if entry_it_global.get() else 10
        params["TolDev"] = float(entry_TolDev.get()) if entry_TolDev.get() else 1e-2
        params["visualizacion"] = int(entry_visualizacion.get()) if entry_visualizacion.get() else 15
        params["in_condct"] = entry_in_condct.get() if entry_in_condct.get() else 'si'
        params["copiado"] = entry_copiado.get() if entry_copiado.get() else 'no'
        params["histl"] = bool(entry_histl.get()) if entry_histl.get() else False
        params["mostrar"] = bool(entry_mostrar.get()) if entry_mostrar.get() else False
        # Parámetros que admiten float, int, bool o None
        params["Sx"] = interpretar_valor(entry_Sx.get())
        params["Sy"] = interpretar_valor(entry_Sy.get())
        params["nodosx"] = interpretar_valor(entry_nodosx.get())
        params["nodosy"] = interpretar_valor(entry_nodosy.get())
        # Parámetros obligatorios (sin valor por defecto)
        required_params = [
            ("R", entry_R),
            ("Vol", entry_Vol),
            ("x_coor", entry_x_coor),
            ("y_coor", entry_y_coor),
        ]

        # Verificar si los parámetros obligatorios están vacíos
        for param_name, entry_widget in required_params:
            if not entry_widget.get():
                messagebox.showerror("Error", f"El parámetro {param_name} no puede estar vacío.")
                return  # Detener ejecución si hay campos vacíos
            params[param_name] = float(entry_widget.get())

        # Convertir los parámetros a JSON
        params_json = json.dumps(params)

        # Ejecutar el script con los parámetros
        subprocess.run(["python", "com55.py", params_json], check=True)
        messagebox.showinfo("Éxito", "El archivo com55.py se ejecutó correctamente.")
    except Exception as e:
        messagebox.showerror("Error", f"Error al ejecutar el script: {e}")

# Configuración de la ventana principal
root = tk.Tk()
root.title("Panel de Control")

# Lista de parámetros y etiquetas asociadas
parametros = [
    ("R", "Radio equivalente (m)", tk.Entry(root)),
    ("Vol", "Voltaje de la línea (V)", tk.Entry(root)),
    ("m", "Factor de rugosidad", tk.Entry(root, textvariable=tk.StringVar(value="1"))),  # Valor por defecto
    ("Pr", "Presión del aire (kPa)", tk.Entry(root, textvariable=tk.StringVar(value="101"))),  # Valor por defecto
    ("Tr", "Temperatura del sitio (K)", tk.Entry(root, textvariable=tk.StringVar(value="303"))),  # Valor por defecto
    ("x_coor", "Coordenada x (m)", tk.Entry(root)),
    ("y_coor", "Coordenada y (m)", tk.Entry(root)),
    ("l", "Altura de interés (m)", tk.Entry(root, textvariable=tk.StringVar(value="1"))),  # Valor por defecto
    ("Sx", "Longitud del plano de tierra (m)", tk.Entry(root)),
    ("Sy", "Altura del área (m)", tk.Entry(root)),
    ("nodosx", "Número de nodos en x", tk.Entry(root)),
    ("nodosy", "Número de nodos en y", tk.Entry(root)),
    ("wndx", "Velocidad del viento x (m/s)", tk.Entry(root, textvariable=tk.StringVar(value="0"))),  # Valor por defecto
    ("wndy", "Velocidad del viento y (m/s)", tk.Entry(root, textvariable=tk.StringVar(value="0"))),  # Valor por defecto
    ("modo", "Modo de simulación", tk.Entry(root, textvariable=tk.StringVar(value="uniforme"))),  # Valor por defecto
    ("coef_drag", "Coeficiente de arrastre", tk.Entry(root, textvariable=tk.StringVar(value="0.3"))),  # Valor por defecto
    ("max_iter_rho", "Máx iter densidad carga", tk.Entry(root, textvariable=tk.StringVar(value="400"))),  # Valor por defecto
    ("max_iter", "Máx iter potencial", tk.Entry(root, textvariable=tk.StringVar(value="250"))),  # Valor por defecto
    ("it_global", "Iteraciones globales", tk.Entry(root, textvariable=tk.StringVar(value="10"))),  # Valor por defecto
    ("TolDev", "Tolerancia convergencia", tk.Entry(root, textvariable=tk.StringVar(value="1e-2"))),  # Valor por defecto
    ("visualizacion", "Ventana visualización", tk.Entry(root, textvariable=tk.StringVar(value="15"))),  # Valor por defecto
    ("in_condct", "Condición de conductor", tk.Entry(root, textvariable=tk.StringVar(value="si"))),  # Valor por defecto
    ("copiado", "Condición de copiado", tk.Entry(root, textvariable=tk.StringVar(value="no"))),  # Valor por defecto
    ("histl", "Historial de iteraciones", tk.Entry(root, textvariable=tk.StringVar(value="False"))),  # Valor por defecto
    ("mostrar", "Mostrar datos convergencia", tk.Entry(root, textvariable=tk.StringVar(value="False"))),  # Valor por defecto
]

# Crear widgets de entrada y etiquetas en una cuadrícula de 4 columnas
for idx, (param, label, entry) in enumerate(parametros):
    row, col = divmod(idx, 4)
    tk.Label(root, text=label).grid(row=row, column=col * 2, sticky="w", padx=5, pady=5)
    entry.grid(row=row, column=col * 2 + 1, padx=5, pady=5)
    globals()[f"entry_{param}"] = entry  # Asignar las entradas a variables globales

# Botón de ejecución
btn_ejecutar = tk.Button(root, text="Ejecutar", command=ejecutar_script)
btn_ejecutar.grid(row=(len(parametros) + 3) // 4, column=0, columnspan=8, pady=10)

# Iniciar la aplicación
root.mainloop()











