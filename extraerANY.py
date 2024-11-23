import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

def extraer_datos_txt(ruta_archivo):
    """
    Extrae datos de un archivo .txt con un formato específico.
    """
    with open(ruta_archivo, "r") as file:
        lines = file.readlines()

    encabezado_clave = "MTR    KV/M    KV/M  nA/M**2  PER CC"
    start_index = None

    for i, line in enumerate(lines):
        if encabezado_clave in line:
            start_index = i + 1  # La siguiente línea después del encabezado
            break

    if start_index is None:
        raise ValueError("No se encontró la fila con los encabezados en el archivo.")

    data_lines = lines[start_index:]

    colspecs = [
        (0, 8),   # MTR
        (9, 17),  # NOMINAL_FIELD (KV/M)
        (18, 26), # ENHANCED_FIELD (KV/M)
        (27, 35), # GROUND_CURRENT (nA/M**2)
        (36, 44), # GROUND_CHARGES (PER CC)
    ]

    df = pd.read_fwf(pd.io.common.StringIO("".join(data_lines)), colspecs=colspecs, header=None)
    df.columns = ["MTR", "NOMINAL_FIELD", "ENHANCED_FIELD", "GROUND_CURRENT", "GROUND_CHARGES"]

    # Convertir columnas a numéricas
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.dropna()

def extraer_datos_xlsx(ruta_archivo):
    """
    Extrae datos de un archivo .xlsx con un formato específico.
    """
    encabezado_clave = [
        "Lateral Distance x[m]", "|E| [kV/m]", "Q [nC/m^3]",
        "J [nA/m^2]", "Ions [k/m^3]", "|E| charge-free"
    ]

    df = pd.read_excel(ruta_archivo, engine='openpyxl', header=None)
    inicio_datos = None

    for i in range(len(df)):
        fila = df.iloc[i].astype(str).str.strip().tolist()
        if all(any(clave in celda for celda in fila) for clave in encabezado_clave):
            inicio_datos = i + 1  # La fila debajo del encabezado
            break

    if inicio_datos is None:
        raise ValueError("Encabezado no encontrado en el archivo Excel.")

    df_datos = pd.read_excel(ruta_archivo, skiprows=inicio_datos, engine='openpyxl', header=None)
    df_datos.columns = encabezado_clave
    return df_datos.apply(pd.to_numeric, errors='coerce')

def extraer_columnas_como_arrays(df, columnas):
    """
    Extrae columnas específicas de un DataFrame como arrays de NumPy.
    """
    arrays = {}
    for columna in columnas:
        if columna in df.columns:
            arrays[columna] = df[columna].dropna().to_numpy()
        else:
            print(f"Advertencia: La columna '{columna}' no existe en el DataFrame.")
    return arrays

def graficar_datos(x, y, xlabel, ylabel, titulo, color, label):
    """
    Grafica los datos proporcionados.
    """
    plt.plot(x, y, label=label, color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titulo)
    plt.legend()
    plt.grid(True)

# Procesar y graficar datos
def procesar_archivo(ruta_archivo):
    """
    Procesa un archivo .txt o .xlsx y genera gráficos.
    """
    if ruta_archivo.endswith('.txt'):
        df = extraer_datos_txt(ruta_archivo)
        posicion = df['MTR'].to_numpy()
        campo_nominal = df['NOMINAL_FIELD'].to_numpy()
        campo_enhanced = df['ENHANCED_FIELD'].to_numpy()
        corriente = df['GROUND_CURRENT'].to_numpy()
        print(len(posicion))

        plt.figure()
        graficar_datos(posicion, campo_enhanced, 'Posición (MTR)', 'Campo Eléctrico (kV/m)',
                       'Campo Eléctrico Mejorado', 'blue', 'Enhanced Field')
        graficar_datos(posicion, campo_nominal, 'Posición (MTR)', 'Campo Eléctrico (kV/m)',
                       'Campo Nominal', 'green', 'Nominal Field')

        plt.figure()
        graficar_datos(posicion, corriente, 'Posición (FT)', 'Densidad de Corriente (nA/m²)',
                       'Densidad de Corriente', 'orange', 'Ground Current')

    elif ruta_archivo.endswith('.xlsx'):
        df = extraer_datos_xlsx(ruta_archivo)
        columnas_a_extraer = ['Lateral Distance x[m]', '|E| [kV/m]', 'J [nA/m^2]']
        arrays = extraer_columnas_como_arrays(df, columnas_a_extraer)

        plt.figure()
        graficar_datos(arrays['Lateral Distance x[m]'], arrays['|E| [kV/m]'],
                       'Distancia Lateral (m)', 'Campo Eléctrico (kV/m)',
                       'Campo Eléctrico vs Distancia', 'blue', '|E|')

        plt.figure()
        graficar_datos(arrays['Lateral Distance x[m]'], arrays['J [nA/m^2]'],
                       'Distancia Lateral (m)', 'Densidad de Corriente (nA/m²)',
                       'Densidad de Corriente vs Distancia', 'orange', 'J')

    else:
        raise ValueError("Formato de archivo no soportado. Solo se aceptan .txt y .xlsx.")

    plt.show()

# Rutas de ejemplo
ruta_txt = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\Anypole1\AP2006_1.txt'  # Ruta a un archivo de texto
ruta_xlsx = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\FACE_results\Cas1_idfp3MARUVADA.xlsx'  # Ruta a un archivo de Excel

try:
    procesar_archivo(ruta_txt)
except Exception as e:
    print(f"Error al procesar el archivo .txt: {e}")

try:
    procesar_archivo(ruta_xlsx)
except Exception as e:
    print(f"Error al procesar el archivo .xlsx: {e}")












