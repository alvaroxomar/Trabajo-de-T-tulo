import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

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

'''
def extraer_datos_xlsx(ruta_archivo, es_face= False):
    """
    Extrae datos de un archivo .xlsx basado en los encabezados específicos que se encuentren.
    Ignora encabezados faltantes en el archivo Excel.

    Args:
        ruta_archivo (str): Ruta al archivo Excel.

    Returns:
        pd.DataFrame: DataFrame con los datos extraídos.
    """
    # Encabezados deseados
    encabezados_deseados = ["Lateral Distance x[m]", "|E| [kV/m]", "J [nA/m^2]"]

    # Leer archivo Excel
    if es_face:
        df = pd.read_excel(ruta_archivo, engine='openpyxl', header=1)  # Encabezados en la fila 2 (índice 1)
    else:
        df = pd.read_excel(ruta_archivo, engine='openpyxl', header=None)
    # Normalizar encabezados deseados para evitar errores por espacios o casos
    encabezados_deseados_normalizados = [encabezado.strip().lower() for encabezado in encabezados_deseados]

    # Buscar fila con encabezados que coincidan
    inicio_datos = None
    encabezados_reales = None

    for i in range(len(df)):
        # Normalizar fila leída
        fila = df.iloc[i].astype(str).str.strip().str.lower().tolist()

        # Imprimir fila para depurar
        #print(f"Analizando fila {i}: {fila}")

        # Verificar si al menos uno de los encabezados deseados está en la fila
        if any(encabezado in fila for encabezado in encabezados_deseados_normalizados):
            inicio_datos = i + 1  # Fila donde comienzan los datos
            encabezados_reales = df.iloc[i].tolist()  # Conservar encabezados originales
            break

    if inicio_datos is None:
        raise ValueError("No se encontraron encabezados válidos en el archivo Excel.")

    # Leer datos desde la fila encontrada
    df_datos = pd.read_excel(ruta_archivo, skiprows=inicio_datos, engine='openpyxl', header=None)

    # Mapear las columnas a los encabezados reales encontrados
    df_datos.columns = encabezados_reales

    # Filtrar columnas que coincidan con los encabezados deseados
    columnas_presentes = [col for col in encabezados_deseados if col in df_datos.columns]

    if not columnas_presentes:
        raise ValueError("Ninguna de las columnas deseadas está presente en el archivo Excel.")

    # Extraer columnas presentes y convertir a numérico
    df_filtrado = df_datos[columnas_presentes].apply(pd.to_numeric, errors='coerce')

    return df_filtrado
'''
def extraer_datos_xlsx(ruta_archivo, es_face=False):
    encabezados_deseados = ["Lateral Distance x[m]", "|E| [kV/m]", "J [nA/m^2]"]
    encabezados_deseados_normalizados = [
        encabezado.strip().lower().replace(" ", "").replace("[", "").replace("]", "")
        for encabezado in encabezados_deseados
    ]

    # Leer el archivo
    df = pd.read_excel(ruta_archivo, engine="openpyxl", header=None)

    inicio_datos = None
    encabezados_reales = None

    for i in range(len(df)):
        fila = df.iloc[i].astype(str).str.strip().str.lower().str.replace(" ", "").str.replace("[", "").str.replace("]", "").tolist()
        print(f"La fila analizada es: {fila}")

        if any(encabezado in fila for encabezado in encabezados_deseados_normalizados):
            inicio_datos = i + 1
            encabezados_reales = df.iloc[i].tolist()  # Guardar encabezados originales
            break

    if inicio_datos is None:
        raise ValueError("No se encontraron encabezados válidos en el archivo Excel.")

    # Leer datos desde la fila encontrada
    df_datos = pd.read_excel(ruta_archivo, skiprows=inicio_datos, engine="openpyxl", header=None)

    # Normalizar encabezados reales
    encabezados_reales_normalizados = [
        str(encabezado).strip().lower().replace(" ", "").replace("[", "").replace("]", "") for encabezado in encabezados_reales
    ]
    df_datos.columns = encabezados_reales_normalizados

    # Filtrar columnas
    columnas_presentes = [
        col for col in encabezados_deseados_normalizados if col in df_datos.columns
    ]

    if not columnas_presentes:
        raise ValueError("Ninguna de las columnas deseadas está presente en el archivo Excel.")

    print(f"Columnas presentes: {columnas_presentes}")
    print(f'datos presentes:  {df_datos[columnas_presentes]}')
    # Extraer columnas y convertir a numérico
    df_filtrado = df_datos[columnas_presentes].apply(pd.to_numeric, errors="coerce")

    # Crear un diccionario que mapea nombres normalizados a los nombres originales
    mapa_columnas = {
        encabezado.strip().lower().replace(" ", "").replace("[", "").replace("]", ""): encabezado
        for encabezado in encabezados_deseados
    }

    # Renombrar las columnas de df_filtrado con los nombres originales
    df_filtrado.columns = [mapa_columnas[col] for col in df_filtrado.columns]

    print(f"Datos filtrados con nombres originales: {df_filtrado}")

    return df_filtrado




def extraer_datos_csv(ruta_archivo):
    """
    Extrae datos de un archivo .csv con un formato específico.
    """
    try:
        df = pd.read_csv(ruta_archivo)
        return df.apply(pd.to_numeric, errors='coerce')
    except Exception as e:
        raise ValueError(f"Error al leer el archivo CSV: {e}")

def extraer_columnas_como_arrays(df, columnas):
    """
    Extrae columnas específicas de un DataFrame como arrays de NumPy.
    """
    df.columns = df.columns.str.strip()  # Elimina espacios en blanco al inicio y al final
    arrays = {}
    for columna in columnas:
        #print(f'columa identificada es: {columna}')
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
    plt.xlabel(xlabel,fontsize=13)
    plt.ylabel(ylabel,fontsize=13)
    plt.title(titulo,fontsize=13)
    plt.legend()
    plt.grid(True)

# Procesar y graficar datos
def procesar_archivo(ruta_archivo):
    """
    Procesa un archivo .txt, .xlsx, .csv o .out y genera gráficos.
    """
    global col
    if ruta_archivo.endswith('.out'):
        nueva_ruta = ruta_archivo + '.txt'
        os.rename(ruta_archivo, nueva_ruta)
        try:
            procesar_archivo(nueva_ruta)
        finally:
            os.rename(nueva_ruta, ruta_archivo)

    elif ruta_archivo.endswith('.txt'):
        df = extraer_datos_txt(ruta_archivo)
        posicion = df['MTR'].to_numpy()
        campo_nominal = df['NOMINAL_FIELD'].to_numpy()
        campo_enhanced = df['ENHANCED_FIELD'].to_numpy()
        corriente = df['GROUND_CURRENT'].to_numpy()
        print(len(posicion))

        plt.figure(1)
        graficar_datos(posicion, campo_enhanced, 'Posición (MTR)', 'Campo Eléctrico (kV/m)',
                       'Campo Eléctrico Mejorado', 'blue', '|E| ANYPOLE')
        graficar_datos(posicion, campo_nominal, 'Posición (MTR)', 'Campo Eléctrico (kV/m)',
                       'Campo Nominal', 'red', 'nominal |E| ANYPOLE')
        #graficar_datos(x[:200], Ei[:200]/1000, 'Posición (MTR)', 'Campo Eléctrico (kV/m)',
        #               'Campo programa', 'red', 'Nominal Field PROGRAMA')

        plt.figure(2)
        graficar_datos(posicion, corriente, 'Posición (FT)', 'Densidad de Corriente (nA/m²)',
                       'Densidad de Corriente', 'blue', 'J ANYPOLE')

    elif ruta_archivo.endswith('.xlsx'):
        df = extraer_datos_xlsx(ruta_archivo)
        encabezados = df.columns.tolist()
        print("Encabezados como lista:")
        print(encabezados)
        columnas_a_extraer = ['Lateral Distance x[m]', '|E| [kV/m]', 'J [nA/m^2]']
        #columnas_a_extraer_normalizados = [
        #encabezado.strip().lower().replace(" ", "").replace("[", "").replace("]", "")
        #for encabezado in columnas_a_extraer]
        arrays = extraer_columnas_como_arrays(df, columnas_a_extraer)
        if col:
            plt.figure(1)
            graficar_datos(arrays['Lateral Distance x[m]'], arrays['|E| [kV/m]'],
                        'Distancia Lateral (m)', 'Campo Eléctrico (kV/m)',
                        'Campo Eléctrico vs Distancia', 'orange', '|E| programa')

            plt.figure(2)
            graficar_datos(arrays['Lateral Distance x[m]'], arrays['J [nA/m^2]'],
                        'Distancia Lateral (m)', 'Densidad de Corriente (nA/m²)',
                        'Densidad de Corriente vs Distancia', 'orange', 'J programa')
            col = False
        else:
            plt.figure(1)
            graficar_datos(arrays['Lateral Distance x[m]'], arrays['|E| [kV/m]'],
                        'Distancia Lateral (m)', 'Campo Eléctrico (kV/m)',
                        'Campo Eléctrico vs Distancia', 'purple', '|E| FACE')

            plt.figure(2)
            graficar_datos(arrays['Lateral Distance x[m]'], arrays['J [nA/m^2]'],
                        'Distancia Lateral (m)', 'Densidad de Corriente (nA/m²)',
                        'Densidad de Corriente vs Distancia', 'purple', 'J FACE')


    elif ruta_archivo.endswith('.csv'):
        df = extraer_datos_csv(ruta_archivo)
        columnas_a_extraer = ['Lateral Distance x[m]', '|E| [kV/m]', 'J [nA/m^2]']
        arrays = extraer_columnas_como_arrays(df, columnas_a_extraer)

        plt.figure(1)
        graficar_datos(arrays['Lateral Distance x[m]'], arrays['|E| [kV/m]'],
                       'Distancia Lateral (m)', 'Campo Eléctrico (kV/m)',
                       'Campo Eléctrico vs Distancia', 'green', '|E| CSV')

        plt.figure(2)
        graficar_datos(arrays['Lateral Distance x[m]'], arrays['J [nA/m^2]'],
                       'Distancia Lateral (m)', 'Densidad de Corriente (nA/m²)',
                       'Densidad de Corriente vs Distancia', 'green', 'J CSV')

    else:
        raise ValueError("Formato de archivo no soportado. Solo se aceptan .txt, .xlsx, .csv y .out.")

    plt.show()

# Rutas de ejemplo
ruta_txt = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\Anypole1\AP2006_4.out'  # Ruta a un archivo de texto
ruta_xlsx1 = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\programa resultados\modelo_250.0_10_600_400.xlsx'  # Ruta a un archivo de Excel
ruta_xlsx2 = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\FACE_results\Cas1_idfp3MARUVADA.xlsx'
#ruta_csv = r'C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\FACE_results\Cas1_idfBIPmaruvada_m05.csv'  # Ruta a un archivo CSV

col=True
try:
    procesar_archivo(ruta_txt)
except Exception as e:
    print(f"Error al procesar el archivo .txt: {e}")

try:
    procesar_archivo(ruta_xlsx1)
except Exception as e:
    print(f"Error al procesar el archivo .xlsx: {e}")

try:
    procesar_archivo(ruta_xlsx2)
except Exception as e:
    print(f"Error al procesar el archivo .xlsx: {e}")

try:
    procesar_archivo(ruta_csv)
except Exception as e:
    print(f"Error al procesar el archivo .csv: {e}")

