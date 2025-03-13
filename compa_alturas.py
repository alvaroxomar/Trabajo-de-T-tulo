import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
print("**********************************")
plt.close("all")
def formatear_numero(valor):
    if np.abs(valor) < 1e-2 or np.abs(valor) > 1e3:
        return f"{valor:.2e}"  # Notación científica para valores extremos
    else:
        return f"{valor:.2f}"  # Notación decimal con 4 decimales
# Lista de archivos Excel
nome1 = r"modeloBIP_800.0_8.0_8.0_524_524.xlsx"
nome2 = r"modeloBIP_800.0_10.0_10.0_524_524.xlsx"
nome3 = r"modeloBIP_800.0_12.0_12.0_524_524.xlsx"
nome4 = r"modeloBIP_700.0_10.0_10.0_672_517.xlsx"
nome = r"modeloBIP_650.0_9.0_9.0_607_607.xlsx"
dir1 = f'C:\\Users\\HITES\\Desktop\\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_800.0_3_8.0_8.0_524_524\\{nome1}'
dir2 = f'C:\\Users\\HITES\\Desktop\\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_800.0_3_10.0_10.0_524_524\\{nome2}'
dir3 = f'C:\\Users\\HITES\\Desktop\\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_800.0_3_12.0_12.0_524_524\\{nome3}'
dir4 = f"C:\\Users\\HITES\\Desktop\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_700.0_3_10.0_10.0_672_517\\Resultados\\{nome4}"
dir5 = f"C:\\Users\\HITES\\Desktop\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_700.0_3_10.0_10.0_672_517\\Resultados_1\\{nome4}"
#dir6 = f"C:\\Users\\HITES\\Desktop\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_700.0_3_10.0_10.0_672_517\\Resultados_new\\{nome4}"
dir7 = f"C:\\Users\\HITES\\Desktop\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_700.0_3_10.0_10.0_672_517\\Resultados_new_new\\{nome4}"
dir8 = f"C:\\Users\\HITES\\Desktop\la uwu\\14vo semestre\\Trabajo de título\\programa resultados\\modeloBIP_700.0_3_10.0_10.0_672_517\\Resultados_new_new_new\\{nome4}"

def ajustar_ruta(ruta, nombre_archivo):
    """
    Reemplaza todas las barras invertidas (\) por doble barra invertida (\\) en una ruta
    y concatena con el nombre de un archivo.

    Args:
        ruta (str): Ruta con barras invertidas.
        nombre_archivo (str): Nombre del archivo a concatenar.

    Returns:
        str: Ruta ajustada con el nombre del archivo.
    """
    # Concatenar la ruta con el nombre del archivo
    ruta_completa = os.path.join(ruta, nombre_archivo)
    # Reemplazar las barras invertidas por doble barra invertida
    #return ruta_completa.replace("\\", "\\\\")
    return ruta_completa
ruta = r"C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\programa resultados\modeloBIP_500.0_3_9.0_9.0_733_611\Resultados"
ruta2 = r"C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\programa resultados\modeloBIP_500.0_3_9.0_9.0_733_611\Resultados_new"
ruta3 = r"C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\programa resultados\modeloBIP_500.0_3_9.0_9.0_733_611\Resultados_new_new"
ruta4 = r"C:\Users\HITES\Desktop\la uwu\14vo semestre\Trabajo de título\programa resultados\modeloBIP_500.0_3_9.0_9.0_733_611\Resultados_new_new_new"
nome = r"modeloBIP_500.0_9.0_9.0_733_611.xlsx"
#nomes = [nome, nome2, nome3]
nomes = nome
rutas = [ruta, ruta2, ruta3, ruta4]
archivos_excel = []
for i in range(len(rutas)):
    dir = ajustar_ruta(rutas[i], nomes)
    archivos_excel.append(dir)

print(archivos_excel)

# Inicializar colores y marcadores para distinguir los archivos
colores = ['tab:blue', 'tab:green', 'tab:purple','tab:orange']
colores_pastel = ['#aec6cf', '#ffb347', '#b39eb5', '#77dd77', '#fdfd96']
marcadores = ['o', 's', '^', '*']  # Cambié a solo los marcadores

def extraer_diccionario_hoja(archivo, hoja):
    """Extrae un diccionario de una hoja de un archivo Excel."""
    try:
        df = pd.read_excel(archivo, sheet_name=hoja)
        return df.to_dict(orient='list')
    except Exception as e:
        print(f"Error al leer la hoja {hoja} en el archivo {archivo}: {e}")
        return None

# Crear figura para Campo Eléctrico
plt.figure(figsize=(8, 6))

plt.xlabel("Lateral Distance x [m]", fontsize=13)
plt.ylabel("Ey [kV/m]", fontsize=13)
plt.grid(True)

# Procesar cada archivo y agregar al gráfico de Campo Eléctrico
supe = []
for idx, archivo in enumerate(archivos_excel):
    hojas = pd.read_excel(archivo, sheet_name=None)
    dic_conductores = extraer_diccionario_hoja(archivo, "Conductores")
    dic_ambientales = extraer_diccionario_hoja(archivo, "Ambientales")
    pos_y = dic_conductores.get("Posición en y (m)", ["N/A"])[0] if dic_conductores else "N/A"
    Esup_izq = dic_conductores.get("Gradiente superficial crítico izq (kV/cm)", ["N/A"])[0] if dic_conductores else "N/A"
    Esup_der = dic_conductores.get("Gradiente superficial crítico der (kV/cm)", ["N/A"])[0] if dic_conductores else "N/A"
    supe.append(Esup_der)
    Wx = dic_ambientales.get("Viento x (m/s)", ["N/A"])[0] if dic_ambientales else "N/A"
    Wy = dic_ambientales.get("Viento y (m/s)", ["N/A"])[0] if dic_ambientales else "N/A"
    alpha = dic_ambientales.get("Rugosidad terreno", ["N/A"])[0] if dic_ambientales else "N/A"
    if idx==0:
        plt.title(f"Campo eléctrico vs Posición, "+r"$E_{on}$"f" = {formatear_numero(supe[0])}", fontsize=15)
    if "Campo y Corriente" in hojas:
        df = hojas["Campo y Corriente"]
        if {"Lateral Distance x[m]", "E_y [kV/m]"}.issubset(df.columns):
            distancia = df["Lateral Distance x[m]"].dropna()
            campo_electrico = df["E_y [kV/m]"].dropna()
            Eprom = df["|E| [kV/m]"].dropna()
            Eav = np.mean(Eprom)
            #plt.plot(distancia[15:-15], campo_electrico[15:-15], label=f" altura = {pos_y} m"+r" , $|E|_{ave}$="f"{formatear_numero(np.mean(Eprom[15:-15]))} kV/m,"+r" $E_{on}$ = "+f"{formatear_numero(supe[idx])} kV/cm",
            #         color=colores[idx], marker=marcadores[idx], linestyle='-', linewidth=1, markersize=5, markevery=10)
            plt.plot(distancia[15:-15], campo_electrico[15:-15], label=f"Wx = {Wx} m/s"+r" , $|E|_{ave}$="f"{formatear_numero(np.mean(Eprom[15:-15]))} kV/m",
                     color=colores[idx], marker=marcadores[idx], linestyle='-', linewidth=1, markersize=5, markevery=10)
            '''
            if idx == len(archivos_excel)-1:
                df2 = pd.read_excel(archivo, sheet_name='Conductores')
                columnas_a_extraer_2 = ['Posición en x (m)', 'Posición en y (m)', 'Posición en x 2(m)', 'Posición en y 2(m)']
                df2 = df2[columnas_a_extraer_2]
                df3 = pd.read_excel(archivo, sheet_name='Dim y discr')
                columnas_a_extraer_3 = ['Altura (m)']
                ymax_ax = float(df3[columnas_a_extraer_3].iloc[0])
                print(f"ymax es {ymax_ax}")
                if any(x < 0 for x in campo_electrico):
                    escala = [-ymax_ax, ymax_ax]
                else:
                    escala = [0, ymax_ax]
                x1 = df2['Posición en x (m)'].iloc[0]
                y1 = df2['Posición en y (m)'].iloc[0]
                x2 = df2['Posición en x 2(m)'].iloc[0]
                y2 = df2['Posición en y 2(m)'].iloc[0]
                plt.legend(loc="lower left")
                ax1 = plt.gca()
                ax2 = ax1.twinx()
                ax2.scatter(x1, y1, color='purple', marker='o', label=f"Polos altura {pos_y} m")  # Agregar el punto
                ax2.scatter(x2, y2, color='purple', marker='o')  # Agregar el punto
                ax2.set_ylabel('Altura (m)', fontsize=11, color='purple')
                ax2.tick_params(axis='y', labelcolor='purple')  # Ajustar el color de las etiquetas del eje secundario
                ax2.set_ylim(escala[0], escala[1])  # Definir el rango de 0 a 20
                ax2.legend(loc="upper right")
            '''
            
        else:
            print(f"Las columnas esperadas no se encontraron en el archivo: {archivo}")
    else:
        print(f"La hoja 'Campo y Corriente' no existe en el archivo: {archivo}")

# Agregar límite, leyenda y mostrar el gráfico
print(f"los gradientes sup son {supe}")

limite = np.ones(len(distancia)) * 25
plt.plot(distancia, limite, color='red', label="Límite CIGRE 25 kV/m")
plt.legend(loc="lower left", fontsize=10)
df2 = pd.read_excel(archivo, sheet_name='Conductores')
columnas_a_extraer_2 = ['Posición en x (m)', 'Posición en y (m)', 'Posición en x 2(m)', 'Posición en y 2(m)']
df2 = df2[columnas_a_extraer_2]
df3 = pd.read_excel(archivo, sheet_name='Dim y discr')
columnas_a_extraer_3 = ['Altura (m)']
ymax_ax = float(df3[columnas_a_extraer_3].iloc[0])
print(f"ymax es {ymax_ax}")
if any(x < 0 for x in campo_electrico):
    escala = [-ymax_ax, ymax_ax]
else:
    escala = [0, ymax_ax]
x1 = df2['Posición en x (m)'].iloc[0]
y1 = df2['Posición en y (m)'].iloc[0]
x2 = df2['Posición en x 2(m)'].iloc[0]
y2 = df2['Posición en y 2(m)'].iloc[0]
plt.legend(loc="lower left")
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.scatter(x1, y1, color='purple', marker='o', label=f"Polos altura {pos_y} m")  # Agregar el punto
ax2.scatter(x2, y2, color='purple', marker='o')  # Agregar el punto
ax2.set_ylabel('Altura (m)', fontsize=11, color='purple')
ax2.tick_params(axis='y', labelcolor='purple')  # Ajustar el color de las etiquetas del eje secundario
ax2.set_ylim(escala[0], escala[1])  # Definir el rango de 0 a 20
ax2.legend(loc="upper right")
plt.tight_layout()
plt.show()

# Crear figura para Densidad de Corriente
plt.figure(figsize=(8, 6))
plt.xlabel("Lateral Distance x [m]", fontsize=13)
plt.ylabel(f"J [nA/$m^2$]", fontsize=13)
plt.grid(True)

# Procesar cada archivo y agregar al gráfico de Densidad de Corriente
for idx, archivo in enumerate(archivos_excel):
    hojas = pd.read_excel(archivo, sheet_name=None)
    dic_conductores = extraer_diccionario_hoja(archivo, "Conductores")
    dic_ambientales = extraer_diccionario_hoja(archivo, "Ambientales")
    pos_y = dic_conductores.get("Posición en y (m)", ["N/A"])[0] if dic_conductores else "N/A"
    Wx = dic_ambientales.get("Viento x (m/s)", ["N/A"])[0] if dic_ambientales else "N/A"
    #Wy = dic_ambientales.get("Viento y (m/s)", ["N/A"])[0] if dic_ambientales else "N/A"
    #alpha = dic_ambientales.get("Rugosidad terreno", ["N/A"])[0] if dic_ambientales else "N/A"
    if idx==0:
        plt.title(f"Densidad de Corriente vs Posición, "+r"$E_{on}$"f" = {formatear_numero(supe[0])}", fontsize=15)
    if "Campo y Corriente" in hojas:
        df = hojas["Campo y Corriente"]
        if {"Lateral Distance x[m]", "J [nA/m^2]"}.issubset(df.columns):
            distancia = df["Lateral Distance x[m]"].dropna()
            corriente = df["J [nA/m^2]"].dropna()
            Jx = df["Jx [nA/m^2]"].dropna()
            Jy = df["Jy [nA/m^2]"].dropna()
            Jav = np.mean(corriente)
            
            #plt.plot(distancia[15:-15], -Jy[15:-15], label=f"Wx = {Wx} m/s,"
            #         +r" $J_{y,ave}$="f"{formatear_numero(np.mean(-Jy[15:-15]))}"+r" $nA/m^2$",
            #         color=colores[idx], marker=marcadores[idx], linestyle='-', linewidth=1, markersize=5, markevery=10)
            plt.plot(distancia[15:-15], corriente[15:-15], label=f"Wx = {Wx} m/s,"
                     +r" $|J|_{ave}$="f"{formatear_numero(np.mean(corriente[15:-15]))}"+r" $nA/m^2$",
                     color=colores[idx], marker=marcadores[idx], linestyle='-', linewidth=1, markersize=5, markevery=10)
            #plt.plot(distancia, Jx, label=f"Componente x de $J$ nA/$m^2$",
            #         color=colores_pastel[idx], marker='*', linestyle='-', linewidth=1, markersize=5)
            #plt.plot(distancia, Jy, label=f"Componente y de $J$ nA/$m^2$",
            #         color=colores_pastel[idx+1], marker='*', linestyle='-', linewidth=1, markersize=5)
        else:
            print(f"Las columnas esperadas no se encontraron en el archivo: {archivo}")
    else:
        print(f"La hoja 'Campo y Corriente' no existe en el archivo: {archivo}")

# Agregar leyenda y mostrar el gráfico
limiteC = np.ones(len(distancia)) * 100
plt.plot(distancia, limiteC, color='red', label=f"Límite CIGRE 100 nA/$m^{2}$")
plt.legend(loc="upper left", fontsize=10)
df2 = pd.read_excel(archivo, sheet_name='Conductores')
columnas_a_extraer_2 = ['Posición en x (m)', 'Posición en y (m)', 'Posición en x 2(m)', 'Posición en y 2(m)']
df2 = df2[columnas_a_extraer_2]
df3 = pd.read_excel(archivo, sheet_name='Dim y discr')
columnas_a_extraer_3 = ['Altura (m)']
ymax_ax = float(df3[columnas_a_extraer_3].iloc[0])
print(f"ymax es {ymax_ax}")
if any(x < 0 for x in campo_electrico):
    escala = [-ymax_ax, ymax_ax]
else:
    escala = [0, ymax_ax]
x1 = df2['Posición en x (m)'].iloc[0]
y1 = df2['Posición en y (m)'].iloc[0]
x2 = df2['Posición en x 2(m)'].iloc[0]
y2 = df2['Posición en y 2(m)'].iloc[0]
plt.legend(loc="lower left")
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.scatter(x1, y1, color='purple', marker='o', label=f"Polos altura {pos_y} m")  # Agregar el punto
ax2.scatter(x2, y2, color='purple', marker='o')  # Agregar el punto
ax2.set_ylabel('Altura (m)', fontsize=11, color='purple')
ax2.tick_params(axis='y', labelcolor='purple')  # Ajustar el color de las etiquetas del eje secundario
ax2.set_ylim(escala[0], escala[1])  # Definir el rango de 0 a 20
ax2.legend(loc="upper right")
plt.tight_layout()
plt.show()