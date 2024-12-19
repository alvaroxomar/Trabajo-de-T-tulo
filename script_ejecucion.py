if __name__ == "__main__":
        ##############################
    ###### LIBRERIAS IMPORTADAS
    ##############################
    import numpy as np
    import matplotlib.pyplot as plt
    import math as ma
    from matplotlib.colors import LogNorm
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.interpolate import griddata
    from scipy.ndimage import laplace
    import json
    import sys
    ########## FUNCIONES ###############
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
        r_eq = distancia * (radio / distancia)**(1 / numero)
        return 1.09 * r_eq if numero == 4 else r_eq
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
        lados = calcula_lados(numero, sepa)
        distancia = distancia_equivalente(lados)
        radio = convierte_a_radio(area_sub, conversion, es_mcm)
        if es_cm:
            radio /= 10  # Convertir de mm a cm si es necesario
        radio_equi = radio_eq(radio, numero, distancia)
        return radio_equi, radio
    
        # Encontrar los índices más cercanos
    def encuentra_nodos(x, y, x0,y0):
        indice_x0 = (np.abs(x - x0)).argmin()
        indice_y0 = (np.abs(y - y0)).argmin()
        return indice_x0, indice_y0

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

    def discretiza(min, nodox=None, nodoy=None, sx=None, sy=None):
        if (sx is None and sy is None) and (nodox is not None and nodoy is not None):
            Sx, Sy = espacio(min, min, nodox, nodoy)
            x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=nodox, nodoy=nodoy)
            distx=np.abs(x_coor)-Sx
            nodos_sx = int(distx/min) + 1
            disty=np.abs(y_coor)-Sy
            nodos_sy = int(disty/min) + 1
            assert np.abs(x_coor) < Sx, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sx+nodosx}, o bien ubica a menor distancia x_coor, la dif x-Sx es {distx}'
            assert np.abs(y_coor) < Sy, f'bajo dimensionamiento, utiliza un radio mayor, o selecciona más nodos: {nodos_sy+nodosy}, o bien ubica a menor distancia y_coor, la dif x-Sx es {disty}'
        elif (sx is not None and sy is not None) and (nodox is None and nodoy is None):
            Sx, Sy = sx, sy
            x, y, nodosx, nodosy = malla(min, Sx, Sy, nodox=None, nodoy=None)
        elif (sx is not None and sy is not None) and (nodox is not None and nodoy is not None):
            Sx, Sy = sx, sy
            x, y, nodosx, nodosy = malla(min, sx, sy, nodox=nodox,nodoy=nodoy)
        else:
            print('Los parámetros no están bien ajustados en la función')
        return x, y, nodosx, nodosy, Sx, Sy

    ############### PARÁMETROS ###################
    # Capturar el diccionario JSON desde los argumentos
    if len(sys.argv) < 2:
        print("Se requiere un diccionario JSON como argumento.")
        sys.exit(1)

    # Cargar el diccionario desde el JSON
    try:
        params = json.loads(sys.argv[1])
    except (IndexError, json.JSONDecodeError) as e:
        print(f"Error al leer los parámetros: {e}")
    # Desestructurar el diccionario en variables específicas
    
    # Extraer parámetros necesarios
    area = float(params["Área (mm²)"])
    cantidad = int(params["Cantidad"])
    sep = float(params["Separación (cm)"])
    es_cm = True # si separacion está en cm entonces es True
    es_mcm = False # si  aárea está en mm2 entonces es False
    conversion = 0.5067  # Conversión de MCM a mm²
    Vol = float(params["Voltaje (kV)"])
    x_coor = float(params["Posición en x (m)"])
    y_coor = float(params["Posición en y (m)"])
    # Obtener los parámetros
    Sx = obtener_parametro(params["Ancho (m)"], float)  # (m) Media longitud del plano de tierra
    Sy = obtener_parametro(params["Altura (m)"], float)  # (m) Altura del área de estudio respecto de tierra
    nodosx = obtener_parametro(params["nodos x"], int)  # Número de nodos en x
    nodosy = obtener_parametro(params["nodos y"], int)  # Número de nodos en y
    mov = float(params["Movilidad iónica (m2/kVs)"])
    Tem = float(params["Temperatura (°C)"])
    Pres = float(params["Presión (Pa)"])
    hum =  float(params["Humedad (%)"])
    viento = float(params["Viento (m/s)"])
    modo = str(params["Modo (str)"])
    rug =  float(params["Rugosidad terreno"])
    m = float(params["factor conductor"])
    l = float(params["Medición (m)"])

    print(str(area))
    print(str(sep))
    Req, R = calculo_radio_eq(cantidad, area, sep, conversion=conversion, es_mcm=es_mcm, es_cm=es_cm) # están en cm
    print(f"radio eq {Req} y radio subconductor {R}")
    Req /= 100 # en m
    R /= 100
    print(f'{Req},{R} en m')
    '''
    Jp_inicial =params.get("Jp_inicial", 1988e-8)
    wndx = params.get("wndx", 0)
    wndy = params.get("wndy", 0)
    modo = params.get("modo", "uniforme")
    coef_drag = params.get("coef_drag", 0.3)
    TolDev = params.get("TolDev", 1e-2)
    in_condct = params.get("in_condct", "si")
    copiado = params.get("copiado", "no")
    mostrar = params.get("mostrar", False)
    visualizacion = params.get("visualizacion", 15)
    l = params.get("l", 1)
    m = params.get("m", 1)
    delta = params.get("delta", 1)
    max_iter_rho = params.get("max_iter_rho", 400)
    max_iter = params.get("max_iter", 250)
    it_global = params.get("it_global", 10)
    '''
    # Parámetros constructivos
    fixed_value = Vol
    P0 =101.3 # (kPa) Presión del aire a nivel de mar
    T0= 298.15  # (Kelvin) Temperatura de 25°C  + 273.15
    delta = Pres*T0/(P0*Tem) # () densidad del aire
    epsilon0 = (1 / (36 * np.pi)) * 10**(-9)  # (F/m) permitividad del vacío
    K = 1 / (2 * np.pi * epsilon0)  # factor de multiplicación
    # Definir coordenadas de la malla en x y y
    x,y, nodosx, nodosy, Sx, Sy = discretiza(Req, nodox=nodosx, nodoy=nodosy, sx=Sx, sy=Sy)  
    X,Y =np.meshgrid(x,y)
    dx = np.abs(x[1] - x[0]) # (m) distancia física entre nodos en cada columna
    dy = np.abs(y[1] - y[0]) # (m) distancia física entre nodos en cada fila
    posx_conductor, posy_conductor = encuentra_nodos(x, y, x_coor, y_coor)
    fixed_point = (posy_conductor, posx_conductor)
    print(f"malla X {X}")
    print("________________")
    print(f"malla Y {Y}")
    print(f"nodos x: {nodosx}, nodos y: {nodosy}")
    print(f"Sx: {Sx}, Sy: {Sy}")
    print(f"dx: {dx}, dy: {dy}")
