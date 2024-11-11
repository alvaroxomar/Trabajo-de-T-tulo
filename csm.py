import numpy as np
import matplotlib.pyplot as plt
import math as ma
from matplotlib.colors import LogNorm 

## Parámetros genericos
R = 6.35 * 10**(-2) # (m) radio equivalente del conductor
epsilon0 = (1/(36*np.pi)) * 10**(-9) # (F/m) permitividad del vacío
Vol = 500000 # (V) voltaje de la línea
K = 1/(2*np.pi*epsilon0) # factor de multiplicación
mov = 1.5*10**(-1) # (m/Vs^2)
m = 1 # (AD) factor de rugosidad
P0 =101.3 # (kPa) Presión del aire a nivel de mar
T0 = 303 # (Kelvin) Temperatura de 25°C
Pr =  90 # Presión del aire
Tr= 290 # (Kelvin) Temperatura del sitio
delta = Pr*T0/(P0*Tr) # () densidad del aire
Sx = 10 # (m) media longitud del plano de tierra 
Sy = 10 # (m) altura del área de estudio respecto de tierra
l = 0.3 # (m) distancia desde el suelo a la altura de interés de influencia de campo
## Definición coordenadas conductores caso bipolar
#coordenada = [(6,7), (-6,7)] # (m) coordenadas cargas en posición de los conductores
coordenada = [(0,7)]
coordenada_im =  [(x, -y) for x, y in coordenada] # (m) coordenadas de las cargas imágenes
h = np.array([y for x,y in coordenada]) # (m) alturas de los conductores
w = np.array([x for x,y in coordenada]) # (m) anchos de los conductores
#Hor = np.linspace(-10,10,100) # (m) Ancho total del área de estudio
#Ver = 1 # (m) Altura donde interesa medir el campo eléctrico
def mod(z1,z2):
    return np.sqrt((z1[0]-z2[0])**2 + (z1[1]-z2[1])**2)

largo = len(coordenada)
D= np.zeros((largo,largo))
D_im= np.zeros((largo,largo))
for i in range(len(coordenada)):
    for j in range(len(coordenada)):
        if i != j:
            D[i][j]=mod(coordenada[i], coordenada[j])
            D_im[i][j]=mod(coordenada[i], coordenada_im[j])

## Coeficientes de potencial
P = np.zeros((largo,largo)) # matriz coeficientes
for i in range(largo):
    for j in range(largo):
        if j==i:
            P[i][j] = K*np.log(2*h[i]/R)
        else:
            if largo !=1 :
                P[i][j] = K*np.log(D_im[i][j]/D[i][j])
#V = [Vol,Vol] # voltajes sobre los conductores reales
V = [Vol]
Q = np.dot(np.linalg.inv(P),V) # Se determinan las cargas de los conductores

### Caso unipolar

## Obtención campo eléctrico
E = np.zeros((100,100))
Exx = np.zeros((100,100))
Eyy = np.zeros((100,100))
#Eu = np.zeros((200,200))
#Exxu = np.zeros((200,200))
#Eyyu = np.zeros((200,200))
x = np.linspace(-Sx,Sx,len(E[0,:]))
y = np.linspace(Sy, 0, len(E[:,0]))
dx = np.abs(x[0]-x[1])
dy = np.abs(y[0]-y[1])
for i in range(len(x)):
    for j in range(len(y)):
        Ex=0
        Ey=0
        Exu=0
        Eyu=0
        for z in range(largo):
            N1x= x[i]-w[z]
            N1y= y[j]-h[z]
            N2y = y[j]+h[z]
            Ex += Q[z]*K*(N1x/(N1x**2+N1y**2)-N1x/(N1x**2+N2y**2))
            Ey += Q[z]*K*(N1y/(N1x**2+N1y**2) - N2y/(N1x**2+N2y**2))
            #Exu += -Q[z]*(K/2)*(N1x/(N1x**2+N2y**2)**(3/2) - N1x/(N1x**2+N1y**2)**(3/2))
            #Eyu += -Q[z]*(K/2)*(N2y/(N1x**2+N2y**2)**(3/2) - N1y/(N1x**2+N1y**2)**(3/2))
        Exx[j][i] = Ex
        Eyy[j][i] = Ey
        E[j][i] = np.sqrt(Ex**2 + Ey**2)
        #Exxu[j][i] = Exu
        #Eyyu[j][i] = Eyu
        #Eu[j][i] = np.sqrt(Exu**2 + Eyu**2)

X,Y =np.meshgrid(x,y)
U = Exx
Ww = Eyy
#Uu = Exxu
#Wu = Eyyu
'''
#######
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Primer gráfico de campo vectorial (con barra de color)
quiver1 = ax1.quiver(X, Y, U, Ww, E, cmap='inferno')
ax1.set_title('Campo Vectorial 1 (Circular)')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
colorbar1 = fig.colorbar(quiver1, ax=ax1)
colorbar1.set_label('Magnitud')

# Segundo gráfico de campo vectorial (con barra de color)
quiver2 = ax2.quiver(X, Y, Uu, Wu, Eu, cmap='plasma')
ax2.set_title('Campo Vectorial 2 (Ondulatorio)')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
colorbar2 = fig.colorbar(quiver2, ax=ax2)
colorbar2.set_label('Magnitud')

# Ajustar el espacio entre los gráficos
plt.tight_layout()
#######
'''
plt.figure(1)
plt.quiver(X, Y, U, Ww, E, cmap='plasma')
#plt.imshow(E, cmap='viridis', interpolation='none')

# Agregar la barra de colores
cbar = plt.colorbar()
cbar.set_label(r'Magnitud campo eléctrico $kV/m$', fontsize=11)

# Obtén los ticks actuales
ticks = cbar.get_ticks()

# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
# Mostrar el gráfico
plt.title('Campo electrostático libre de iones', fontsize=15)
#plt.show()



##### Resolución ecuación de continuidad
Jp = 2.4*10**(-3) # (A/m^2) Densidad de corriente iónica promedio sobre el plano de tierra (Se debe adivinar este valor)
# Condiciones de borde
rho_i = (Sx*Jp)*10**(-3)/(np.pi*mov*(100*R)*m*(30*delta + 9*np.sqrt(delta/(100*R)))) # radio está en cm
rho = np.zeros((100,100))*10**(-8)
posx_conductor = int(Sx/dx)
posy_conductor = int((Sy - coordenada[0][1])/dy)
#rho[posy_conductor][posx_conductor] = rho_i
rho[-1][:] = np.abs(Jp/(mov*E[-1][:]))
## Resolución para densidad de carga en el espacio
con = 1
c = 0
rho_n = rho.copy() # se copia el primer valor
#while np.all(con > 10**(-7)):
#    print(rho_n)
while True:
    if c<200:
        for i in range(1, len(rho[:,0])-1):
            for j in range(1, len(rho[0,:])-1):
                if i == posy_conductor and j == posx_conductor:
                    rho_n[i][j] = rho_i
                else:
                    #alpha = 0.5*(np.abs(Exx[i][j])/dx + np.abs(Eyy[i][j])/dy)*epsilon0
                    #beta = epsilon0*(np.abs(Eyy[i][j])*rho_n[i-1][j]/dy + np.abs(Exx[i][j])*rho_n[i][j-1]/dx)
                    #alpha = 0.5*(Exx[i][j]/dx + Eyy[i][j]/dy)*epsilon0
                    #beta = epsilon0*(Eyy[i][j]*rho_n[i-1][j]/dy + Exx[i][j]*rho_n[i][j-1]/dx)
                    #rho_n[i][j] = -alpha + np.sqrt(np.abs(alpha**2 + beta)) # Se toma solución positiva dado que se asume que todas los iones son de carga positiva
                    a = rho_n[i][j+1]- rho_n[i][j-1]
                    b = rho_n[i+1][j]- rho_n[i-1][j]
                    rho_n[i][j] = -np.sqrt(0.5*epsilon0*(np.abs(Exx[i][j]*(-a))/dx + np.abs(Eyy[i][j]*(-b))/dy))
                    #if alpha**2 + beta < 0:
                    #    print(alpha**2 + beta)
                    #    break
        #print(r'rho_n: '+str(rho_n))
        #print(r'rho: '+str(rho))
        #print(r'iteracion: '+str(c))
        con = np.abs(rho_n - rho)
        #print(r'conv: '+str(con))
        c +=1 
        rho = rho_n.copy()
    else:
        break
# Mostrar la matriz en una gráfica con un mapa de colores
plt.figure(2)
#plt.imshow(rho,extent=[x[0], x[-1], y[-1], y[0]], cmap='viridis', interpolation='none',norm=LogNorm())
plt.pcolormesh(X, Y, rho, cmap='viridis', shading='auto',norm=LogNorm())
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
plt.title('Densidad de carga',fontsize=15)
# Añadir una barra de colores para mostrar la escala
cbar = plt.colorbar()
cbar.set_label(r'Densidad de carga $C/m^3$',fontsize=11)
#ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
#cbar.set_ticks(ticks)
#cbar.set_ticklabels([f'{tick/1000000:.1f}' for tick in ticks]) 
#plt.show()


#### Calculo del potencial con FMG

Vm = np.zeros((100,100))
## Condiciones de borde
Vm[posy_conductor][posx_conductor] = Vol
Vm[-1][:] = 0
Vm[0][:] = Q*(0.5*K/np.sqrt((x-w)**2+(y[0]-h)**2))
Vm[0:-1,0] = Q*(0.5*K/np.sqrt((x[0]-w)**2+(y[0:-1]-h)**2))
Vm[0:-1,-1] = Q*(0.5*K/np.sqrt((x[-1]-w)**2+(y[0:-1]-h)**2))
print(Vm)
Vm_n = Vm.copy()
cd=0
while True:
    if cd < 200:
        for i in range(1, len(Vm[0][:])-1):
            for j in range(1, len(Vm[:][0])-1):
                if i==posy_conductor and j == posx_conductor:
                    Vm_n[i][j] = Vol
                    #Vm_n[i][j] = ((dx**2)*(Vm_n[i+1][j]+Vm_n[i-1][j])+(dy**2)*(Vm_n[i][j+1]+Vm_n[i][j-1])+((dx*dy)**2)*rho[i][j]/epsilon0)/(2*(dx**2+dy**2))
                else:
                    Vm_n[i][j] = ((dx**2)*(Vm_n[i+1][j]+Vm_n[i-1][j])+(dy**2)*(Vm_n[i][j+1]+Vm_n[i][j-1])+((dx*dy)**2)*rho[i][j]/epsilon0)/(2*(dx**2+dy**2))
                    #Vm_n[i][j] = ((dx**2)*(Vm_n[i+1][j]+Vm_n[i-1][j])+(dy**2)*(Vm_n[i][j+1]+Vm_n[i][j-1]))/(2*(dx**2+dy**2))
        Vm = Vm_n.copy()
        cd+=1
    else:
        break
print(r'nuevo V: '+str(Vm))
plt.figure(3)
#plt.imshow(Vm,extent=[x[0], x[-1], y[-1], y[0]], cmap='plasma', interpolation='none',norm=LogNorm())
plt.pcolormesh(X, Y, Vm, cmap='plasma', shading='auto')
plt.title('Potencial',fontsize=15)
plt.xlabel('Distancia horizontal (m)',fontsize=11)
plt.ylabel('Distancia vertical (m)',fontsize=11)
# Añadir una barra de colores para mostrar la escala
cbar = plt.colorbar()
cbar.set_label(r'Potencial $kV$')
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
#plt.show()


##### Cálculo de campo eléctrico definitivo
##### Cálculo de densidad de corriente iónica
## Ocupa el negativo del gradiente del potencial calculado

Edef = np.zeros((100,100))
Edefx = np.zeros((100,100))
Edefy = np.zeros((100,100))
J = np.zeros((100,100))
Ei = np.zeros(100)
Ji = np.zeros(100)
Jpp = []
for i in range(1, len(Edef[0][:])-1):
    for j in range(1,len(Edef[:][0])-1):
        xx = -(Vm[i][j+1]-Vm[i][j-1])/(2*dx)
        yy = (Vm[i+1][j]-Vm[i-1][j])/(2*dy)
        Edefx[i][j] = xx
        Edefy[i][j] = yy
        Edef[i][j] = np.sqrt(xx**2 + yy**2)
        J[i][j] = rho[i][j]*mov*Edef[i][j]
        if i == int((Sy - l)/dy):
            Ei[j] = Edef[i][j]
            Ji[j] = J[i][j]
        elif i == int(len(Edef[0][:])-2):
            Jpp.append(J[i][j])
Jave  = np.mean(Jpp)
print(r'Jp promedio calculado: '+str(Jave))
Xe,Ye =np.meshgrid(x[1:-1],y[1:-1])
plt.figure(4)
plt.quiver(Xe, Ye, Edefx[1:-1, 1:-1], Edefy[1:-1, 1:-1], Edef[1:-1, 1:-1], cmap='plasma')
# Agregar la barra de colores
cbar = plt.colorbar()
cbar.set_label(r'Magnitud campo eléctrico $kV/m$', fontsize=11)
# Obtén los ticks actuales
ticks = cbar.get_ticks()
# Cambia las etiquetas a una magnitud diferente (por ejemplo, dividiendo por 1000)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{tick/1000:.1f}' for tick in ticks]) 
# Mostrar el gráfico
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Distancia vertical (m)',fontsize=11)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo (V/m)', fontsize=13)
plt.tight_layout()


plt.figure(5)
plt.plot(x[1:-1], Ji[1:-1])
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Densidad de corriente iónica ($A/m^2$)',fontsize=11)
plt.title(r'Magnitud de corriente iónica a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.grid(True)

plt.figure(6)
plt.plot(x[1:-1], Ei[1:-1])
plt.xlabel(r'Distancia horizontal (m)',fontsize=11)
plt.ylabel(r'Campo eléctrico (V/m)',fontsize=11)
plt.title(r'Magnitud de campo eléctrico a nivel de suelo, $l=$'+str(l)+' m', fontsize=13)
plt.tight_layout()
plt.grid(True)

plt.show()


