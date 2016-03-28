import numpy as np
import scipy.optimize as scipy
import matplotlib.pyplot as plt

#variables
gamma = 1.4

#condiciones iniciales Left and Right - con dimensiones
crhol = 1.0
cul = 0.0
cpl = 1.0
cel = cpl / (gamma - 1.0) + crhol * cul**2 / 2.0
cal = np.sqrt( gamma * cpl / crhol )
crhor = 0.125
cur = 0.0
cpr = 0.1
cer = cpr / (gamma - 1.0) + crhor * cur**2 / 2.0
car = np.sqrt( gamma * cpr / crhor)
 
#Condiciones iniciales sin dimensiones
rhol = crhol/crhor
ul = 0.0
pl = cpl / (gamma * cpr)
el = pl / (gamma - 1.0) + rhol * ul**2 / 2.0
al = cal/car 
rhor = 1.0
ur = 0.0
pr = 1.0 / gamma 
er = pr / (gamma - 1.0) + rhor * ur**2 / 2.0
ar = 1.0 

#variables para la ecuacion de ms:
v1 = (gamma + 1) / (gamma - 1)
v2 = 2 * gamma / (gamma + 1)
v4 = (gamma - 1) / (2 * gamma)

#funcion para hallar las raices (ms):
def funcion(x):
    return  1/x - x + al * v1 * ( 1.0 - ( pr/pl * ( v2 * x**2 - 1/v1 ) )** v4  )

#encontrar raices de la ecuacion
solucion = scipy.fsolve(funcion, 5)
ms = solucion[0]

#El tubo se separa en 5 intervalos: L, E, 2, 1, R
#En los intervalos L y R se mantienen las condiciones iniciales

#Utilizar el valor ms para hallar variables de los intervalos 1 y 2
u1 = ( ms - 1/ms ) * 2 / (gamma + 1)
u2 = u1
a2 = al - v4 * gamma * u2
p1 = pl * ( a2 / al )**( 1/v4 )
p2 = p1
rho2 = rhol * ( p2 / pl )**( 1.0 / gamma )
rho1 = rhor / ( v2 / (gamma * ms**2) + 1/v1 ) 
e1 = p1 / (gamma - 1) + rho1 * u1**2 /2
e2 = p2 / (gamma - 1) + rho2 * u2**2 /2

#El intervalo espacial se separa en (0, x1, x2, x3, x4)
# 0:x1 left,  x1:x2 expansion fan,  x2:x3 region 2,  x3:x4 region 1,  x4:1  right
x0 = 0.5
# t=1
t = 0.2
x1 = x0 - al*t
x2 = x0 + (u2 - a2)*t
x3 = x0 + u2*t
x4 = x0 + ms*t
print x1, x2, x3, x4

#Arreglos con informacion x,rho, u, p, e
N = 5000
x = np.linspace(0,1,5000)
posicion = np.linspace(1,4,5000)
rho = np.zeros(N)
u = np.zeros(N)
p = np.zeros(N)
e = np.zeros(N)

#escribir solucion
for i in range(N):
    #Left region
    if(x[i]<x1):
        rho[i] = crhol
        e[i] = cel
        u[i] = cul
        p[i] = cpl
    #Expansion fan entre x1 y x2
    if(x1 <= x[i] and x[i] < x2):
        #solucion dentro de la expansion
        ue = v2 / gamma * (al + x[i] -x0)
        ae = al - (gamma - 1.0) * ue /2.0
        pe = pl * (ae / al)**(1/v4)
        rhoe = gamma * pe / ae**2
        ee= pe / (gamma - 1.0) + rhoe * ue**2 /2.0
        rho[i] = rhoe * crhor
        e[i] = ee * crhor * ar**2
        u[i] = ue * car
        p[i] = pe * gamma * cpr
    #Intervalo 2
    if(x2 <= x[i] and x[i] < x3):
        rho[i] = rho2 * crhor
        e[i] = e2 * crhor * ar**2
        u[i] = u2 * car
        p[i] = p2 * gamma * cpr
    #Intervalo 1
    if(x3 <= x[i] and x[i] < x4):
        rho[i] = rho1 * crhor
        e[i] = e1 * crhor * ar**2
        u[i] = u1 * car
        p[i] = p1 * gamma * cpr
    # right region
    if(x4 <= x[i]):
        rho[i] = crhor
        e[i] = cer
        u[i] = cur
        p[i] = cpr
        


plt.figure(figsize=(10,10))
data=np.genfromtxt("UpwindGodunov_step_4.dat")
for j in range(1,5):
    plt.subplot(2,2,j)
    if(j==1):
        plt.plot(posicion,data[:,j],label="Godunov")
        plt.plot(posicion,rho,label=(u'Analitica'))
        plt.title("Densidad Rho")
        plt.legend()
    if(j==2):
        plt.plot(posicion,data[:,j],label="Godunov")
        plt.plot(posicion,u,label="Analitica")
        plt.title("Velocidad u")
        plt.legend()
    if(j==3):
        plt.plot(posicion,data[:,j],label="Godunov")
        plt.plot(posicion,e,label="Analitica")
        plt.title("Energia e")
        plt.legend()
    if(j==4):
        plt.plot(posicion,data[:,j],label="Godunov")
        plt.plot(posicion,p,label="Analitica")
        plt.title("Presion p")
        plt.legend()
plt.savefig("Comparacion_graficas")
plt.show()

#Graficas de la evolucion temporal de la solucion usando Godunov
for i in (range(5)):
    plt.figure(figsize=(10,10))
    data=np.genfromtxt("UpwindGodunov_step_" + str(i) + ".dat")
    for j in range(1,5):
        plt.subplot(2,2,j)
        plt.plot(data[:,0],data[:,j])
        if(j==1):
            plt.title("God rho")
        if(j==2):
            plt.title("God u")
        if(j==3):
            plt.title("God e")
        if(j==4):
            plt.title("God p")
    plt.savefig("Tiempo_" + str(i) + "_Godunov")
