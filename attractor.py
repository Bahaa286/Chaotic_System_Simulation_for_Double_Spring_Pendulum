from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import math

dt:float = 1e-3
steps:int=100000

#-------------------------------------------------------------
#Lorentze 
sigma:float=10
rho:float=28
beta:float=8/3

def xLorentz(x:float, y:float, z:float) -> float:
    return x + sigma*(y-x)*dt

def yLorentz(x:float, y:float, z:float) -> float:
    return y + (-x*z+rho*x-y)*dt

def zLorentz(x:float, y:float, z:float) -> float:
    return z + (x*y-beta*z)*dt
#-------------------------------------------------------------
#Rossler
a_Rossler=0.2 
b_Rossler=0.2  
c_Rossler=5.7

def xRossler(x:float, y:float, z:float) -> float:
    return x + (-y-z)*dt

def yRossler(x:float, y:float, z:float) -> float:
    return y + (x+a_Rossler*y)*dt

def zRossler(x:float, y:float, z:float) -> float:
    return z + (b_Rossler + z*(x-c_Rossler))*dt

#-------------------------------------------------------------

xl:list[float] = []
yl:list[float] = []
zl:list[float] = []

def findAttractor(x, y, z, x0:float, y0:float, z0:float):
    xn:float = x0
    yn:float = y0
    zn:float = z0

    tempx:float=x0
    tempy:float=y0
    tempz:float=z0

    xl.clear()
    yl.clear()
    zl.clear()

    for i in range(steps):
        tempx=xn
        tempy=yn
        tempz=zn

        xn = x(tempx, tempy, tempz)
        yn = y(tempx, tempy, tempz)
        zn = z(tempx, tempy, tempz)

        xl.append(xn)
        yl.append(yn)
        zl.append(zn)

findAttractor(xLorentz, yLorentz, zLorentz, 3, 1.4, 10.2)
findAttractor(xRossler, yRossler, zRossler, 10, 0, 10)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(xl, yl, zl)
plt.title('Rossler attractor', size=40)
plt.show()

    

