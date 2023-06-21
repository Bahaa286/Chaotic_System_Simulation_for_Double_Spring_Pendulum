import matplotlib.pyplot as plt
import math

#constants ----------------------------------------------------
M1:float = 1
M2:float = 1
L1:float = 50.0
L2:float = 50.0
K1:float = 1000
K2:float = 1000
g:float = 9.81
h:float = 1e-2
N:int = 2e3

#initial values ------------------------------------------------
#particle 1:
x1_0:float = 30
y1_0:float = 0
z1_0:float = -40
vx1_0:float = 0
vy1_0:float = 0
vz1_0:float = 0
#particle 2:
x2_0:float = 60
y2_0:float = 30
z2_0:float = -90
vx2_0:float = 0
vy2_0:float = 0
vz2_0:float = 0

# ------------------------------------------------ functions ----------------------------------------------------
# functions for the equation vi' = fi(x1,y1,z1,x2,y2,z2)
def vx1(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    #print(x1, y1, z1, x2, y2, z2)
    return (K1*x1/M1)*(L1/(math.sqrt((x1)**2+(y1)**2+(z1)**2))-1)+(K2*(x1-x2)/M1)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1)

def vy1(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    return (K1*y1/M1)*(L1/(math.sqrt((x1)**2+(y1)**2+(z1)**2))-1)+(K2*(y1-y2)/M1)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1)

def vz1(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    return (K1*z1/M1)*(L1/(math.sqrt((x1)**2+(y1)**2+(z1)**2))-1)+(K2*(z1-z2)/M1)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1) - g

def vx2(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    return (K2*(x2-x1)/M2)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1)

def vy2(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    return (K2*(y2-y1)/M2)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1)

def vz2(x1:float, y1:float, z1:float, x2:float, y2:float, z2:float) ->float:
    return (K2*(z2-z1)/M2)*(L2/(math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))-1) - g

#functions for the equation qi' = vi
def fx1(vx1:float) ->float:
    return vx1

def fy1(vy1:float) ->float:
    return vy1

def fz1(vz1:float) ->float:
    return vz1

def fx2(vx2:float) ->float:
    return vx2

def fy2(vy2:float) ->float:
    return vy2

def fz2(vz2:float) ->float:
    return vz2

#--------------------------------- numerical solution -------------------------------------
# function for Runge Kutta method 
# r_list[x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2]
# k_list[kx1, ky1, kz1, kx2, ky2, kz2, kvx1, kvy1, kvz1, kvx2, kvy2, kvz2]
def k1(r_list:list[float]):
    # for vi
    kvx1:float = vx1(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])
    kvy1:float = vy1(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])
    kvz1:float = vz1(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])

    kvx2:float = vx2(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])
    kvy2:float = vy2(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])
    kvz2:float = vz2(r_list[0], r_list[1], r_list[2], r_list[3], r_list[4], r_list[5])

    # for xi
    kx1:float = fx1(r_list[6])
    ky1:float = fy1(r_list[7])
    kz1:float = fz1(r_list[8])

    kx2:float = fx2(r_list[9])
    ky2:float = fy2(r_list[10])
    kz2:float = fz2(r_list[11])

    return kx1, ky1, kz1, kx2, ky2, kz2, kvx1, kvy1, kvz1, kvx2, kvy2, kvz2

def k2(r_list:list[float],k_list:list[float]):

    # for vi
    kvx1:float = vx1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvy1:float = vy1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvz1:float = vz1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)

    kvx2:float = vx2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvy2:float = vy2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvz2:float = vz2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)

    # for xi 
    kx1:float = fx1(r_list[6]+h*k_list[6]/2)
    ky1:float = fy1(r_list[7]+h*k_list[7]/2)
    kz1:float = fz1(r_list[8]+h*k_list[8]/2)
    

    kx2:float = fx2(r_list[9]+h*k_list[9]/2)
    ky2:float = fy2(r_list[10]+h*k_list[10]/2)
    kz2:float = fz2(r_list[11]+h*k_list[11]/2)

    return kx1, ky1, kz1, kx2, ky2, kz2, kvx1, kvy1, kvz1, kvx2, kvy2, kvz2

def k3(r_list:list[float],k_list:list[float]):

    # for vi
    kvx1:float = vx1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvy1:float = vy1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvz1:float = vz1(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)

    kvx2:float = vx2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvy2:float = vy2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)
    kvz2:float = vz2(r_list[0]+h*k_list[0]/2, r_list[1]+h*k_list[1]/2, r_list[2]+h*k_list[2]/2, r_list[3]+h*k_list[3]/2, r_list[4]+h*k_list[4]/2, r_list[5]+h*k_list[5]/2)

    # for xi
    kx1:float = fx1(r_list[6]+h*k_list[6]/2)
    ky1:float = fy1(r_list[7]+h*k_list[7]/2)
    kz1:float = fz1(r_list[8]+h*k_list[8]/2)

    kx2:float = fx2(r_list[9]+h*k_list[9]/2)
    ky2:float = fy2(r_list[10]+h*k_list[10]/2)
    kz2:float = fz2(r_list[11]+h*k_list[11]/2)

    return kx1, ky1, kz1, kx2, ky2, kz2, kvx1, kvy1, kvz1, kvx2, kvy2, kvz2

def k4(r_list:list[float],k_list:list[float]):

    # for vi
    kvx1:float = vx1(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])
    kvy1:float = vy1(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])
    kvz1:float = vz1(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])

    kvx2:float = vx2(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])
    kvy2:float = vy2(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])
    kvz2:float = vz2(r_list[0]+h*k_list[0], r_list[1]+h*k_list[1], r_list[2]+h*k_list[2], r_list[3]+h*k_list[3], r_list[4]+h*k_list[4], r_list[5]+h*k_list[5])

    # for xi
    kx1:float = fx1(r_list[6]+h*k_list[6])
    ky1:float = fy1(r_list[7]+h*k_list[7])
    kz1:float = fz1(r_list[8]+h*k_list[8])

    kx2:float = fx2(r_list[9]+h*k_list[9])
    ky2:float = fy2(r_list[10]+h*k_list[10])
    kz2:float = fz2(r_list[11]+h*k_list[11])

    return kx1, ky1, kz1, kx2, ky2, kz2, kvx1, kvy1, kvz1, kvx2, kvy2, kvz2

# function to find next point
def next_position(r_list:list[float]):
    
    k1_list:list[float] = k1(r_list)
    k2_list:list[float] = k2(r_list, k1_list)
    k3_list:list[float] = k3(r_list, k2_list)
    k4_list:list[float] = k4(r_list, k3_list)

    new_r_list = []
    for i in range(12):
        new_r_list.append(r_list[i] + (h/6)*( k1_list[i] + 2*k2_list[i] + 2*k3_list[i] + k4_list[i]))

    return new_r_list

#--------------------------------- results and solution -------------------------------------
xl1:list[float] = []
yl1:list[float] = []
zl1:list[float] = []
xl2:list[float] = []
yl2:list[float] = []
zl2:list[float] = []

vxl1:list[float] = []
vyl1:list[float] = []
vzl1:list[float] = []
vxl2:list[float] = []
vyl2:list[float] = []
vzl2:list[float] = []

def clear_lists():
    xl1.clear()
    yl1.clear()
    zl1.clear()
    xl2.clear()
    yl2.clear()
    zl2.clear()

    vxl1.clear()
    vyl1.clear()
    vzl1.clear()
    vxl2.clear()
    vyl2.clear()
    vzl2.clear()

def result():
    # r_list[x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2]
    r_list:list[list[float]] = []
    r_list = [x1_0, y1_0, z1_0, x2_0, y2_0, z2_0, vx1_0, vy1_0, vz1_0, vx2_0, vy2_0, vz2_0]

    for i in range(int(N)):
        # position lists
        xl1.append(r_list[0])
        yl1.append(r_list[1])
        zl1.append(r_list[2])
        xl2.append(r_list[3])
        yl2.append(r_list[4])
        zl2.append(r_list[5])

        # velocities list
        vxl1.append(r_list[6])
        vyl1.append(r_list[7])
        vzl1.append(r_list[8])
        vxl2.append(r_list[9])
        vyl2.append(r_list[10])
        vzl2.append(r_list[11])

        r_list = next_position(r_list=r_list)

# 3d plot-----------------------------------------------------------
def plot3d(particle1=False, particle2=False):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    if particle1:
        ax.plot3D(xl1, yl1, zl1, label='particle 1')
    if particle2:
        ax.plot3D(xl2, yl2, zl2, label='particle 2')
    ax.set_xlabel('x', size=18)
    ax.set_ylabel('y', size=18)
    ax.set_zlabel('z', size=18)
    ax.set_title('two particle motion', size=20)
    plt.legend()
    plt.show()

# 2d plot-----------------------------------------------------------
def plot2d(particle1=False, particle2=False):
    if particle1:
        plt.plot(xl1,zl1, label='particle 1')
    if particle2:
        plt.plot(xl2,zl2, label='particle 2')
    plt.title('z vs x', size=20)
    plt.xlabel('x', size=18)
    plt.ylabel('z', size=18)
    plt.legend()
    plt.show()


# find the result and fill the lists
result()
# plot the result
plot3d(particle1=False, particle2=True)
clear_lists()

def change_in_particle1():
    global x1_0
    x1_0 = 30.01
    global y1_0
    y1_0 = 0
    global z1_0
    z1_0 = -40

def change_in_particle2():
    global x2_0
    x2_0 = 70.01
    global y2_0
    y2_0 = 30
    global z2_0
    z2_0 = -90


# plot for twe results
def plot_twe_results():
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    #plot 
    result()
    ax.plot3D(xl2, yl2, zl2, label='result 1')
    clear_lists()
    change_in_particle1()
    result()
    ax.plot3D(xl2, yl2, zl2, label='result 2')

    ax.set_xlabel('x', size=18)
    ax.set_ylabel('y', size=18)
    ax.set_zlabel('z', size=18)
    ax.set_title('one particle two different paths', size=20)
    plt.legend()
    plt.show()

plot_twe_results()
