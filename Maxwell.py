import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

def minmod(a,b):
    global k
    if abs(a) < abs(b) and (a*b)>0:
        k = a
    elif abs(b)<abs(a) and (a*b)>0:
        k = b
    elif (a*b)<=0:
        k = 0
    return k

def gridfunction(nx):
    x = np.linspace(0, 10, nx)
    y = np.zeros((len(x)-1,1))
    for i in range(0, len(x)-1):
        y[i] = (x[i+1]+x[i])/2
    return x,y

def stepsizes(nx, nt, endtime):
    dx = 10/(nx-1)
    dt = endtime/nt
    return dx, dt  

def initialf(y):
    f = np.zeros((len(y),1))
    for i in range(0, len(y)):
        if y[i]>2 and y[i]<4:
            f[i] = 1
        if y[i]>6 and y[i]<8:
            f[i] = (1-np.cos(np.pi*y[i]))/2
    return f

def extendedinitial(y, time):
    #f1:f(x-1/2t)
    f1 = np.zeros((len(y),1))
    for i in range(0, len(y)):
        if y[i]>(2+int(time*0.5))%10 and y[i]<(4+int(time*0.5)-1e-15)%10:
            f1[i] = 1
        if y[i]>(6+int(time*0.5))%10 and y[i]<(8+int(time*0.5)-1e-15)%10:
            f1[i] = (1-np.cos(np.pi*y[i]))/2
    #f2:f(x+1/2t)
    f2 = np.zeros((len(y),1))
    for i in range(0, len(y)):
        if y[i]>(2-int(time*0.5))%10 and y[i]<(4-int(time*0.5)-1e-15)%10:
            f2[i] = 1
        if y[i]>(6-int(time*0.5))%10 and y[i]<(8-int(time*0.5)-1e-15)%10:
            f2[i] = (1-np.cos(np.pi*y[i]))/2
    return f1, f2   

#initial for upwind
def initialw(y):
    f = initialf(y)
    w10 = 0.5*np.ones((len(f),1))+f
    w20 = 0.5*np.ones((len(f),1))-f
    return w10, w20

#initial for fromm
def initialwfromm(y):
    f = initialf(y)
    w10 = 0.5*np.ones((len(f),1))+f
    w10 = np.vstack((w10, w10[1,:]))
    w10 = np.vstack((w10[-3,:], w10))
    w10 = np.vstack((w10[-4,:], w10))
    w20 = 0.5*np.ones((len(f),1))-f
    w20 = np.vstack((w20[-2,:], w20))
    w20 = np.vstack((w20, w20[2,:]))
    w20 = np.vstack((w20, w20[3,:]))
    return w10, w20    

#initial for maxwell
def initialwhighres(y):
    f = initialf(y)
    w10 = 0.5*np.ones((len(f),1))+f
    w10 = np.vstack((w10, w10[1,:]))
    w10 = np.vstack((w10[-3,:], w10))
    w10 = np.vstack((w10[-4,:], w10))
    w20 = 0.5*np.ones((len(f),1))-f
    w20 = np.vstack((w20[-2,:], w20))
    w20 = np.vstack((w20, w20[2,:]))
    w20 = np.vstack((w20, w20[3,:]))
    return w10, w20  




# upwind

def upwindmaxwell(nx, nt, endtime):
    dx, dt = stepsizes(nx, nt, endtime)
    x,y = gridfunction(nx)
    w1 = np.zeros((len(y),nt+1))
    w2 = np.zeros((len(y),nt+1))
    w10, w20 = initialw(y)
    cfl = 0.5*dt/dx
    for i in range(0,len(y)):
        w1[i,0] = w10[i]
        w2[i,0] = w20[i]
    for i in range(0, nt):
        for j in range(1, len(y)):
            w1[j,i+1] = w1[j,i] - 0.5*dt/dx*(w1[j,i]-w1[j-1,i])
            w1[0, i+1] = w1[len(y)-1,i+1]
        for j in range(0, len(y)-1):
            w2[j,i+1] = w2[j,i] + 0.5*dt/dx*(w2[j+1,i]-w2[j,i])
            w2[len(y)-1, i+1] = w2[0,i+1]
    # correct up untill this part
    E2 = 0.5*(w1-w2)
    B3 = w1+w2
    return y, cfl, E2, B3

def frommmaxwell(nx, nt, endtime):
    dx, dt = stepsizes(nx, nt, endtime)
    x,y = gridfunction(nx)
    w1 = np.zeros((len(y)+3,nt+1))
    w2 = np.zeros((len(y)+3,nt+1))
    w10, w20 = initialwfromm(y)
    cfl = 0.5*dt/dx
    for i in range(0,len(w10)):
        w1[i,0] = w10[i]
        w2[i,0] = w20[i]
    for i in range(0, nt):
        for j in range(2, len(w10)-1):
            w1[j, i+1] = w1[j, i] - 0.5*dt/(4*dx)*(w1[j+1, i] + 3*w1[j, i]-5*w1[j-1, i]+w1[j-2, i]) + (0.5*dt/dx)**2/4*(w1[j+1, i] - w1[j, i] - w1[j-1, i] + w1[j-2, i])
            w1[0, i+1] = w1[len(w1)-4,i+1]
            w1[1, i+1] = w1[len(w1)-3,i+1]
            w1[len(w1)-1,i+1] = w1[3,i+1]
        for j in range(1, len(w20)-2):
            w2[j, i+1] = w2[j, i] - 0.5*dt/(4*dx)*(w2[j-1, i] + 3*w2[j, i]-5*w2[j+1, i]+w2[j+2, i]) + (0.5*dt/dx)**2/4*(w2[j-1, i] - w2[j, i] - w2[j+1, i] + w2[j+2, i])
            w2[0, i+1] = w2[len(w2)-4, i+1]
            w2[len(w2)-2, i+1] = w2[2, i+1]
            w2[len(w2)-1, i+1] = w2[3, i+1]
    w1=np.delete(w1, [0,1, len(w1)-1], axis=0)
    w2=np.delete(w2, [0, len(w2)-1, len(w2)-2], axis=0)
    
    # correct up untill this part
    E2 = 0.5*(w1-w2)
    B3 = w1+w2
    return y, cfl, E2, B3

def highresmaxwell(nx, nt, endtime):
    dx, dt = stepsizes(nx, nt, endtime)
    x,y = gridfunction(nx)
    w1 = np.zeros((len(y)+3,nt+1))
    w2 = np.zeros((len(y)+3,nt+1))
    w10, w20 = initialwhighres(y)
    cfl = 0.5*dt/dx
    for i in range(0,len(w10)):
        w1[i,0] = w10[i]
        w2[i,0] = w20[i]
    for i in range(0, nt):
        for j in range(2, len(w10)-1):
            w1[j, i+1] = w1[j,i] - 0.5*dt/(dx)*(w1[j,i]-w1[j-1,i]) - 0.5*0.5*dt/(dx)*(dx-0.5*dt)*(minmod((w1[j, i]-w1[j-1, i])/dx,(w1[j+1, i]-w1[j, i])/dx) - minmod((w1[j-1, i]-w1[j-2, i])/dx,(w1[j, i]-w1[j-1, i])/dx))
            w1[0, i+1] = w1[len(w1)-4,i+1]
            w1[1, i+1] = w1[len(w1)-3,i+1]
            w1[len(w1)-1,i+1] = w1[3,i+1]
        for j in range(1, len(w20)-2):
            w2[j, i+1] = w2[j,i] - 0.5*dt/(dx)*(w2[j,i]-w2[j+1,i]) - 0.5*0.5*dt/(dx)*(dx-0.5*dt)*(minmod((w2[j, i]-w2[j-1, i])/dx,(w2[j+1, i]-w2[j, i])/dx) - minmod((w2[j+1, i]-w2[j, i])/dx,(w2[j+2, i]-w2[j+1, i])/dx))
            w2[0, i+1] = w2[len(w2)-4, i+1]
            w2[len(w2)-2, i+1] = w2[2, i+1]
            w2[len(w2)-1, i+1] = w2[3, i+1]
    w1=np.delete(w1, [0,1, len(w1)-1], axis=0)
    w2=np.delete(w2, [0, len(w2)-1, len(w2)-2], axis=0)
    
    # correct up untill this part
    E2 = 0.5*(w1-w2)
    B3 = w1+w2
    return y, cfl, E2, B3, w1, w2

def plot_maxwell(E2, B3,y,T, dt, Eexact, Bexact):
   plt.figure(1)
   moment = int(T/dt)
   plt.plot(y,E2[:,moment],'r')
   plt.xlabel('x')
   plt.ylabel('E2')
   plt.ylim([0,1.5])
   plt.title('Electric field over time')
   plt.hold(True)
   plt.plot(y,Eexact,'b')
   
   
   plt.figure(2)
   moment = int(T/dt)
   plt.plot(y,B3[:,moment],'r')
   plt.xlabel('x')
   plt.ylabel('B3')
   plt.ylim([-2,2.5])
   plt.title('Magnetic field over time')
   plt.hold(True)
   plt.plot(y,Bexact,'b')
   
def exactsolution(y, T):
   f1, f2 = extendedinitial(y, T)
   Eexact = 0.5*(f1+f2)
   Bexact = np.ones((len(f1),1)) + f1-f2
   return Eexact, Bexact

def errors(E2, B3, Eexact, Bexact, moment):
    e1 = np.zeros((len(Eexact),1))
    e2 = np.zeros((len(Bexact),1))
    for i in range(0, len(e1)):
        e1[i] = E2[i, moment] - Eexact[i]
        e2[i] = B3[i, moment] - Bexact[i] 
    error1 = LA.norm(e1)
    error2 = LA.norm(e2)
    return error1, error2

def TV(E2, B3,y, moment):
    TV1 = 0
    TV2 = 0
    for i in range(0, len(y)-1):
        TV1 = TV1 + abs(E2[i+1, moment]-E2[i, moment])
        TV2 = TV2 + abs(B3[i+1, moment] - B3[i, moment])
    return TV1, TV2

endtime = 100
nt = 7501
nx = 501
dx, dt = stepsizes(nx, nt, endtime)


y, cfl, E2, B3, w1, w2 = highresmaxwell(nx, nt, endtime)
T=100
moment = int(T/dt)
Eexact, Bexact = exactsolution(y, T)
plot_maxwell(E2, B3, y,T, dt,Eexact, Bexact)

error1, error2 = errors(E2, B3, Eexact, Bexact, moment)
TV1, TV2 = TV(E2, B3,y, moment)


print('The error of the electric field is', error1)
print('The error of the magnetic field is', error2)
print('The TV of the electric field solution is', TV1)
print('The TV of the magnetic field solution is', TV1)
print('The cfl number is', cfl)