from numpy import linalg as LA
import numpy as np

def convection(nt, nx, tmax, xmax, c):
   """
   Returns the velocity field and distance for 1D linear convection
   """
   # Increments
   dt = tmax/(nt-1)
   dx = xmax/(nx-1)

   # Initialise data structures
   import numpy as np
   u = np.zeros((nx-1,nt))
   x = np.zeros(nx-1)
         
   # Initial conditions
   for i in range(1,nx-1):
      if(i > (nx-1)/10*2 and i < (nx-1)/10*4):
         u[i,0] = 1
      elif(i > (nx-1)/10*6 and i < (nx-1)/10*8):
          u[i,0] = 1/2 - 1/2*np.cos(np.pi*(i-1/2)*dx)

   # Loop
   for n in range(0,nt-1):
      for i in range(1,nx-1):
         u[i,n+1] = u[i,n]-c*(dt/dx)*(u[i,n]-u[i-1,n])
    # Periodic boundary condition
      u[0, n+1] = u[nx-2, n+1]   

   # X Loop
   for i in range(0,nx-1):
      x[i] = (i+0.5)*dx

   return u, x, dt, dx

def plot_convection(u,x,nt,title, f):
   """
   Plots the 1D velocity field
   """

   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,200):
      plt.plot(x,u[:,i],'r')
      plt.xlabel('x')
      plt.ylabel('q')
      plt.ylim([0,2])
      plt.title(title)
      plt.hold(True)
      plt.plot(x,f, 'b')
      plt.show()
      
def exactsolution(dx):
    import numpy as np
    x = np.arange(0, 10, dx)
    f1 = np.zeros((len(x), 1))
    f2 = np.zeros((len(x), 1))
    f3 = np.zeros((len(x), 1))
    for i in range(0,len(x)):
        if x[i] > 4 and x[i] < 6:
            f1[i] = 1
        if x[i] > 8 and x[i] < 10:
            f1[i] = (1 - np.cos(np.pi*x[i]))/2
    for i in range(0,len(x)):   
        if x[i] > 6 and x[i] < 8:
            f2[i] = 1
        if x[i] > 0 and x[i] < 2:
            f2[i] = (1 - np.cos(np.pi*x[i]))/2    
    for i in range(0,len(x)):   
        if x[i] > 2 and x[i] < 4:
            f3[i] = 1
        if x[i] > 6 and x[i] < 8:
            f3[i] = (1 - np.cos(np.pi*x[i]))/2   
            
    return f1, f2, f3
            
            
u,x, dt, dx = convection(1001, 101, 10, 10, 1)
cfl = dt/dx
f1, f2, f3 = exactsolution(dx)
plot_convection(u, x, 1001, 'Advection Equation over time, upwind method', f3)


g1 = np.zeros((len(f1), 1))
g2 = np.zeros((len(f1), 1))
g3 = np.zeros((len(f1), 1))
Q1 = 0
Q2 = 0
Q3 = 0

for i in range(0, len(f1)):
    g1[i] = f1[i] - u[i, 100]
    g2[i] = f2[i] - u[i, 400]
    g3[i] = f3[i] - u[i, 1000]

for i in range(0, len(f1)-1):
    Q1 = Q1 + abs(u[i+1, 100] - u[i, 100])
    Q2 = Q2 + abs(u[i+1, 400] - u[i, 400])
    Q3 = Q3 + abs(u[i+1, 1000] - u[i, 1000])
    
error1 = LA.norm(g1, ord=2)
error2 = LA.norm(g2, ord=2)
error3 = LA.norm(g3, ord=2)