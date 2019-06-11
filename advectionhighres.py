from numpy import linalg as LA
import numpy as np

def minmod(a,b):
    global r
    if abs(a) < abs(b) and (a*b)>0:
        r = a
    if abs(b)<abs(a) and (a*b)>0:
        r = b
    if (a*b)<=0:
        r = 0
    return r

def gridfunction(nx):
    # the input nx should be the number of nodes on the interval [0,1]
    # minimum value of nx is 3
        import numpy as np
        x = np.linspace(0,10,nx+9*(nx-1))
        y = np.zeros(len(x)-1)
        DX = np.zeros(len(x)-1)
        for i in range(0, len(x)-1):
            y[i]=(x[i+1]+x[i])/2
        for i in range(0, len(x)-1):
            DX[i]=(x[i+1]-x[i])/2  
        DX = np.append(DX[0],(np.append(DX[0], DX)))
        return y

def uniformadvectionhighres(nt, nx, tmax, c):
   # for nx give number of nodes on interval [0,1]
   # minimum value allowed is 3
   # Increments
   dt = tmax/(nt-1)
   y = gridfunction(nx)
   DX = 10/(nx+9*(nx-1)-1)

   # Initialise data structures
   import numpy as np
   zeros = np.zeros((1, nt))
   Q = np.zeros((len(y),nt))
         
   # Initial conditions
   for i in range(0,len(y)):
      if y[i] > 2 and y[i]<4:
         Q[i,0] = 1
      elif y[i] > 6 and y[i] < 8:
          Q[i,0] = 1/2 - 1/2*np.cos(np.pi*(y[i]))
   
   Q = np.vstack((Q, Q[0]))
   Q = np.vstack((Q[len(Q)-2], Q))
   Q = np.vstack((Q[len(Q)-3], Q))
   
   # Loop
   for n in range(0,nt-1):
      for i in range(2,len(Q)-1):
         Q[i, n+1] = Q[i,n] - c*dt/(DX)*(Q[i,n]-Q[i-1,n]) - 0.5*c*dt/(DX)*(DX-c*dt)*(minmod((Q[i, n]-Q[i-1, n])/DX,(Q[i+1, n]-Q[i, n])/DX) - minmod((Q[i-1, n]-Q[i-2, n])/DX,(Q[i, n]-Q[i-1, n])/DX))
    # Periodic boundary condition
      Q[0, n+1] = Q[len(Q)-3, n+1]
      Q[1, n+1] = Q[len(Q)-2, n+1]
      Q[len(Q)-1, n+1] = Q[2, n+1]
      

   Q = np.delete(Q, [0,1,len(Q)-1], axis=0)
   return y, Q, DX, dt

def plot_convection(u,x,nt,title, f):
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
      

      
def exactsolution(dx, y):
    import numpy as np
    f1 = np.zeros((len(y), 1))
    f2 = np.zeros((len(y), 1))
    f3 = np.zeros((len(y), 1))
    for i in range(0,len(y)):
        if y[i] > 4 and y[i] < 6:
            f1[i] = 1
        if y[i] > 8 and y[i] < 10:
            f1[i] = (1 - np.cos(np.pi*y[i]))/2
    for i in range(0,len(y)):   
        if y[i] > 6 and y[i] < 8:
            f2[i] = 1
        if y[i] > 0 and y[i] < 2:
            f2[i] = (1 - np.cos(np.pi*y[i]))/2    
    for i in range(0,len(y)):   
        if y[i] > 2 and y[i] < 4:
            f3[i] = 1
        if y[i] > 6 and y[i] < 8:
            f3[i] = (1 - np.cos(np.pi*y[i]))/2   
            
    return f1, f2, f3
      
y, Q, DX, dt = uniformadvectionhighres(1001, 11, 10, 1)
f1, f2, f3 = exactsolution(DX, y)
plot_convection(Q, y, 1001, 'Advection Equation over time, High Resolution', f3)
cfl = dt/DX

g1 = np.zeros((len(f1), 1))
g2 = np.zeros((len(f1), 1))
g3 = np.zeros((len(f1), 1))
Q1 = 0
Q2 = 0
Q3 = 0

for i in range(0, len(f1)):
    g1[i] = f1[i] - Q[i, 200]
    g2[i] = f2[i] - Q[i, 400]
    g3[i] = f3[i] - Q[i, 1000]

for i in range(0, len(f1)-1):
    Q1 = Q1 + abs(Q[i+1, 200] - Q[i, 200])
    Q2 = Q2 + abs(Q[i+1, 400] - Q[i, 400])
    Q3 = Q3 + abs(Q[i+1, 1000] - Q[i, 1000])
    
error1 = LA.norm(g1, ord=2)
error2 = LA.norm(g2, ord=2)
error3 = LA.norm(g3, ord=2)