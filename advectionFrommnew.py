from numpy import linalg as LA
import numpy as np

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
        return DX,y

def uniformconvectionfromm(nt, nx, tmax, c):
   # for nx give number of nodes on interval [0,1]
   # minimum value allowed is 3
   # Increments
   dt = tmax/(nt-1)
   DX,y = gridfunction(nx)

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

   # Loop
   Q = np.append(Q, zeros, axis=0)
   Q = np.vstack((Q[(len(y)-2)], Q))
   Q = np.vstack((Q, Q[3]))
   for n in range(0,nt-1):
      for i in range(2,len(Q)-1):
         Q[i, n+1] = Q[i,n] - c*dt/(4*DX[i])*(Q[i+1,n] + 3*Q[i,n]- 5*Q[i-1, n]+Q[i-2, n])+c**2*dt**2/(4*(DX[i])**2)*(Q[i+1, n] - Q[i, n] - Q[i-1, n]+Q[i-2, n])
    # Periodic boundary condition
      Q[0, n+1] = Q[len(Q)-3, n+1]
      Q[1, n+1] = Q[len(Q)-2, n+1]
      Q[len(Q)-1, n+1] = Q[3, n+1]
      

   Q = np.delete(Q, [0,1,len(Q)-1], axis=0)
   return y, Q, DX, dt

def plot_convection(u,x,nt,title, f):
   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,100):
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
      
y, Q, DX, dt = uniformconvectionfromm(1001, 101, 10, 1)
f1, f2, f3 = exactsolution(DX[0], y)
plot_convection(Q, y, 1001, 'Advection Equation over time, Fromm Method', f1)
cfl = dt/DX[0]

g1 = np.zeros((len(f1), 1))
g2 = np.zeros((len(f1), 1))
g3 = np.zeros((len(f1), 1))
Q1 = 0
Q2 = 0
Q3 = 0

for i in range(0, len(f1)):
    g1[i] = f1[i] - Q[i, 100]
    g2[i] = f2[i] - Q[i, 400]
    g3[i] = f3[i] - Q[i, 1000]

for i in range(0, len(f1)-1):
    Q1 = Q1 + abs(Q[i+1, 200] - Q[i, 200])
    Q2 = Q2 + abs(Q[i+1, 400] - Q[i, 400])
    Q3 = Q3 + abs(Q[i+1, 1000] - Q[i, 1000])
    
error1 = LA.norm(g1, ord=2)
error2 = LA.norm(g2, ord=2)
error3 = LA.norm(g3, ord=2)
