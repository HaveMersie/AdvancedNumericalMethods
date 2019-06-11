import numpy as np
from numpy import linalg as LA

def gridfunction(nx):
    # the input nx should be the number of nodes on the interval [0,1]
    # minimum value of nx is 3
        import numpy as np
        x_1 = np.linspace(0,5,nx+4*(nx-1))
        x_2 = np.linspace(5,10, 5*(4+2*(nx-3))+1)
        x_2 = np.delete(x_2, 0,0)
        x = np.append(x_1, x_2)
        y = np.zeros(len(x)-1)
        DX = np.zeros(len(x)-1)
        for i in range(0, len(x)-1):
            y[i]=(x[i+1]+x[i])/2
        for i in range(0, len(x)-1):
            DX[i]=(x[i+1]-x[i])/2  
        DX = np.append(DX[0],(np.append(DX[0], DX)))
        return DX,y

def nonuniformconvectionfromm(nt, nx, tmax, c):
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

def plot_convection(Q,x,nt,title):
   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,100):
      plt.plot(x,Q[:,i],'r')
      plt.xlabel('x (m)')
      plt.ylabel('q (m/s)')
      plt.ylim([0,2.2])
      plt.title(title)
      plt.hold(True)
      plt.plot(x, exactsolution3(y))
      plt.show()
      
def exactsolution1(y):
    # exact solution after 1 period
    import numpy as np
    f = np.zeros((len(y), 1))
    for i in range(0, len(y)):
        if y[i]>4 and y[i]<6:
            f[i]=1
        if y[i] >8 and y[i]<10:
            f[i] = 0.5-0.5*np.cos(np.pi*y[i])
    return f

def exactsolution2(y):
    # exact solution after 2 period
    import numpy as np
    f = np.zeros((len(y), 1))
    for i in range(0, len(y)):
        if y[i]>6 and y[i]<8:
            f[i]=1
        if y[i] >0 and y[i]<2:
            f[i] = 0.5-0.5*np.cos(np.pi*y[i])
    return f

def exactsolution3(y):
    # exact solution after 5 period
    import numpy as np
    f = np.zeros((len(y), 1))
    for i in range(0, len(y)):
        if y[i]>2 and y[i]<4:
            f[i]=1
        if y[i] >6 and y[i]<8:
            f[i] = 0.5-0.5*np.cos(np.pi*y[i])
    return f
      
y, Q, DX, dt = nonuniformconvectionfromm(1001, 25, 10, 1)
plot_convection(Q, y, 1001, 'Advection Equation over time with non-uniform grid, Fromm')
cfl2 = dt/DX[len(DX)-1]

exact = exactsolution3(y)
error = np.zeros((len(exact), 1))

for i in range(0, len(exact)):
    error[i] = exact[i]-Q[i,500]


errornorm = LA.norm(error, ord=2)
P = 0
for i in range(0, len(exact)-1):
    P = P + abs(Q[i+1, 500] - Q[i, 500])