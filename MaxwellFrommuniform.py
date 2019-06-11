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

def uniformmaxwellfromm(nt, nx, tmax):
   # for nx give number of nodes on interval [0,1]
   # minimum value allowed is 3
   # Increments
   dt = tmax/(nt-1)
   DX,y = gridfunction(nx)

   # Initialise data structures
   import numpy as np
   zeros = np.zeros((1, nt))
   Q1 = 0.5+np.ones((len(y),nt))
   Q2 = 0.5+np.ones((len(y),nt))
         
   # Initial conditions
   for i in range(0,len(y)):
      if y[i] > 2 and y[i]<4:
         Q1[i,0] = 1+0.5
      elif y[i] > 6 and y[i] < 8:
          Q1[i,0] = 0.5 - 1/2*np.cos(np.pi*(y[i])) + 0.5
          
   for i in range(0,len(y)):
      if y[i] > 2 and y[i]<4:
         Q2[i,0] = -(1) + 0.5
      elif y[i] > 6 and y[i] < 8:
          Q2[i,0] = -(0.5 - 1/2*np.cos(np.pi*(y[i]))) + 0.5

   # Loop
   Q1 = np.vstack((Q1[(len(Q1)-1)], Q1))
   Q1 = np.vstack((Q1[(len(Q1)-2)], Q1))
   Q1 = np.vstack((Q1, Q1[2]))
   for n in range(0,nt-1):
      for i in range(2,len(Q1)-1):
         Q1[i, n+1] = Q1[i,n] - 0.5*dt/(4*DX[i])*(Q1[i+1,n] + 3*Q1[i,n]- 5*Q1[i-1, n]+Q1[i-2, n])+0.5**2*dt**2/(4*(DX[i])**2)*(Q1[i+1, n] - Q1[i, n] - Q1[i-1, n]+Q1[i-2, n])
    # Periodic boundary condition
      Q1[0, n+1] = Q1[len(Q1)-2, n+1]
      Q1[1, n+1] = Q1[len(Q1)-1, n+1]
      Q1[len(Q1)-1, n+1] = Q1[2, n+1]

      
   Q2 = np.vstack((Q2[(len(Q2)-1)], Q2))
   Q2 = np.vstack((Q2, Q2[1]))
   Q2 = np.vstack((Q2, Q2[2]))
   for n in range(0,nt-1):
      for i in range(1,len(Q2)-3):
         Q2[i, n+1] = Q2[i,n] - 0.5*dt/(4*DX[i])*(Q2[i-1,n] + 3*Q2[i,n]- 5*Q2[i+1, n]+Q2[i+2, n])+0.5**2*dt**2/(4*(DX[i])**2)*(Q2[i-1, n] - Q2[i, n] - Q2[i+1, n]+Q2[i+2, n])
    # Periodic boundary condition
      Q2[0, n+1] = Q2[len(Q2)-3, n+1]
      Q2[len(Q1)-1, n+1] = Q2[2, n+1]
      Q2[len(Q1)-2, n+1] = Q2[1, n+1]
      

   Q1 = np.delete(Q1, [0,1,len(Q1)-1], axis=0)
   Q2 = np.delete(Q2, [0,len(Q2)-1, len(Q2)-2], axis=0)
   return y, Q1, Q2

def plot_convection(Q,x,nt,title):
   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,10):
      plt.plot(x,Q[:,i],'r')
      plt.xlabel('x (m)')
      plt.ylabel('q (m/s)')
      plt.ylim([0,5])
      plt.title(title)
      plt.show()
      
y, Q1, Q2 = uniformmaxwellfromm(1001, 10, 100)
E2=0.5*(Q1-Q2)
B3 = Q1+Q2
plot_convection(E2, y, 1001, 'Electric field over time with uniform grid')
plot_convection(B3, y, 1001, 'Magnetic field over time with uniform grid')