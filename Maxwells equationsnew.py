   # Initialise data structures
import numpy as np
f = np.zeros(nx)
x = np.zeros(nx)
w1 = np.zeros((nx, nt))
w2 = np.zeros((nx, nt))
   
    for i in range(0, nx):
       x[i] = (i+0.5)*dx

            
   # Creating function f
   for i in range(1,nx-1):
      if(i > (nx-1)/10*2 and i < (nx-1)/10*4):
         f[i] = 1
      elif(i > (nx-1)/10*6 and i < (nx-1)/10*8):
          f[i] = 1/2 - 1/2*np.cos(np.pi*(i-1/2)*dx)
          
   # Initial conditions
   w1[:,0] = np.ones(nx)*0.5 + f
   w2[:,0] = np.ones(nx)*0.5 - f

   # Loop
   for n in range(0,nt-1):
      for i in range(1,nx-1):
         w1[i,n+1] = w1[i,n]-0.5*(dt/dx)*(w1[i,n]-w1[i-1,n])
         w2[i,n+1] = w2[i,n]+0.5*(dt/dx)*(w2[i,n]-w2[i-1,n])
    # Periodic boundary condition
      w1[0, n+1] = w1[nx-1, n+1]   
      w2[0, n+1] = w2[nx-1, n+1]    
   return w1, w2, x

def plot_convection(u,x,nt,title):
   """
   Plots the 1D field
   """

   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,400):
      plt.plot(x,u[:,i],'r')
      plt.xlabel('x')
      plt.ylabel('q')
      plt.ylim([0,2])
      plt.title(title)
      plt.show()
      
w1, w2, x = Maxwell(10001, 1001, 50, 10)
E2 = 0.5*(w1 - w2)
B3 = w1 + w2
plot_convection(B3, x, 10001, 'Electric field over time')
plot_convection(E2, x, 10001, 'Magnetic field over time')