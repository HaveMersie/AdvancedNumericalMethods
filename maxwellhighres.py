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
        return y

def uniformmaxwellhighres(nt, nx, tmax):
   # for nx give number of nodes on interval [0,1]
   # minimum value allowed is 3
   # Increments
   dt = tmax/(nt-1)
   y = gridfunction(nx)
   DX = 10/(nx+9*(nx-1)-1)

   # Initialise data structures
   import numpy as np
   zeros = np.zeros((1, nt))
   Q1 = 0.5*np.ones((len(y),nt))
   Q2 = 0.5*np.ones((len(y),nt))
         
   # Initial conditions
   for i in range(0,len(y)):
      if y[i] > 2 and y[i]<4:
         Q1[i,0] = 1.5
      elif y[i] > 6 and y[i] < 8:
          Q1[i,0] = 1 - 0.5*np.cos(np.pi*(y[i]))
          
   for i in range(0,len(y)):
      if y[i] > 2 and y[i]<4:
         Q2[i,0] = -0.5
      elif y[i] > 6 and y[i] < 8:
          Q2[i,0] = 0.5*np.cos(np.pi*(y[i]))

   # Loop
   Q1 = np.vstack((Q1[(len(Q1)-1)], Q1))
   Q1 = np.vstack((Q1[(len(Q1)-2)], Q1))
   Q1 = np.vstack((Q1, Q1[2]))
   for n in range(0,nt-1):
      for i in range(2,len(Q1)-1):
         Q1[i, n+1] = Q1[i,n] - 0.5*dt/(DX)*(Q1[i,n]-Q1[i-1,n]) - 0.5*0.5*dt/(DX)*(DX-0.5*dt)*(minmod((Q1[i, n]-Q1[i-1, n])/DX,(Q1[i+1, n]-Q1[i, n])/DX) - minmod((Q1[i-1, n]-Q1[i-2, n])/DX,(Q1[i, n]-Q1[i-1, n])/DX))
    # Periodic boundary condition
      Q1[0, n+1] = Q1[len(Q1)-2, n+1]
      Q1[1, n+1] = Q1[len(Q1)-1, n+1]
      Q1[len(Q1)-1, n+1] = Q1[2, n+1]

      
   Q2 = np.vstack((Q2[(len(Q2)-1)], Q2))
   Q2 = np.vstack((Q2, Q2[1]))
   Q2 = np.vstack((Q2, Q2[2]))
   for n in range(0,nt-1):
      for i in range(1,len(Q2)-2):
          Q2[i, n+1] = Q2[i,n] - 0.5*dt/(DX)*(Q2[i,n]-Q2[i+1,n]) - 0.5*0.5*dt/(DX)*(DX-0.5*dt)*(minmod((Q2[i, n]-Q2[i+1, n])/DX,(Q2[i-1, n]-Q2[i, n])/DX) - minmod((Q2[i+1, n]-Q2[i+2, n])/DX,(Q2[i, n]-Q2[i+1, n])/DX))
    # Periodic boundary condition
      Q2[0, n+1] = Q2[len(Q2)-3, n+1]
      Q2[len(Q1)-1, n+1] = Q2[2, n+1]
      Q2[len(Q1)-2, n+1] = Q2[1, n+1]
      

   Q1 = np.delete(Q1, [0,1,len(Q1)-1], axis=0)
   Q2 = np.delete(Q2, [0,len(Q2)-1, len(Q2)-2], axis=0)
   return y, Q1, Q2, dt

def plot_convection(Q,x,nt,title, exact):
   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,100):
      plt.plot(x,Q[:,i],'r')
      plt.xlabel('x')
      plt.ylabel('q')
      plt.ylim([0,4])
      plt.title(title)
      plt.hold(True)
      plt.plot(x,exact,'b')
      plt.show()
      
def minmod(a,b):
    global k
    if abs(a) < abs(b) and (a*b)>0:
        k = a
    elif abs(b)<abs(a) and (a*b)>0:
        k = b
    elif (a*b)<=0:
        k = 0
    return k

def exactsolution1(y):
    #after 1 period
    Eexact = np.zeros((len(y), 1))
    Bexact = np.ones((len(y), 1))
    for i in range(0, len(y)):
        if y[i] > 1 and y[i]<3:
            Eexact[i] = 1/2
            Bexact[i] = Bexact[i] - 1
        if y[i] > 3 and y[i]<5:
            Eexact[i] = 1/2
            Bexact[i] = Bexact[i] + 1
        if y[i] > 5 and y[i] < 7:
            Eexact[i] = (1 - np.cos(np.pi*y[i]))/4
            Bexact[i] = Bexact[i] + (np.cos(np.pi *y[i])-1)/2
        if y[i] > 7 and y[i] < 9:
            Eexact[i] = (1 - np.cos(np.pi*y[i]))/4    
            Bexact[i] = Bexact[i] - (np.cos(np.pi *y[i])-1)/2
    return Eexact, Bexact
            
y, Q1, Q2, dt = uniformmaxwellhighres(10001,40 , 50)
E2=0.5*(Q1-Q2)
B3 = Q1+Q2
Eexact, Bexact = exactsolution1(y)
#plot_convection(E2, y, 10001, 'Electric field over time with uniform grid', Eexact)
plot_convection(B3, y, 10001, 'Magnetic field over time with uniform grid', Bexact)

