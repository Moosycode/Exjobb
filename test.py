import numpy as np
import matplotlib.pyplot as plt

L = 1
Lprime = 4 
Nx = 5 # Number of spatial points
Nt = 4

x = np.linspace(0, Lprime, 4*Nx)

# Initialize solution matrix
C = np.zeros((Nt, Nx))
print(C)
# Apply initial condition

C[0, :] = [1, 2 ,3, 4, 5]
print(C)

newC = np.hstack((C, np.zeros((Nt,(-L+Lprime)*Nx))))
print(newC)

plt.plot(x,newC[0,:])
plt.show()