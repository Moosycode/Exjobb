import numpy as np
import datetime
import matplotlib.pyplot as plt

def histog(data, length):
    L = length*1000000 #Make sure unit is correct
    n, bins = np.histogram(data, bins = 100,range = (0,L))
    width = bins[1]-bins[0]
    # x = np.linspace(min(data), L,100)
    # x = x*0.0001 #rescale to micrometer, if file is in Å
    y = n/sum(n) #Normalize
    return width,y


root = '/Users/nilsw/Srim Results/Fe_in_ZrO2_300keV.txt'
data = np.loadtxt(root,skiprows=18,usecols=1,unpack=True,encoding='cp437',dtype= 'U25')
data = [d.replace(",","") for d in data]
data = [float(d) for d in data]
data = [int(d) for d in data]

binwidth,y_hist = histog(data, 1)

plt.hist(data, bins=100, edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram')
plt.grid(True)
plt.show()