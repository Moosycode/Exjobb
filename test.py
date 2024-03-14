import numpy as np
import matplotlib.pyplot as plt

root = '/Users/nilsw/Python/Exjobb/verticalfurnaceupto1600.txt'
data = np.loadtxt(root,usecols=1,unpack=True, dtype='U25',skiprows=1 )


x_vals = np.linspace(0,len(data),len(data))
print(x_vals)