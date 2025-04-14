# Define the function to compute T2 given G2/G1 ratio
import matplotlib.pyplot as plt
import numpy as np


plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=14)

def func(T1, k, Emv, Efv, ratio):
    numerator = k * T1**2/(Emv + 2*Efv) * np.log(ratio)
    denominator = 1 - k * T1/(Emv + 2 * Efv)* np.log(ratio)
    return T1 + numerator / denominator

# Fixed parameters
ratios = np.linspace(1.0, 1000, 500)  # Avoid ratio=1 to prevent log(1)=0 causing issues
k = 8.6e-5
Emv = 1.5
Efv = 2

Ts = [723.15,773.15, 823.15]  # Fixed T1 in Kelvin
plt.figure(figsize=(8, 6))

G1 = 1e-6 #Dpa/s
G2 = 1e-4 #Dpa/s
highlight_ratio = 100
for T1 in Ts:
    T2 = func(T1, k, Emv, Efv, ratios)
    highlight_T2 = func(T1, k, Emv, Efv, highlight_ratio)
    # Plotting
    plt.plot( ratios, T2-273.15, label=r'T1 = ' + str(T1-273.15)+ r'\textdegree C')
    plt.plot(highlight_ratio,highlight_T2 - 273.15, 'ro')

plt.plot()
plt.ylabel(r'T$_2$')
plt.xlabel(r'G$_2$/G$_1$')
plt.xscale('log')
plt.title(r'T$_2$ as a Function of G$_2$/G$_1$ Ratio')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
