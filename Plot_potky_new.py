import pandas as pd
import matplotlib.pyplot as plt
import os
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)

# Load your CSV
folder_path = 'Data/ToFERDA/Data/'
files = os.listdir(folder_path)

rho = 6.025# density of target [g/cm^3]
m_a = 123.218# atomic mass of target in [g/mole]
Na = 6.023e23
n_atoms = rho*Na/m_a #atomic density of target [atoms/cm^3]
for file in files:
    df = pd.read_csv(folder_path + file)
    df = df[df['depth'] <= 3000]
    
    # Get x-axis (depth)
    x = df['depth']*3*1e21/n_atoms

    color_dict = {'Zr':'r','O':'b', 'Fe':'gray', 'Xe': 'c', 'Kr':'g', 'Hf': 'y', 'N':'m', 'C':'k', 'Cr': 'silver'}

    # Loop through the columns two at a time (measurement + error)
    columns = df.columns[1:-1]  # skip 'depth'
    plt.figure(figsize=(11, 4))

    for i in range(0, len(columns), 2):
        y_col = columns[i]
        err_col = columns[i+1]
        plt.plot(x, df[y_col]*100, label=y_col, color = color_dict[columns[i]])
        plt.errorbar(x, df[y_col]*100, yerr=df[err_col]*100, color = color_dict[columns[i]],fmt='none')
    if 'Xe' in file:
        plt.axvline(75, color='black', linestyle='--',label= 'SRIM Peak')
        plt.title(r'\textbf{Xe}')
    elif 'Kr' in file:
        plt.axvline(105,color='black', linestyle='--',label= 'SRIM Peak')
        plt.title(r'\textbf{Kr}')
    plt.xlabel('Depth [nm]')
    plt.ylabel(r'Concentration [at.\%]')
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.ylim((0.1,100))
    plt.grid(True)
    plt.tight_layout()

plt.show()
