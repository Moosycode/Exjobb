import numpy as np

# Your atomic fraction data
atomic_fraction = np.array([
    2.61750265e-03, 2.61750265e-03, 3.38660623e-03, 4.29431902e-03,
    5.03211049e-03, 5.40650986e-03, 5.37971547e-03, 5.01723524e-03,
    4.44077461e-03, 3.77465759e-03, 3.09383531e-03, 2.43225549e-03,
    1.82329115e-03, 1.28419184e-03, 8.13330683e-04, 4.46595741e-04,
    2.32256978e-04, 1.44705664e-04, 6.10703962e-05, 1.43093528e-05,
    1.99133600e-06, 1.82928802e-07, 1.22001442e-08, 6.28596853e-10,
    2.68697986e-11, 9.96906560e-13, 3.27647541e-14, 9.66464656e-16,
    2.58548991e-17, 6.32724645e-19, 1.42667706e-20, 2.98212691e-22,
    5.80889287e-24, 1.05953791e-25, 1.81832623e-27, 0.00000000e+00
])

# Constants
atomic_density_zro2 = 8.31e22  # atoms/cm³
depth_step_cm = 1.6857142857142858e-06 # 10 nm = 1e-6 cm

# Convert to atoms/cm³
kr_concentration_cm3 = atomic_fraction * atomic_density_zro2

# Integrate over depth to get atoms/cm²
areal_density = np.sum(kr_concentration_cm3 * depth_step_cm)

print(f"Areal density of Kr = {areal_density:.2e} atoms/cm²")
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd

# ----------------- Physical Constants -----------------
k = 8.6e-5         # Boltzmann constant [eV/K]
Na = 6.022e23      # Avogadro's number [atoms/mol]

# ----------------- Utility Functions ------------------

def D(D0, Ea, T):
    """Temperature-dependent diffusion coefficient."""
    return D0 * np.exp(-Ea / (k * T))

def hist_integral(n, width):
    """Integral of histogram (sum of counts * bin width)."""
    return np.sum(n) * width

def V(x):
    """Optional velocity field (placeholder)."""
    return a

# ------------------ Data Classes ----------------------

class SDTrimSP_output:
    """Parser for SDTrimSP output file (E0_31_target.dat)."""

    def __init__(self, filename):
        d_marker = 0
        h_marker = 0
        d_data = []
        s_data = []

        for line in reversed(list(open(filename))):
            if d_marker == 0 and line[0] == ' ':
                d_marker = 1
            elif d_marker == 1 and line[0] != ' ':
                d_marker = 0
                h_marker += 1
            if h_marker >= 1:
                break
            if d_marker == 0:
                d_data.append([float(i) for i in line.rstrip().split()])
            else:
                s_data.append(line.rstrip().split())

        d_data = np.flipud(np.array(d_data)).T
        self.s = s_data

        # Parse scalar properties
        self.nhist   = float(s_data[-1][0])
        self.srem    = float(s_data[-1][1]) / 10
        self.areald  = [1e16 * float(i) for i in s_data[-4][:2]]
        self.npart   = float(s_data[-5][0])
        self.ibackp  = [float(i) for i in s_data[-6][:2]]
        self.ibackr  = [float(i) for i in s_data[-10][:2]]

        # Depth-resolved data
        self.depth   = d_data[0] / 10
        self.density = d_data[1] * 1e24
        self.afracs  = d_data[2:-2]
        self.Pdam    = d_data[-2]
        self.Pdam_dn = self.Pdam / np.sum(self.Pdam)
        self.binw_nm = self.depth[1] - self.depth[0]
        self.binw_cm = self.binw_nm * 1e-7
        self.something = d_data[-1]  # Unknown last column

class data_table:
    """Simple parser for column-based table data."""

    def __init__(self, filename, delimiter='\t'):
        with open(filename, 'r') as file:
            self.headers = file.readline().strip().split(delimiter)
        self.ncols = len(self.headers)
        self.data = np.loadtxt(filename, delimiter=delimiter, skiprows=1).T

    def get_data(self, header):
        return self.data[self.headers.index(header)]

def get_trim_res(path):
    """Load and report SDTrimSP results."""
    d = SDTrimSP_output(path + 'E0_31_target.dat')
    print('-------------------------------------------\n')
    print(f'number of histories: {int(d.nhist)}')
    print(f'number of particles: {int(d.npart)}')
    print(f'surface removal:     {np.round(d.srem, 1)} nm')
    print(f'areal density of Kr: {d.areald[0]:.2e} atoms/cm^2')
    print(f'total area of Pdam:  {np.sum(d.Pdam)}\n')
    return d

# ---------------- Element Database --------------------

elementdict = {
    'Fe_ZrO2': {'D0': 2.26e-6, 'Ea': 2.3,  'rho': 6.08, 'Ma': 123.218},
    'Kr_ZrO2': {'D0': 8.11e-7, 'Ea': 2.53, 'rho': 6.08, 'Ma': 123.218},
    'Xe_ZrO2': {'D0': 1.83e-6, 'Ea': 2.91, 'rho': 6.08, 'Ma': 123.218},
    'Zr_UN':   {'D0': 6.9e-7,  'Ea': 2.7,  'rho': 14.05,'Ma': 252.036},
    'Kr_UN':   {'D0': 8.11e-7, 'Ea': 2.53, 'rho': 14.05,'Ma': 252.036},
    'Xe_UN':   {'D0': 1.83e-6, 'Ea': 2.91, 'rho': 14.05,'Ma': 252.036},
    'Zr_UO2':  {'D0': 6.9e-7,  'Ea': 2.7,  'rho': 10.6, 'Ma': 270.02},
    'Kr_UO2':  {'D0': 8.11e-7, 'Ea': 2.53, 'rho': 10.6, 'Ma': 270.02}
}

# ---------------- Simulation Setup --------------------

# File paths
pre_data_path = 'Data/ToFERDA/Data/pre-anneal_Kr-imp_corrected.csv'
post_data_path = 'Data/ToFERDA/Data/post-anneal_Kr-imp_corrected.csv'
d = get_trim_res('Data/SDTrim/')
interp_func = interp1d(d.depth, d.Pdam, kind='linear', fill_value="extrapolate")

# Experiment conditions
Times_in = [9, 9]  # in hours
Temp = 1473.15     # in Kelvin
fluence = 1e17
element = 'Kr_ZrO2'
sample = 'Kr'
Times = [t * 3600 for t in Times_in]  # to seconds
Concentrations = []

# Plot settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

# ---------------- Main Diffusion Loop -----------------

for T in Times:
    # --- Data prep ---
    df_pre = pd.read_csv(pre_data_path)
    df_pre = df_pre[(df_pre['depth'] >= 0) & (df_pre['depth'] <= 3000)]
    df_post = pd.read_csv(post_data_path)
    df_post = df_post[(df_post['depth'] >= 0) & (df_post['depth'] <= 3000)]

    c_pre = df_pre[sample].values
    c_post = df_post[sample].values
    x_pre = df_pre['depth'].values
    x_post = df_post['depth'].values
    binwidth = (x_pre[1] - x_pre[0]) * 1e15

    I_initial = hist_integral(c_pre, binwidth)
    I_diffused = hist_integral(c_post, binwidth)
    print(f'Initial: {I_initial:.2e}, Diffused: {I_diffused:.2e}, Ratio: {I_initial/I_diffused:.2f}')

    # --- Physical Parameters ---
    D0 = elementdict[element]['D0'] * 1e8
    Ea = elementdict[element]['Ea']
    a = 1e-5
    d = 1
    rho = elementdict[element]['rho']
    Ma = elementdict[element]['Ma']
    n_atoms = rho * Na / Ma
    Nx = len(x_pre)
    Nt = int(T / 60) * 2
    x_um = 3e21 * x_pre / n_atoms
    Extendby = 2
    L = int(x_um[-1]) * Extendby
    dx = L / (L * Nx - 1)
    dt = T / Nt

    x = np.linspace(0, L, Nx * Extendby)
    C = np.zeros((Nt, Nx * Extendby))
    C[0, :Nx] = c_pre

    # Interpolated damage profile and diffusivity
    rebinned = interp_func(x)
    rebinned = (rebinned * d / np.sum(rebinned)) + 1
    D_x = D(D0, Ea, Temp) * rebinned
    V_x = np.array([V(i * dx) for i in range(Nx * Extendby)])

    # --- Diffusion Loop ---
    for n in range(Nt - 1):
        Cn = C[n]
        diffusion_term = (D_x[2:] * (Cn[2:] - Cn[1:-1]) - D_x[:-2] * (Cn[1:-1] - Cn[:-2])) / dx**2
        C[n+1, 1:-1] = Cn[1:-1] + dt * diffusion_term
        C[n+1] -= C[n] * V_x * dt
        C[n+1, -1] = 0
        C[n+1, 0] = C[n+1, 1]

    Concentrations.append(C[-1])

# ---------------- Plotting ----------------------------

plt.figure(figsize=(8, 6))
bin_w = (x[1] - x[0]) * 1e-7
postint = hist_integral(Concentrations[-1][:21], bin_w)
print(f'Modeled integral: {postint * n_atoms * 3:.2e} atoms/cm^2')

scale = 100
plt.errorbar(x_pre, c_pre * scale, yerr=df_pre[f'{sample}_sd'] * scale,
             label='As implanted', capsize=3, color='blue')
plt.errorbar(x_post, c_post * scale, yerr=df_post[f'{sample}_sd'] * scale,
             label='Post-anneal', capsize=3, color='orange')

for i, C_ in enumerate(Concentrations):
    plt.plot(x, C_ * scale, linestyle='--',
             label=f'{i+1}$^{{th}}$ model', color='green' if i else 'red')

plt.xlabel('Depth [nm]')
plt.ylabel('Concentration [at.\%]')
plt.title(r'\textbf{Kr}')
plt.grid(True)
plt.xlim([0, 300])
plt.ylim([0, 5.5])
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(x, D_x)
plt.title('Spatially Resolved Diffusivity')
plt.xlabel('Depth [nm]')
plt.ylabel(r'D(x) [$\mu$ m$^2$/s]')
plt.grid(True)
plt.tight_layout()

plt.show()
