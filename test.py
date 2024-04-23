import pandas as pd

filepath = '/Users/niwi9751/Desktop/Data Zr In UN/'

# Load the data into a DataFrame, assuming space-separated columns without headers
data = pd.read_csv(file_path, delim_whitespace=True, header=None)

# Select columns of interest based on the user's description
# Column 3: Depth (nm), Column 4: Concentration (atomic fraction)
data_depth_concentration = data[[2, 3]]
data_depth_concentration.columns = ['Depth (nm)', 'Concentration (atomic fraction)']

# Update the atomic density for UN (uranium nitride)
atomic_density_UN = 3.36e22  # atoms/cm³

# Recalculate the total implanted dose (ions/cm²) using the atomic density for UN
total_implanted_dose_UN = (data_depth_concentration['Concentration (atomic fraction)'][:-1] * atomic_density_UN * 
                           data_depth_concentration['Delta Depth (cm)'][1:]).sum()

total_implanted_dose_UN
