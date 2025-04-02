import numpy as np 
import yt
from CIC_package import cloud_in_cell_3D
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="path of info_*.txt file")
parser.add_argument("-lm","--levelmin",help="Base level of simulation",default=8)
parser.add_argument("-o","--output", help="folder for outputs")

args = parser.parse_args()
info = args.input
lev=int(args.levelmin)
direc=args.output

N = 2**lev

ds=yt.load(info)
h = ds.hubble_constant
print("Hubble Constant is",h)
z = ds.current_redshift
rounded_redshift=round(z,2)
L=ds.length_unit.in_units('Mpccm').d*h  #length in cMpc/h
print("Lemgth in units cMpc/h is",L)

# --- Dark Matter Deposition using CICDeposit_3 ---
# Extract dark matter particle positions and masses, converting to Mpc/h and Msol.
l=ds.length_unit.in_units('m').d	#length in m
dm_mass = ds.all_data()[("all", "particle_mass")].in_units("kg").to_ndarray()
dm_px = ds.all_data()[("all", "particle_position_x")].in_units("m").to_ndarray()
dm_py = ds.all_data()[("all", "particle_position_y")].in_units("m").to_ndarray()
dm_pz = ds.all_data()[("all", "particle_position_z")].in_units("m").to_ndarray()
npositions = dm_px.shape[0]  # total number of particles
print("Total number of particles is", npositions)

# Prepare an empty grid (field) to deposit DM mass.
dm_grid = np.zeros((N, N, N), dtype=np.float64)

# Call cloud_in_cell_3D:
dm_grid=cloud_in_cell_3D(dm_px,dm_py,dm_pz, dm_mass, N, l)

avg_density=np.mean(dm_grid)

print("Avg. value of Dark Matter Density from simulation is",avg_density, "kg m^(-3)")

overdensity=(dm_grid/avg_density)-1

###### Theoretical Overdensity is ###############
G=6.67e-11
Om_m = ds.omega_matter
print("Omega_matter is", Om_m)
Mpc2km=3.086e19
Ho=100*h
critical_density=3*(Ho/Mpc2km)**2/(8*np.pi*G)
dark_matter_density=critical_density*Om_m*(1+z)**3
print("Theoretical value for dark matter density is", dark_matter_density, " kg/m^3 \n")

os.makedirs(direc,exist_ok=True)
# Define output filename
output_file_name=direc+"dm_grid_"+info[-4:-1]+".bin"

# Save overdensity as a binary file
with open(output_file_name, "wb") as f:
    overdensity.flatten().astype(np.float64).tofile(f)

print(f"Overdensity field saved as binary file: {output_file_name}")



