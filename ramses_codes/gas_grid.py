import numpy as np 
import yt
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
print("Length in units cMpc/h is",L)

cg= ds.covering_grid(level=0, left_edge=[0, 0.0, 0.0], dims=[2**lev, 2**lev, 2**lev])
quan=cg['gas','density']
quan=quan.to('kg/m**3').d
print("Avg. value of Gas Density from simulation is",np.mean(quan), "kg m^(-3)")
quan=quan/np.mean(quan)

overdensity=quan-1

######## Theoretical Value ############
G=6.67e-11
Omb=0.049
print("Omega_b is", Omb)
Mpc2km=3.086e19
Ho=100*h
critical_density=3*(Ho/Mpc2km)**2/(8*np.pi*G)
gas_density=critical_density*Omb*(1+z)**3
print("Theoretical value for gas density is", gas_density, " kg/m^3 \n")

os.makedirs(direc,exist_ok=True)
# Define output filename
output_file_name=direc+"gas_grid_"+info[-4:-1]+".bin"

# Save overdensity as a binary file
with open(output_file_name, "wb") as f:
    overdensity.flatten().astype(np.float64).tofile(f)

print(f"Overdensity field saved as binary file: {output_file_name}")