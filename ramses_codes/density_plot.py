import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from matplotlib.colors import LogNorm
import yt
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter

parser = argparse.ArgumentParser(description="Plot 1+overdensity field slices and projections.")
parser.add_argument("-lm", "--levelmin", type=int, required=True, help="Base level (example: 8 for 256^3).")
parser.add_argument("-len", "--length", type=int, required=True, help="Box length in cMpc/h.")
parser.add_argument("-t", "--type", choices=["dm", "gas"], default="dm", help="Type of field: 'dm' or 'gas'.")
parser.add_argument("-haloes", "--haloes", type=int, choices=[0, 1], default=0, help="0: No haloes, 1: Add haloes.")
parser.add_argument("-n", "--number", type=int, required=True, help="Snapshot number, e.g., '004' or '035'.")
args = parser.parse_args()

# Read arguments
levelmin = args.levelmin
box_length = args.length
field_type = args.type
add_haloes = args.haloes
snapshot_number = args.number

z_red=np.loadtxt(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/redshift.txt")
rounded_redshift=round(z_red[snapshot_number-1],2)

N = 2**levelmin
print(f"Grid size: {N} x {N} x {N}")
print(f"Box length: {box_length} cMpc/h")
print(f"Field type: {field_type}")
print(f"Add haloes overlay: {'Yes' if add_haloes==1 else 'No'}")

# Location of grid file
if field_type=="dm":
    bin_file=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/DM_grid/dm_grid_{snapshot_number:03}.bin"
elif field_type=="gas":
    bin_file=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/Gas_grid/gas_grid_{snapshot_number:03}.bin"
else:
    print("Invalid Field Type")
    exit()
# Load overdensity field
print("Reading file : ",bin_file)
overdensity = np.fromfile(bin_file, dtype=np.float64)
overdensity = overdensity.reshape((N, N, N))

# 1 + overdensity
one_plus_delta = 1.0 + overdensity

# Middle slice
z_index = N // 2
slice_data = one_plus_delta[:, :, z_index]

# Projection along z-axis
projection_data = np.mean(one_plus_delta, axis=2)

# --- Plot Projection ---
fig, ax = plt.subplots(figsize=(8, 6))
c = ax.imshow(projection_data.T, cmap="viridis", aspect='equal', origin='lower',
              extent=[-box_length/2, box_length/2, -box_length/2, box_length/2], norm=colors.LogNorm(vmin=np.min(projection_data), vmax=np.max(projection_data)))

# ---- Add Haloes if requested ----
if add_haloes == 1:
    print("Adding Haloes on Projection Plot")
    hp = f"/user1/poojarani/Lya_Comparison/ramses_data/lev{levelmin:03}_len{box_length}/Halo_data/info_{snapshot_number:05}/info_{snapshot_number:05}.0.h5"

    # Check if halo catalog exists
    if not os.path.exists(hp):
        print('Pooja: Halo catalog does not exist! Please generate it first. Bye!')
        sys.exit()

    dsh = yt.load(hp)
    adh = dsh.all_data()
    # Check if there are any haloes
    mass = adh['halos', 'particle_mass'].d
    if mass.shape[0] == 0:
        print("No haloes found in the catalog. Skipping plot.")
        sys.exit()
    else:
        print(f"Found {mass.shape[0]} haloes.")
        posi = adh['halos', 'particle_position'].in_units('Mpccm/h').d - box_length/2
        m_min = min(mass)

        # Plot haloes on the graph
        ax.scatter(posi[:, 0], posi[:, 1], facecolors='none', edgecolors='white', s=20*(mass/m_min), linewidth=1)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c, cax=cax)
cb.ax.tick_params(labelsize=10, length=6, direction='out')
cb.set_label(r'$1+\delta$', fontsize=15)
cb.ax.tick_params(labelsize=15)
ax.set_xlabel("$x\,$(cMpc$h^{-1}$)", fontsize=14)
ax.set_ylabel("$y\,$(cMpc$h^{-1}$)", fontsize=14)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=15,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

# Create a ScalarFormatter for the colorbar
# formatter = ScalarFormatter(useMathText=True)
# formatter.set_powerlimits((0, -3))  # Adjust power limits for scientific notation

# # Access the colorbar's axis and apply the formatter
# cb.ax.yaxis.set_major_formatter(formatter)

if field_type=="dm":
    ax.set_title(fr"$z\,$ = {rounded_redshift}, (DM PROJ)", fontsize=15)
elif field_type=="gas":
    ax.set_title(fr"$z\,$ = {rounded_redshift}, (Gas PROJ)", fontsize=15)
plt.tight_layout()


if field_type=="dm":
    if add_haloes == 1:
        save_proj=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/DM_density/Projection/Haloes/"
    elif add_haloes == 0:
        save_proj=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/DM_density/Projection/No_Haloes/"
    else:
        print("Invalid value for Haloes")
elif field_type=="gas":
    if add_haloes == 1:
        save_proj=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/Gas_density/Projection/Haloes/"
    elif add_haloes == 0:
        save_proj=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/Gas_density/Projection/No_Haloes/"
    else:
        print("Invalid values for Haloes")
else:
    print("Invalid Field Type")

os.makedirs(save_proj,exist_ok=True)
plt.savefig(save_proj+f"{snapshot_number:03}.png",dpi=200)


# --- Plot Slice ---
fig, ax = plt.subplots(figsize=(8, 6))
c=ax.imshow(slice_data.T, cmap="viridis", aspect='equal', origin='lower',
              extent=[-box_length/2, box_length/2, -box_length/2, box_length/2], norm=colors.LogNorm(vmin=np.min(slice_data), vmax=np.max(slice_data)))

if add_haloes == 1:
    print("Adding Haloes on Slice Plot")
    # Convert position to grid units:
    # Shift origin back to 0
    posi_grid = (posi + box_length/2) * (N / box_length)  # now between [0, N]

    # Select haloes within slice thickness:
    # We want haloes with z between (N//2 - 1) and (N//2 + 1)
    z_pos = posi_grid[:, 2]
    halo_mask = (z_pos > (N//2 - 1)) & (z_pos < (N//2 + 1))

    if np.sum(halo_mask) == 0:
        print('No haloes found in the slice region. Skipping slice plot.')
        sys.exit()
    else:
        print(f"Found {np.sum(halo_mask)} haloes in slice region. Plotting them.")
        posi_slice = posi[halo_mask]
        # Now plot these haloes on top of slice
        ax.scatter(posi_slice[:, 0], posi_slice[:, 1], facecolors='none', edgecolors='white', 
               s=20*(mass[halo_mask]/m_min), linewidth=1)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(c, cax=cax)
cb.ax.tick_params(labelsize=10, length=6, direction='out')
cb.set_label(r'$1+\delta$', fontsize=15)
cb.ax.tick_params(labelsize=15)
ax.set_xlabel("$x\,$(cMpc$h^{-1}$)", fontsize=14)
ax.set_ylabel("$y\,$(cMpc$h^{-1}$)", fontsize=14)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=15,direction='in')
ax.tick_params(axis='both', which='minor', length=3, width=1, direction='in')
ax.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

# Create a ScalarFormatter for the colorbar
# formatter = ScalarFormatter(useMathText=True)
# formatter.set_powerlimits((0, -3))  # Adjust power limits for scientific notation

# # Access the colorbar's axis and apply the formatter
# cb.ax.yaxis.set_major_formatter(formatter)

if field_type=="dm":
    ax.set_title(fr"$z\,$ = {rounded_redshift}, (DM SLC)", fontsize=15)
elif field_type=="gas":
    ax.set_title(fr"$z\,$ = {rounded_redshift}, (Gas SLC)", fontsize=15)
plt.tight_layout()

if field_type=="dm":
    if add_haloes == 1:
        save_slc=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/DM_density/Slice/Haloes/"
    elif add_haloes == 0:
        save_slc=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/DM_density/Slice/No_Haloes/"
    else:
        print("Invalid value for Haloes")
elif field_type=="gas":
    if add_haloes == 1:
        save_slc=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/Gas_density/Slice/Haloes/"
    elif add_haloes == 0:
        save_slc=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{levelmin:03}_len{box_length}/Gas_density/Slice/No_Haloes/"
    else:
        print("Invalid values for Haloes")
else:
    print("Invalid Field Type")

os.makedirs(save_slc,exist_ok=True)
plt.savefig(save_slc+f"{snapshot_number:03}.png",dpi=200)

        
