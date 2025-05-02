import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

lev_hi=10
lev_lo=8
length=100
snap=8
exp_of_photons=8
path= f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev_hi}/{length}cMpc_{2**lev_lo}/output_{snap:03}/p1e{exp_of_photons}/"

############### 21cm Signal ####################
# Load the data
T21 = np.load(path+"T21.npy")
print(np.shape(T21))
# Compute min and max values
T21_avg = np.mean(T21, axis=2)  # Shape will be (Nx, Ny)
vmin, vmax = np.min(T21), np.max(T21)
print("In 3D")
print("Min value:", vmin,"mK")
print("Max value:", vmax, "mK")
print("Average value of T21 :", np.mean(T21))
print("After averaging in z-direction")
vmin, vmax = np.min(T21_avg), np.max(T21_avg)
print("Min value:", vmin, "mK")
print("Max value:", vmax, "mK")
T21_avg=T21_avg/1000
# Plot the x-y plane at the chosen z-index
fig, ax = plt.subplots(figsize=(6, 6))
# Create custom colormap
custom_cmap = LinearSegmentedColormap.from_list("custom", ["blue", "red", "green"])

im = ax.imshow(T21_avg[:, :], cmap=custom_cmap, origin='lower',
               extent=[-length/2, length/2, -length/2, length/2])
# Divider for colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)  # control size and spacing

# Colorbar
cbar = plt.colorbar(im, cax=cax)
cbar.set_label(r'$T_{21}(K)$', fontsize=12)

# Labels
ax.set_xlabel(r'$x\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)
ax.set_ylabel(r'$y\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)

plt.tight_layout()
# Optional: adjust tick font size
ax.tick_params(labelsize=10)
cbar.ax.tick_params(labelsize=10)

save_path=f"/user1/poojarani/Lya_Comparison/rascas_analysis/{length}cMpc_{2**lev_hi}/{length}cMpc_{2**lev_lo}/output_{snap:03}/p1e{exp_of_photons}/"
import os
os.makedirs(save_path,exist_ok=True)
# Save the figure
plt.savefig(save_path+"T21.png")


plt.clf()  # Clears the current figure

############## lyman-a coupling #########################
from matplotlib.colors import LogNorm

# Load the data
xa = np.load(path+"xa.npy")
xa_avg = np.mean(xa, axis=2)  # Shape will be (Nx, Ny)
# Compute min and max values
vmin, vmax = np.min(xa), np.max(xa)
print("In 3D")
print("Min value:", vmin)
print("Max value:", vmax)
print("Average value of xa :", np.mean(xa))

print("After averaging in z-direction")
vmin, vmax = np.min(xa_avg), np.max(xa_avg)
print("Min value:", vmin)
print("Max value:", vmax)


# Plot the x-y plane at the chosen z-index
fig, ax = plt.subplots(figsize=(6, 6))
im = ax.imshow(xa_avg[:, :], cmap='cool', origin='lower',
               extent=[-length/2, length/2, -length/2, length/2],norm=LogNorm(vmin=np.min(xa_avg), vmax=np.max(xa_avg)))
# Divider for colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)  # control size and spacing

# Colorbar
cbar = plt.colorbar(im, cax=cax)
cbar.set_label(r'$x_a$', fontsize=12)

# Labels
ax.set_xlabel(r'$x\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)
ax.set_ylabel(r'$y\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)

plt.tight_layout()
# Optional: adjust tick font size
ax.tick_params(labelsize=10)
cbar.ax.tick_params(labelsize=10)
# Save the figure
plt.savefig(save_path+"xa.png")

plt.clf()  # Clears the current figure

########## For Semelin's Comparison #############
xas=xa/(1+xa)
# Compute min and max values
vmin, vmax = np.min(xas), np.max(xas)
print("In 3D")
print("Min value:", vmin)
print("Max value:", vmax)
print("Average value of xa_semelin :", np.mean(xas))

fig, ax = plt.subplots(figsize=(6, 6))

# Plotting
im = ax.imshow(xas[:, :,2**(lev_lo-1)], cmap='turbo', origin='lower',
               extent=[-length/2, length/2, -length/2, length/2])

# Divider for colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)  # control size and spacing

# Colorbar
cbar = plt.colorbar(im, cax=cax)
cbar.set_label(r'${x_a}/{(1 + x_a)}$', fontsize=12)

# Labels
ax.set_xlabel(r'$x\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)
ax.set_ylabel(r'$y\ \mathrm{(cMpc\ h^{-1})}$', fontsize=12)

plt.tight_layout()
# Optional: adjust tick font size
ax.tick_params(labelsize=10)
cbar.ax.tick_params(labelsize=10)

# Save the figure
plt.savefig(save_path+"semelin_xa.png")


