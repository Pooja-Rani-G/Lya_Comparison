import numpy as np

n_ph=8
no_of_photons=10**n_ph

snapnum=8
lev=10
lev_lo=8
length=200
snap_str=f"output_{snapnum:03}"
file_location=f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/xa.npy"
# Load xa field
xa = np.load(file_location)

# Compute delta_xa = xa / (1 + xa)
delta_xa = xa / (1.0 + xa)

# Save the result as raw binary file (float64), assuming your C++ code reads binary
save_location=f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/delta_xa.bin"
delta_xa.astype(np.float64).tofile(save_location)

# Alternatively, save as .npy if you plan to load it back in Python
# np.save("delta_xa.npy", delta_xa)
