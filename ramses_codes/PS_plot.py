import numpy as np
import matplotlib.pyplot as plt
import os

lev=8
length=100

z=np.loadtxt(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/redshift.txt")
n_outputs=len(z)

# Theoretical
# ================================
# Now, compute the theoretical linear power spectrum using CAMB

from camb import model, initialpower
import camb

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.143)
pars.InitPower.set_params(ns=0.965)
# For linear theory only:
pars.NonLinear = model.NonLinear_none

k_f=2*np.pi/length
k_n=(np.pi/length)*(2**lev)
for i in range(n_outputs):
    kd,powd,deltad=np.loadtxt(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_DM/psdm_{i+1:03}.txt",unpack=True)
    plt.loglog(kd,deltad,label="Simulation (DM)")
    kg,powg,deltag=np.loadtxt(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_gas/psgas_{i+1:03}.txt",unpack=True)
    plt.loglog(kg,deltag,label="Simulation (Gas)")
    # Use the simulation redshift:
    pars.set_matter_power(redshifts=[z[i]], kmax=k_n)
    results = camb.get_results(pars)
    k_camb, zs, P_k_linear = results.get_matter_power_spectrum(minkh=k_f, maxkh=k_n, npoints=2**(lev-1))
    Delta2_theory = (k_camb**3 / (2 * np.pi**2)) * P_k_linear[0]
    plt.loglog(k_camb, Delta2_theory, '-', label=f"Linear Theory (CAMB)")
    # ================================
    plt.title(f"z={z[i]:.2f}")
    plt.xlabel(r'$k$ [$h$/Mpc]')
    plt.ylabel(r'$\Delta^2(k)$')
    plt.legend()
    plt.axvline(x=k_n, color='g', linestyle='--', linewidth=1.5)  # Dashed red line
    plt.axvline(x=k_f, color='g', linestyle='--', linewidth=1.5)  # Dashed red line
    save_path=f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_plot/"
    os.makedirs(save_path,exist_ok=True)
    plt.savefig(f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_plot/psdm_{i+1:02}",dpi=150)
    plt.clf()  # Clear the current figure
