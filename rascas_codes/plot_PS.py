import numpy as np
import matplotlib.pyplot as plt
import os

n_ph=8
no_of_photons=10**n_ph

snapnum=8
lev=10
for lev_lo in [6,8]:
    for length in [100,200]:

        k_f=2*np.pi/length
        k_n=(np.pi/length)*(2**lev_lo)

        kd,powd,deltad=np.loadtxt(f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/output_{snapnum:03}/p1e{n_ph}/xaPS.txt",unpack=True)
        plt.loglog(kd,deltad,label=f"{length}-{2**lev} to {length}-{2**lev_lo}")

        plt.title(fr"PS for $\frac{{x_a}}{{1 + x_a}}$ with $10^{{{n_ph}}}$ photons", fontsize=12)
        plt.xlabel(r'$k$ [$h$/Mpc]')
        plt.ylabel(r'$\Delta^2(k)$')
        plt.legend()
        # plt.legend()
        # plt.axvline(x=k_n, color='g', linestyle='--', linewidth=1.5)  # Dashed red line
        # plt.axvline(x=k_f, color='g', linestyle='--', linewidth=1.5)  # Dashed red line
        # save_path=f"/user1/poojarani/Lya_Comparison/rascas_analysis/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/output_{snapnum:03}/p1e{n_ph}/"
        # os.makedirs(save_path,exist_ok=True)
        # plt.savefig(f"/user1/poojarani/Lya_Comparison/rascas_analysis/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/output_{snapnum:03}/p1e{n_ph}/PSxa.png",dpi=150)
        # plt.clf()  # Clear the current figure
os.makedirs("/user1/poojarani/Lya_Comparison/rascas_analysis/PS_xa/",exist_ok=True)
plt.savefig("/user1/poojarani/Lya_Comparison/rascas_analysis/PS_xa/All_in_one.png")
