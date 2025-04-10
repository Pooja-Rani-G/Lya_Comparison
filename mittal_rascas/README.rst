My version of public code RASCAS used for the paper `arXiv:2311.03447 <https://arxiv.org/abs/2311.03447>`_. Compared to original version introduced by `Michel-Dansac et al. (2020) <https://www.aanda.org/articles/aa/full_html/2020/03/aa34961-18/aa34961-18.html>`_ my version has Hubble flow implemented, stopping criterion is decided by medium transparency and finally, computes the scattering rate, :math:`P_{\alpha}=P_{\alpha}(x,y,z)` throughout the box.

The main branch to use is 'master'. It has Hubble flow implemented, stopping criterion is decided by medium transparency and computes the scattering rate, :math:`P_{\alpha}=P_{\alpha}(x,y,z)` throughout the box.

The other branches have one thing or the other implemented differently. I created these during the developement process. You may ignore these branches. Nevertheless, following is a description of each of these.

1. 'Pa-by-number': it computes :math:`P_{\alpha}` using point-based technique. Everything else is same as in the master branch. To read more on this see Sec. 2.3.3 and App. B in the paper.

2. 'all-but-bulk', 'all-but-nhi', and 'all-but-temp' have everything same as in the master branch except they adopt a uniform bulk velocity, uniform hydrogen number density, or uniform temperature, respectively, throught the cosmological box. The uniform value needs to be supplied and entered manually in modules. (Contact me for the details.)

3. 'no-source-vel-or-mass': the master version distributes the intial photons according to the source mass and Doppler shifts the emitting photon according to the source velocity. However, in this branch the initial photon packet distribution is completely independent of source properties.

4. 'non-hubble': as the name suggests, this version does not implement the Hubble flow. Everything else is the same as in the master branch.

5. 'trilinear': in this version the bulk velocity field comes from an external file. In other versions, RASCAS reads everything from the same cosmological box (which is generated from RAMSES). But here, only the velocity field comes from a possibly higher-resolution cosmological box. Before using this you need to run the python code 'vel_gas.py'. The main scientific difference is as follows. In the main branch, throughout the cell the photon sees the same gas velocity. However, in this branch the velocity is dependent on position - approximated using a trilinear interpolation according the velocity of neighbouring cells.

6. 'zero-kelvin': in this version the thermal motion of gas particles is not accounted for. Consequently, a Lorentzian line profile is used. Note that temperature plays no role in this version. Also, note that the effect of source verlocity and mass on the emitted photons is not accounted for.

In case of confusion, please contact me.
