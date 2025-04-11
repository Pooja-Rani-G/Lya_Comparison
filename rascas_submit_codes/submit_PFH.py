import numpy as np
import os
import subprocess

n_ph=8
no_of_photons=10**n_ph

snapnum=8
lev=10
length=200
lev_lo=6
snap_str=f"output_{snapnum:03}"
# Output folder structure
output_dir = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/"
os.makedirs(output_dir, exist_ok=True)

run_id=f"{snapnum}_p1e{n_ph}"
# Template and generated file paths
cfg_template_path = "/user1/poojarani/Lya_Comparison/mittal_rascas/PFH.cfg"
cfg_output_path = os.path.join(output_dir, "PFH.cfg")
job_script_path = os.path.join(output_dir, f"{run_id}.pbs")
output_log = os.path.join(output_dir, f"{run_id}.o")
error_log = os.path.join(output_dir, f"{run_id}.e")

# Dynamic parameter injection
haloesinfo = f"/user1/poojarani/Lya_Comparison/ramses_data/lev{lev:03}_len{length}/Halo_data/haloes_{snapnum:05}.dat"
outputfile = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/p1e{n_ph}.dat"

# Step 1: create the config
with open(cfg_template_path, "r") as f:
    cfg_template = f.read()

cfg_filled = cfg_template.format(
    haloesinfo=haloesinfo,
    outputfile=outputfile,
    nPhotonPackets=no_of_photons
)

with open(cfg_output_path, "w") as f:
    f.write(cfg_filled)

# Step 2: create the job script
job_script = f"""#!/bin/bash -l
#PBS -N {run_id}
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -o {output_log}
#PBS -e {error_log}

cd $PBS_O_WORKDIR
module load mpi/openmpi-x86_64

/usr/bin/time -v /user1/poojarani/Lya_Comparison/mittal_rascas/PhotonsFromHaloes {cfg_output_path}
"""

with open(job_script_path, "w") as f:
    f.write(job_script)

# Step 3: cleanup and submission
for log in [output_log, error_log]:
    if os.path.exists(log):
        os.remove(log)

print(f"Submitting job: {job_script_path}")
res = subprocess.run(f"qsub {job_script_path}", shell=True, text=True, capture_output=True)

if res.returncode != 0:
    print(f"❌ Job submission failed:\n{res.stderr}")
else:
    print(f"✅ Job submitted successfully: {res.stdout.strip()}")
