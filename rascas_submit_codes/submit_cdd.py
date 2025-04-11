import numpy as np
import os
import subprocess

#################### Finding Appropriate lesser resolution box ##################
# Step 1: Original high-resolution parameters
lev_hi = 10
length = 200  # in cMpc
snapnum_hi = 8  # <- this is probably 1-based, so subtract 1 for index

# Path to redshift file of high-resolution level
redshift_hi_path = f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev_hi:03}_len{length}/redshift.txt"
z_hi = np.loadtxt(redshift_hi_path)
# print("z_hi",z_hi)
z_target = z_hi[snapnum_hi - 1]  # Corrected indexing here

# Step 2: Lower-resolution level to find best match
lev_lo = 6
redshift_lo_path = f"/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev_lo:03}_len{length}/redshift.txt"
z_lo = np.loadtxt(redshift_lo_path)

# Step 3: Find snapnum in lower level that has closest redshift
snapnum_lo = int(np.argmin(np.abs(z_lo - z_target)))
snapnum_lo=snapnum_lo+1
print("snapnum_lo",snapnum_lo)
z_matched = z_lo[snapnum_lo-1]

print(f"📌 Matching redshift z = {z_target:.5f} from lev {lev_hi} snap {snapnum_hi}")
print(f"👉 Closest match in lev {lev_lo}: snap {snapnum_lo} with z = {z_matched:.5f}")


# Derived parameters
run_id = f"snap{snapnum_hi}_{length}cMpc_{2**lev_hi}"
snap_str = f"output_{snapnum_hi:03}"

# Output folder structure
output_dir = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev_hi}/{length}cMpc_{2**lev_lo}/{snap_str}/"
os.makedirs(output_dir, exist_ok=True)

########### Creating info.txt #########################
info_path = os.path.join(output_dir, "info.txt")

with open(info_path, "w") as f:
    f.write("# High resolution: lev_hi, snapnum_hi, z_hi\n")
    f.write(f"{lev_hi}\t{snapnum_hi}\t{z_target:.6f}\n\n")

    f.write("# Matched low resolution: lev_lo, snapnum_lo, z_lo\n")
    f.write(f"{lev_lo}\t{snapnum_lo}\t{z_matched:.6f}\n\n")

    f.write("# Final run parameters: lev, length [cMpc], snapnum\n")
    f.write(f"{lev_lo}\t{length}\t{snapnum_lo}\n")

# Template and generated file paths
cfg_template_path = "/user1/poojarani/Lya_Comparison/mittal_rascas/CDD.cfg"
cfg_output_path = os.path.join(output_dir, "CDD.cfg")
job_script_path = os.path.join(output_dir, f"{run_id}.pbs")
output_log = os.path.join(output_dir, f"{run_id}.o")
error_log = os.path.join(output_dir, f"{run_id}.e")

# Dynamic parameter injection
DomDumpDir = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev_hi}/{length}cMpc_{2**lev_lo}/{snap_str}/"
repository = f"/user1/poojarani/Lya_Comparison/ramses_data/lev{lev_lo:03}_len{length}/"

# Step 1: create the config
with open(cfg_template_path, "r") as f:
    cfg_template = f.read()

cfg_filled = cfg_template.format(
    DomDumpDir=DomDumpDir,
    repository=repository,
    snapnum=snapnum_lo
)

with open(cfg_output_path, "w") as f:
    f.write(cfg_filled)

# Step 2: create the job script
job_script = f"""#!/bin/bash -l
#PBS -N {run_id}
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -o {output_log}
#PBS -e {error_log}

cd $PBS_O_WORKDIR
module load mpi/openmpi-x86_64

/usr/bin/time -v /user1/poojarani/Lya_Comparison/mittal_rascas/CreateDomDump {cfg_output_path}
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
