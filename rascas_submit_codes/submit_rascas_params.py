import numpy as np
import os
import subprocess
import math

n_ph=8
no_of_photons=10**n_ph

snapnum=8
lev=10
lev_lo=6
length=200
cores = 6
cpus=32*cores
snap_str=f"output_{snapnum:03}"
# Output folder structure
output_dir = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/"
os.makedirs(output_dir, exist_ok=True)

run_id=f"prascas_p1e{n_ph}"
# Template and generated file paths
cfg_template_path = "/user1/poojarani/Lya_Comparison/mittal_rascas/params_rascas.cfg"
cfg_output_path = os.path.join(output_dir, "params_rascas.cfg")
job_script_path = os.path.join(output_dir, f"{run_id}.pbs")
output_log = os.path.join(output_dir, f"{run_id}.o")
error_log = os.path.join(output_dir, f"{run_id}.e")

# Dynamic parameter injection
DomDumpDir = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/"
PhotonICFile = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/p1e{n_ph}.dat"
fileout = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/irrelevant.dat"
pa_fileout = f"/user1/poojarani/Lya_Comparison/rascas_data/{length}cMpc_{2**lev}/{length}cMpc_{2**lev_lo}/{snap_str}/p1e{n_ph}/Pa.dat"
nbundle = math.ceil(no_of_photons/(10*(cpus-1)))
# Step 1: create the config
with open(cfg_template_path, "r") as f:
    cfg_template = f.read()

cfg_filled = cfg_template.format(
    DomDumpDir=DomDumpDir,
    PhotonICFile=PhotonICFile,
    fileout=fileout,
    pa_fileout=pa_fileout,
    nbundle=nbundle
)

with open(cfg_output_path, "w") as f:
    f.write(cfg_filled)

# Step 2: create the job script
job_script = f"""#!/bin/bash -l
#PBS -N {run_id}
#PBS -l nodes={cores}:ppn=32
#PBS -l walltime=72:00:00
#PBS -q temp
#PBS -l cput=1000:00:00
#PBS -o {output_log}
#PBS -e {error_log}

cd $PBS_O_WORKDIR
module load mpi/openmpi-x86_64

/usr/bin/time -v /usr/lib64/openmpi/bin/mpirun -np {cpus} -machinefile $PBS_NODEFILE  /user1/poojarani/Lya_Comparison/mittal_rascas/rascas  {cfg_output_path}
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