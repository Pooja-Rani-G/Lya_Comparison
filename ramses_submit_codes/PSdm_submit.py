import os
import subprocess

# Base command for running the script
base_command = ("/usr/bin/time -v /user1/poojarani/Lya_Comparison/ramses_codes/psdm {lev} {length} {output_folder} ")

lev = 8 # Refinement level
lengths = 100  # Box length(s)

# List of output numbers formatted correctly
output_numbers = [x for x in range(1,10)]

# Define the main job script directory
job_script_dir = "/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_DM/job_scripts/".format(
    lev=lev, length=lengths)
os.makedirs(job_script_dir, exist_ok=True)

# Define the log directory for .o and .e files
log_dir = "/user1/poojarani/Lya_Comparison/ramses_analysis/lev{lev:03}_len{length}/PS_DM/logs/".format(
    lev=lev, length=lengths)
os.makedirs(log_dir, exist_ok=True)  # Create log directory if it doesn't exist

# Iterate over each length and output number
for output_number in output_numbers:
    # Generate the command
    command = base_command.format(output_folder=output_number, length=lengths, lev=lev)

    # Create a unique job script name. For clarity, we use the same base for log files.
    job_script_name = f"PS_{output_number}_{lev}_{lengths}.pbs"
    job_script_path = os.path.join(job_script_dir, job_script_name)

    # Use the base name (without extension) for log files
    job_log_base = os.path.splitext(job_script_name)[0]
    output_log = os.path.join(log_dir, f"{job_log_base}.o")
    error_log = os.path.join(log_dir, f"{job_log_base}.e")

    # Delete old log files if they exist
    if os.path.exists(output_log):
        os.remove(output_log)
    if os.path.exists(error_log):
        os.remove(error_log)

    # Define log file names in the PBS script with full paths
    job_script_content = f"""#!/bin/bash -l
#PBS -N {job_log_base}
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -o {output_log}
#PBS -e {error_log}

cd $PBS_O_WORKDIR

{command}
"""

    # Write the job script to a file
    with open(job_script_path, "w") as script_file:
        script_file.write(job_script_content)

    # Submit the job script using qsub
    print(f"Submitting job script: {job_script_path}")
    result = subprocess.run(f"qsub {job_script_path}", shell=True, text=True, capture_output=True)

    # Check if the submission was successful
    if result.returncode != 0:
        print(f"Failed to submit job: {job_script_path}")
        print(f"Error: {result.stderr}")
        break
    else:
        print(f"Job submitted successfully: {job_script_name}")
