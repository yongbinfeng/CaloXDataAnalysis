submissiontext = """#!/bin/bash
#SBATCH -J JOBNAME
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o LOGDIR/%x.%j.out
#SBATCH -e LOGDIR/%x.%j.err
#SBATCH -p nocona
"""

singularity_cmd = "singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_sbox/"
run_cmd = "python SCRIPTDIR/convertData.py -i ROOTFILE"

import os

current_dir = os.getcwd()
print("Current directory: ", current_dir)

script_dir = os.path.abspath(os.path.join(current_dir, ".."))
log_dir = os.path.abspath(os.path.join(current_dir, "log"))
jobname_prefix = "convertDataJob"


if not os.path.exists(log_dir):
    print(f"Creating log directory: {log_dir}")
    os.makedirs(log_dir)
    

def generate_submission_script(input_files):
    fnames = []
    for i, input_file in enumerate(input_files):
        jobname = f"{jobname_prefix}_{i}"
        run_cmd_tmp = run_cmd.replace("SCRIPTDIR", script_dir).replace("ROOTFILE", input_file)
        run_cmd_tmp = singularity_cmd + " bash -c \"" + run_cmd_tmp + "\""


        submissiontext_tmp = submissiontext.replace("JOBNAME", jobname).replace("LOGDIR", log_dir)

        # write the submission script
        fname = f"{log_dir}/submit_{jobname}.sh"
        with open(fname, "w") as f:
            f.write(submissiontext_tmp)
            f.write("\n\n")
            f.write(run_cmd_tmp)

        fnames.append(fname)

    submit_sh = f"{current_dir}/submit_all.sh"
    with open(submit_sh, "w") as f:
        f.write("#!/bin/bash\n")
        for fname in fnames:
            f.write(f"sbatch {fname}\n")

    os.system(f"chmod +x {submit_sh}")

    print(f"Submission script written to {submit_sh}")
    print("To submit jobs, run:")
    print(f"bash {submit_sh}")
    
    
if __name__ == "__main__":
    import json
    data = open("data_org_TB.json", "r")
    input_files = json.load(data)
    generate_submission_script(input_files)
