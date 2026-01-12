submissiontext = """#!/bin/bash
#SBATCH -J JOBNAME
#SBATCH -N 8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -o LOGDIR/%x.%j.out
#SBATCH -e LOGDIR/%x.%j.err
#SBATCH -p nocona
"""

env_cmd = """echo "This job is running on the following nodes:"
echo $SLURM_NODELIST
"""
JOBTMP="/tmp/root_cache_RUNNUMBER"

singularity_cmd = "singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_v3.sif"
run_cmd = "mkdir -p JOBTMP; export PCM_CACHE_DIR=JOBTMP; export TMPDIR=JOBTMP; cd SCRIPTDIR && python scripts/make_dqm_hists.py --run RUNNUMBER && python scripts/make_dqm_plots.py --run RUNNUMBER && python scripts/check_service_drs.py --run RUNNUMBER && python scripts/make_fers_energy_plots.py --run RUNNUMBER && python scripts/check_beam_composition.py --run RUNNUMBER"

import os

current_dir = os.getcwd()
print("Current directory: ", current_dir)

script_dir = os.path.abspath(os.path.join(current_dir, ".."))
log_dir = os.path.abspath(os.path.join(current_dir, "log"))
jobname_prefix = "calox"


if not os.path.exists(log_dir):
    print(f"Creating log directory: {log_dir}")
    os.makedirs(log_dir)
    

def generate_submission_script(runlist):
    fnames = []
    for i, run_number in enumerate(runlist):
        jobname = f"{jobname_prefix}_{run_number}"
        run_cmd_tmp = run_cmd.replace("SCRIPTDIR", script_dir).replace("JOBTMP",JOBTMP).replace("RUNNUMBER", str(run_number))
        run_cmd_tmp = singularity_cmd + " bash -c \"" + run_cmd_tmp + "\""


        submissiontext_tmp = submissiontext.replace("JOBNAME", jobname).replace("LOGDIR", log_dir)

        # write the submission script
        fname = f"{log_dir}/submit_{jobname}.sh"
        with open(fname, "w") as f:
            f.write(submissiontext_tmp)
            f.write("\n\n")
            f.write(env_cmd)
            f.write(run_cmd_tmp)

        fnames.append([run_number, fname])

    submit_sh = f"{current_dir}/submit_all.sh"
    with open(submit_sh, "w") as f:
        f.write("#!/bin/bash\n")
        for run_number, fname in fnames:
            f.write(f"sbatch --job-name=calox{run_number} {fname}\n")

    os.system(f"chmod +x {submit_sh}")

    print(f"Submission script written to {submit_sh}")
    print("To submit jobs, run:")
    print(f"bash {submit_sh}")
    
    
if __name__ == "__main__":
    runList = range(1350, 1528)
    generate_submission_script(runlist=runList)
