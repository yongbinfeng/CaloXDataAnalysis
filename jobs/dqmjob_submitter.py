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

singularity_cmd = "singularity run --cleanenv --bind /lustre:/lustre -B /usr/lib64/libgif.so.7:/usr/lib64/libgif.so.7 -B /usr/lib64/libtiff.so.5:/usr/lib64/libtiff.so.5 -B /usr/lib64/libjpeg.so.62:/usr/lib64/libjpeg.so.62 -B /usr/lib64/libjbig.so.2.1:/usr/lib64/libjbig.so.2.1 /lustre/work/yofeng/SimulationEnv/alma9forgeant4_sbox/"
run_cmd = "cd SCRIPTDIR && python makeDQMHists.py --run RUNNUMBER && python makeDQMPlots.py --run RUNNUMBER && python checkServiceDRS.py --run RUNNUMBER && python makeFERSEnergyPlots.py --run RUNNUMBER && python checkMCP.py --run RUNNUMBER && python makeDRSPeakTS.py --run RUNNUMBER"

import os

current_dir = os.getcwd()
print("Current directory: ", current_dir)

script_dir = os.path.abspath(os.path.join(current_dir, ".."))
log_dir = os.path.abspath(os.path.join(current_dir, "log"))
jobname_prefix = "makeDQMPlots"


if not os.path.exists(log_dir):
    print(f"Creating log directory: {log_dir}")
    os.makedirs(log_dir)
    

def generate_submission_script(runlist):
    fnames = []
    for i, run_number in enumerate(runlist):
        jobname = f"{jobname_prefix}_{i}"
        run_cmd_tmp = run_cmd.replace("SCRIPTDIR", script_dir).replace("RUNNUMBER", str(run_number))
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
    runList = range(1350, 1448)
    generate_submission_script(runlist=runList)
