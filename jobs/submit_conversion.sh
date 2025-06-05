#!/bin/bash
#SBATCH -J dataconvertion
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o /home/yofeng/CaloX/DataAnalysis/jobs/log/%x.%j.out
#SBATCH -e /home/yofeng/CaloX/DataAnalysis/jobs/log/%x.%j.err
#SBATCH -p nocona


singularity run --cleanenv --bind /lustre:/lustre /lustre/work/yofeng/SimulationEnv/alma9forgeant4_v1_sbox/ bash -c "cd /home/yofeng/CaloX/DataAnalysis && python convertData.py"
