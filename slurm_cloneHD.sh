#$!/user/bin/bash
#slurm_cloneHD.sh
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V

TUMOR=$1

sbatch -A spellmanlab --get-user-env ./run_cloneHD.sh -t $TUMOR
