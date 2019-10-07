#$!/user/bin/bash
#slurm_cloneHD.sh
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V

sbatch -A spellmanlab --get-user-env ./run_cloneHD.sh $INDIR $TUMOR
