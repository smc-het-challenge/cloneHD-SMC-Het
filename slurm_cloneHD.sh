#$!/user/bin/bash
#slurm_cloneHD.sh
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -V

INDIR=$1
TUMOR=$2

sbatch -A spellmanlab --get-user-env ./run_cloneHD.sh -i $INDIR -t $TUMOR
