#!/bin/bash
# run_cloneHD.sh
#SBATCH --partition=exacloud
#SBATCH --output=cloneHD-%j.out
#SBATCH --error=cloneHD-%j.err
#SBATCH --job-name=run_smchet-clonehd
#SBATCH --gres disk:1024
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:45:00

function usage()
{
        echo "run_cloneHD.sh"    " [a.k.a. *this* script] "
        echo "Author: Kami E. Chiotti "
        echo "Date: 10.06.19"
        echo
        echo "A wrapper for the SMC-Het Dream Challenge submission, 'ivazquez/smchet-challenge'. It "
        echo "assembles and executes the command to run cloneHD_tool.sh, then collects and returns "
        echo "the output files as a tarball. "
        echo
        echo "NOTE #1: The 'run_cloneHD.cwl' CWLtool and the 'run_cloneHD.json~' template file must "
        echo "         be in the same directory as *this* script. "
        echo "NOTE #2: The output will be written to a subdirectory called 'outputs' in the same directory"
        echo "         as *this* script. "
        echo "NOTE #3: As assembled now, cloneHD_tool.sh must be run via *this* script."
        echo
        echo "Usage: $0 [ -t TUMOR ]"
        echo
        echo " [-t TUMOR]    - Full path to the directory *containing* the  'tumors' subdirectory; where 'tumors' "
        echo "                 holds a subdirectory for each tumor ID;  each within which resides VCF and CNA "
        echo "                 data specific to that tumor [e.g., /<full>/<path>/<to>/tumors/T0_noXY]."
        exit
}

TUMOR=""
while getopts ":t:h" Option
        do
        case $Option in
                t ) TUMOR="$OPTARG" ;;
                h ) usage ;;
                * ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
        esac
done
shift $(($OPTIND - 1))

if [[ "$TUMOR" == "" ]]
then
        usage
fi

source /home/groups/EllrottLab/activate_conda

ALPHA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $ALPHA
TMPJSON=$ALPHA/run_cloneHD.json~
CWL=$ALPHA/run_cloneHD.cwl
VCF=$TUMOR/*mutect.vcf
CNA=$TUMOR/*battenberg.txt
OUTDIR=$ALPHA/cloneHD_outputs

if [ ! -e $TUMOR ];
then
    echo "Input data directory for $TUMOR does not exist"; exit -1 ;
fi
if [ ! -e $OUTDIR ];
then
    mkdir -p $OUTDIR
fi

WORKDIR=`mktemp -d -p /mnt/scratch/ cloneHD.XXX`
chmod 775 $WORKDIR
chmod g+s $WORKDIR
JSON=$WORKDIR/run_cloneHD.json

sed -e "s|input_vcf|$WORKDIR\/`basename $VCF`|g" -e "s|input_cna|$WORKDIR\/`basename $CNA`|g" -e  "s|sample_name|`basename $TUMOR`|g" -e "s|output_dir|$WORKDIR|g" $TMPJSON > $JSON
cp $VCF $CNA $CWL $WORKDIR

cd $WORKDIR
cwltool `basename $CWL` `basename $JSON`

if [ ! -z $WORKDIR/`basename $TUMOR`.1A.txt ]; then mv $WORKDIR/`basename $TUMOR`.1A.txt $WORKDIR/cellularity.predfile; fi
if [ ! -z $WORKDIR/`basename $TUMOR`.1B.txt ]; then mv $WORKDIR/`basename $TUMOR`.1B.txt $WORKDIR/population.predfile; fi
if [ ! -z $WORKDIR/`basename $TUMOR`.1C.txt ]; then mv $WORKDIR/`basename $TUMOR`.1C.txt $WORKDIR/proportion.predfile; fi
if [ ! -z $WORKDIR/`basename $TUMOR`.2A.txt ]; then mv $WORKDIR/`basename $TUMOR`.2A.txt $WORKDIR/cluster_assignment.predfile; fi
if [ ! -z $WORKDIR/`basename $TUMOR`.2B.txt.gz ]; then mv $WORKDIR/`basename $TUMOR`.2B.txt.gz $WORKDIR/cocluster_assignment.predfile; fi

tar -czf cloneHD.tar.gz *predfile
rsync -a cloneHD.tar.gz $OUTDIR

cd $ALPHA
rm -rf $WORKDIR

