#!/bin/bash
# run_cloneHD.sh
#SBATCH --partition=exacloud
#SBATCH --output=cloneHD-%j.out
#SBATCH --error=cloneHD-%j.err
#SBATCH --job-name=run_smchet-clonehd
#SBATCH --gres disk:50
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
        echo "Usage: $0 [ -i INDIR -t TUMOR ]"
        echo
        echo " [-t TUMOR]   - Full path to the directory *containing* the  'tumors' subdirectory; where 'tumors' "
        echo "                holds a subdirectory for each tumor ID;  each within which resides VCF and CNA "
        echo "                data specific to that tumor [e.g., /<full>/<path>/<to>/tumors/T0_noXY]."
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

ALPHA="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $ALPHA
OUTDIR=$ALPHA/cloneHD_outputs
#The following file should exist and must remain unchanged
TMPJSON=$ALPHA/run_cloneHD.json~
#The following file may not exist, but will be overwritten if it does
JSON=$ALPHA/run_cloneHD.json
CWL=$ALPHA/run_cloneHD.cwl
VCF=$TUMOR/*mutect.vcf
CNA=$TUMOR/*battenberg.txt
ID=`basename $TUMOR`
TIMING=$OUTDIR/timing

chmod -R 770 ./

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

sed -e "s|input_vcf|$WORKDIR\/`basename $VCF`|g" -e "s|input_cna|$WORKDIR\/`basename $CNA`|g" -e  "s|sample_id|$ID|g" -e "s|output_dir|$WORKDIR|g" $TMPJSON > $JSON
chmod 775 $JSON
cp $VCF $CNA $WORKDIR

{ time cwltool $CWL $JSON } 2> ${TIMING}/cloneHD_runtime.txt

if [ ! -z $WORKDIR/$ID.1A.txt ]; then mv $WORKDIR/$ID.1A.txt $WORKDIR/cellularity.predfile; fi
if [ ! -z $WORKDIR/$ID.1B.txt ]; then mv $WORKDIR/$ID.1B.txt $WORKDIR/population.predfile; fi
if [ ! -z $WORKDIR/$ID.1C.txt ]; then mv $WORKDIR/$ID.1C.txt $WORKDIR/proportion.predfile; fi
if [ ! -z $WORKDIR/$ID.2A.txt ]; then mv $WORKDIR/$ID.2A.txt $WORKDIR/cluster_assignment.predfile; fi
if [ ! -z $WORKDIR/$ID.2B.txt.gz ]; then mv $WORKDIR/$ID.2B.txt.gz $WORKDIR/cocluster_assignment.predfile; fi

tar -czf ../$OUTDIR/cloneHD.tar.gz $WORKDIR/*predfile
rm -rf $WORKDIR
