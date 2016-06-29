#!/bin/bash

# set an initial value for the flag
input_vcf=false
input_cna=false
output_dir=false
sample_name=false
show_help=false
debug=false

if [ $# -eq 0 ];
then
    show_help=true
fi

# read the options
TEMP=`getopt -o v:c:s:o:hd --long vcf:,cna:,sample:,output:,help,debug -n 'smchet_workflow.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables
while true ; do
    case "$1" in
        -v|--vcf) input_vcf=$2 ; shift 2 ;;
        -c|--cna) input_cna=$2 ; shift 2 ;;
        -s|--sample) sample_name=$2 ; shift 2 ;;
        -o|--output) output_dir=$2 ; shift 2 ;;
        -h|--help) show_help=true ; shift ;;
        -d|--debug) debug=true ; shift ;;
        --) shift ; break ;;
        *) echo $2 "Internal error!" ; exit 1 ;;
    esac
done

if [[    $input_vcf == false
      || $input_cna == false
      || $sample_name == false
      || $output_dir == false
    ]];
then
    show_help=true
fi

if [ $show_help == true ];
then
    echo "usage: $0 [options]"
    echo "options:"
    echo "    -v, --vcf arg       VCF file describing the SNV data"
    echo "    -c, --cna arg       Battenberg-format segmentation file describing the CNA data"
    echo "    -s, --sample arg    sample name"
    echo "    -o, --output arg    write output into this directory named dir/sample.cloneHD.gz"
    echo "    -d, --debug         turns on debugging"
    echo "    -h, --help          this text"
    echo "    add --debug for debugging output"
    echo "Run cloneHD and save results in the output directory."
    exit
fi

prefix=$sample_name
snv=$prefix.snv.txt
mean_tcn=$prefix.mean_tcn.txt
avail_cn=$prefix.avail_cn.txt

### SNV parser ###
python /opt/snv_parser.py \
	--variant-type 'mutect-smchet' \
	--output-snvs ${snv} \
	${input_vcf}

### CNA parser ###
python /opt/cna_parser.py \
	--cna-format 'battenberg-smchet' \
	--cellularity 1.0 \
	--mean-tcn ${mean_tcn} \
	--avail-cn ${avail_cn} \
	${input_cna}

### cloneHD ###
declare -A clones_to_clusters=( [1]=0 [2]=1 [3]=3 )
n_clones=1

while [ $n_clones -le 3 ]
do
	
	summary[$n_clones]=$prefix.Nc$n_clones.summary.txt
	snv_posterior[$n_clones]=$prefix.Nc$nclones.snv.posterior.txt
	
	n_clusters=${clones_to_clusters[$n_clones]}
	
	snv_fprate=`wc -l $snv | awk '{print 800/$1}'`

	/opt/cloneHD/build/cloneHD \
		--pre $prefix.Nc$n_clones \
		--snv $snv \
		--seed 123 \
		--trials 10 \
		--force $n_clones \
		--max-tcn 8 \
		--restarts 10 \
		--mean-tcn $mean_tcn \
		--avail-cn $avail_cn \
		--snv-rnd 1E-2 \
		--snv-fpfreq 5E-2 \
		--snv-fprate $snv_fprate \
		--learn-cluster-w $n_clusters \
		--snv-pen-high 3E-1 \
		--print-all 0
	
	n_clones=`expr $n_clones + 1`
	
done

# ### Model selection and SMC-Het conversion ###
# python ${dir}/cloneHD-tools/clonehd/report.py \
# 	--summary ${summary[*]} \
# 	--snv-posterior ${snv_posterior[*]} \
# 	--output $prefix

## Model selection ###
perl /opt/subclone_model_selection.pl \
	-i $summary[1] -j $snv_posterior[1] \
	-k $summary[2] -l $snv_posterior[2] \
	-m $summary[3] -n $snv_posterior[3] \
	-a 10.0 -s 50.0 \
	-o $prefix

### SMC-Het conversion ###
assignment=$prefix.mutation_assignment.txt
perl /opt/smchet_conversion.pl -i $assignment -o $prefix
#/opt/run_metrics $assignment | gzip > $prefix.2B.txt.gz
