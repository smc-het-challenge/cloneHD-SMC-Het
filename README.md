# SMC-Het Challenge

This is a repository for the [cloneHD](http://www.sanger.ac.uk/science/tools/clonehd) submission to the [SMC-Het Challenge](http://dreamchallenges.org/project/home-upcoming/dream-9-5-icgc-tcga-dream-somatic-mutation-calling-tumor-heterogeneity-challenge-smc-het/). You first have to follow [these instructions](https://www.synapse.org/#!Synapse:syn2813581/wiki/303161) to sign up on Google Compute Engine, set up your working environment and set up your VM instance.

## Testing

To test scoring 

    $ smchet_workflow.sh --vcf Tumour1.mutect.vcf --cna Tumour1.battenberg.txt --sample Tumour1

## Cloud

If you have the [Google Cloud SDK](https://cloud.google.com/sdk/) installed locally, you can SSH into the VM instance with gcloud:

    $ gcloud compute ssh ubuntu@planemo

To clone this repository, run the following command in a local directory:

    $ git clone https://github.com/ivazquez/smchet-challenge.git

To set up the workflow in a container using `docker`:

    $ cd smchet-challenge
    $ planemo docker_build .

This will build and compile cloneHD and cloneHD-tools, plus all dependencies.

* Sub-challenge 1A - Predicting normal contamination
* Sub-challenge 1B - Predicting number of subclones
* Sub-challenge 1C - Predicting subclone proportions
* Sub-challenge 2A - Determining mutation assignments to subclones
* Sub-challenge 2B - Co-clustering matrix