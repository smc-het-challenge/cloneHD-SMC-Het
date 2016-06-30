# cloneHD - SMC-Het Challenge

This is a repository for the [cloneHD](http://www.sanger.ac.uk/science/tools/clonehd) submission to the [SMC-Het Challenge](http://dreamchallenges.org/project/home-upcoming/dream-9-5-icgc-tcga-dream-somatic-mutation-calling-tumor-heterogeneity-challenge-smc-het/). It solves the following sub-challenges:

* Sub-challenge 1A - Predicting normal contamination
* Sub-challenge 1B - Predicting number of subclones
* Sub-challenge 1C - Predicting subclone proportions
* Sub-challenge 2A - Determining mutation assignments to subclones
* Sub-challenge 2B - Co-clustering matrix

The method is described in a published [manuscript](http://www.cell.com/cell-reports/abstract/S2211-1247(14)00373-8) and on the [project wiki](https://www.synapse.org/#!Synapse:syn6148310/wiki/400535). The results of the challenge will become available on th [leaderboards](https://www.synapse.org/#!Synapse:syn2813581/wiki/303141) in the coming months.

## Local testing

To test the workflow locally, run the shell script:

    $ smchet_workflow.sh --vcf Tumour1.mutect.vcf --cna Tumour1.battenberg.txt --sample Tumour1

You can modify the default number of trials providing the optional flag `--trials <int>`, and the number of independent optimizations with the optional flag `--restarts <int>`. You can also set the seed with `--seed <int>`.

## Cloud testing

You first have to follow [these instructions](https://www.synapse.org/#!Synapse:syn2813581/wiki/303161) to sign up on Google Compute Engine, set up your working environment and set up your VM instance.

If you have the [Google Cloud SDK](https://cloud.google.com/sdk/) installed locally, you can SSH into the VM instance with gcloud:

    $ gcloud compute ssh ubuntu@planemo

To clone this repository, run the following command in a local directory. For example:

    $ cd /opt/galaxy/tools
    $ git clone https://github.com/ivazquez/cloneHD-SMC-Het.git

To set up the workflow in a container using `docker`:

    $ cd cloneHD-SMC-Het
    $ planemo docker_build .

This will build and compile cloneHD plus all dependencies.