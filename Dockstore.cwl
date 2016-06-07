#!/usr/bin/env cwl-runner

class: CommandLineTool
description: "A Docker container for the cloneHD command. See the [cloneHD](http://www.sanger.ac.uk/science/tools/clonehd) website for more information."
id: "cloneHD"
label: "cloneHD tool"

description: |
    The cloneHD subclonal reconstruction workflow for the ICGC PanCancer Analysis of Whole Genomes (PCAWG) 
		project. For more information see the PCAWG project [page](https://dcc.icgc.org/pcawg) and our 
		GitHub [page](https://github.com/ICGC-TCGA-PanCancer) for our code including the source for 
		[this workflow](https://github.com/ICGC-TCGA-PanCancer/CGP-Somatic-Docker).
    ```
    Usage:
    # fetch CWL
    $> dockstore cwl --entry quay.io/ivazquez/pcawg-clonehd-workflow:2.0.0 > Dockstore.cwl
    # make a runtime JSON template and edit it
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/ivazquez/pcawg-clonehd-workflow:2.0.0 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "http://orcid.org/0000-0003-0427-2639"
  foaf:name: Ignacio Vazquez-Garcia
  foaf:mbox: "mailto:ivg@sanger.ac.uk"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/ivazquez/pcawg-clonehd-workflow:1.25-2"
  - { import: node-engine.cwl }

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000
    description: "the process requires at least 4G of RAM"

inputs:
  - id: "#mem_gb"
    type: int
    default: 4
    description: "The memory, in GB, for the reporting tool"
    inputBinding:
      position: 1

  - id: "#bam_input"
    type: File
    description: "The BAM file used as input, it must be sorted."
    inputBinding:
      position: 2

outputs:
  - id: "#bamstats_report"
    type: File
    outputBinding:
      glob: bamstats_report.zip
    description: "A zip file that contains the HTML report and various graphics."

baseCommand: ["bash", "/usr/local/bin/cloneHD"]