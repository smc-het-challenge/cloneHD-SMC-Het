cwlVersion: v1.0
class: CommandLineTool
label: cloneHD Driver
baseCommand: ["bash", "/opt/cloneHD_tool.sh"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/smchet-clonehd:0.1

inputs:    
  vcf:
    type: File
    inputBinding:
      prefix: --vcf
  cna:
    type: File
    inputBinding:
      prefix: --cna
  sample:
    type: string
    default: sample_name
    inputBinding:
      prefix: --sample
  outdir:
    type: Directory
    inputBinding:
      prefix: --output
  trials:
    type: int?
    default: 10
    inputBinding:
      prefix: --trials
  restarts:
    type: int?
    default: 10
    inputBinding:
      prefix: --restarts
  seed:
    type: int?
    default: 123
    inputBinding:
      prefix: --seed
  debug:
    type: boolean
    default: false
    inputBinding:
      prefix: --debug
  snv_pen_high:
    type: float
    default: 0.3 
    inputBinding:
      prefix: --snv-pen-high
  snv_pen_tree:
    type: float
    default: 0.001
    inputBinding:
      prefix: --snv-pen-tree
  llh_diff:
    type: int
    default: 20
    inputBinding:
      prefix: --llh-diff
  snv_min_per_cluster: 
    type: int
    default: 20
    inputBinding:
      prefix: --snv-min-per-cluster

outputs:
  cellularity:
    type: File
    outputBinding:
      glob: $(inputs.sample).1A.txt
  population:
    type: File
    outputBinding:
      glob: $(inputs.sample).1B.txt
  proportion:
    type: File
    outputBinding:
      glob: $(inputs.sample).1C.txt
  cluster_assignment:
    type: File
    outputBinding:
      glob: $(inputs.sample).2A.txt
  cocluster_assignment:
    type: File
    outputBinding:
      glob: $(inputs.sample).2B.txt.gz
