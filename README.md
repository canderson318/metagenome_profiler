
# Metagenome profiler
<!-- # Assignment -->
<!-- Environmental monitoring is a critical component of pandemic readiness. The objective of the programming assignment is to develop a general method for characterizing the viral communities in  metagenomic sample. The input to your method is a short-read fastq file (examples below) containing sequences from a metagenomic sample. The output is a characterization of the virome that would be useful to someone working in pandemic preparedness (researcher, clinician, and/or public health official). Is up to you to decide the content and format of this report, both of which should be justified in your written and recorded reports. -->

<!-- 
# Data
Data was downloaded from these accession numbers using the SRA toolkit [here](https://github.com/ncbi/sra-tools). I used a script to parallelize the process here `src/download_sra.sh`.
- [https://www.ncbi.nlm.nih.gov/sra/SRX8958936%5Baccn%5D](https://www.ncbi.nlm.nih.gov/sra/SRX8958936%5Baccn%5D)
- (https://www.ncbi.nlm.nih.gov/sra/SRX8927795%5Baccn%5D)[https://www.ncbi.nlm.nih.gov/sra/SRX8927795%5Baccn%5D]

The viral genome reference database was download was downloaded from [here](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz) and queried using Kraken2 ((installation guide)[https://github.com/DerrickWood/kraken2/wiki/Manual#installation]). -->
This is a 


# Input Files
This package expects sample files in the `.fastq` format and reference genomes in the `.fna` format. 

# Code
## Environment setup
use `environment.yaml` to setup your python environment with conda

```bash
conda env create -n environment.yaml
conda activate metagenome_profiler
```

# Directory Format
This package depends on the following directory structure:
```
.
  |-- environment.yaml
  |-- in
  |   |-- host_genome
  |   |   |-- GCF...
  |   |   |   `-- GCF....fna
  |   |-- kraken_2_viral_db
  |   `-- sample.fastq
  |-- notes.typ
  |-- out
  |   `-- <<outputs>>
  `-- src 
```

# Usage
## QC

