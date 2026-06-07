
# Metagenome profiler


Environmental monitoring is a critical component of pandemic readiness. The objective of the programming 
assignment is to develop a general method for characterizing the viral communities in a metagenomic sample. 
The input to your method is a short-read fastq file (examples below) containing sequences from a metagenomic 
sample. The output is a characterization of the virome that would be useful to someone working in pandemic 
preparedness (researcher, clinician, and/or public health official). Is up to you to decide the content and format 
of this report, both of which should be justified in your written and recorded reports.

Data was downloaded from these accession numbers using the SRA toolkit [here](https://github.com/ncbi/sra-tools). I used a script to parallelize the process here `src/download_sra.sh`.

Viral genome was downloaded from [here](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz)


# Environment setup
use `environment.yaml` to setup your python environment with conda

```bash
conda env create -n environment.yaml
conda activate basic
```
