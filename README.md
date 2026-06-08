![tests](https://github.com/canderson318/metagenome_profiler/actions/workflows/tests.yml/badge.svg)

# Metagenome profiler
<!-- # Assignment -->
<!-- Environmental monitoring is a critical component of pandemic readiness. The objective of the programming assignment is to develop a general method for characterizing the viral communities in  metagenomic sample. The input to your method is a short-read fastq file (examples below) containing sequences from a metagenomic sample. The output is a characterization of the virome that would be useful to someone working in pandemic preparedness (researcher, clinician, and/or public health official). Is up to you to decide the content and format of this report, both of which should be justified in your written and recorded reports. -->

The objective of this package is to provide a general method for characterizing the viral communities in  metagenomic sample. The inputs include a short-read fastq file containing sequences from a metagenomic sample, a reference fna file, and a viral reference database (examples below). The output is a characterization of the virome that may be useful to someone working in monitoring viral communities in the environment. The tools include filtering methods like kmer-based host-alignment and phred-score filtering, and read classification using Kraken2's viral genome reference database. 

__Core components:__
- Host alignment based read filtering
- Phred score filtering
- Virome characterization (read classification)
- Abundance analysis
- Read edit distance analysis


# Input Files
This package expects sample files in the `.fastq` format and reference genomes in the `.fna` format.

The data used in the example below was downloaded from these accession links using the SRA toolkit ([sra-tools guide](https://github.com/ncbi/sra-tools)) ([src/download_sra.sh](src/download_sra.sh)).

- [https://www.ncbi.nlm.nih.gov/sra/SRX8958936%5Baccn%5D](https://www.ncbi.nlm.nih.gov/sra/SRX8958936%5Baccn%5D)
- [https://www.ncbi.nlm.nih.gov/sra/SRX8927795%5Baccn%5D](https://www.ncbi.nlm.nih.gov/sra/SRX8927795%5Baccn%5D)

The viral genome reference database was download was downloaded from [this link](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz) ([src/download_viral_db.sh](src/download_viral_db.sh)) and queried using Kraken2 ([installation guide](https://github.com/DerrickWood/kraken2/wiki/Manual#installation)). 

The reference genome, in this case a bat, was downloaded using the ncbi `datasets` tool ([src/download_bat_genome.sh](src/download_bat_genome.sh)).


# Directory - Format
This package depends on the following directory structure:
```
.
  |-- environment.yaml
  |-- in
  |   |-- host_genome # ncbi download "GCF...directory/GCF...file.fna "
  |   |   |-- GCF...
  |   |   |   `-- GCF....fna
  |   |-- kraken_2_viral_db 
  |   `-- sample.fastq
  |-- notes.typ
  |-- out
  |   `-- # script outputs including filtered files, figures, and intermediate steps
  `-- src 
      `-- # qc and analysis scripts
```

# Code
## Environment setup
Use `environment.yaml` to setup your python environment with conda

```bash
conda env create -n environment.yaml
conda activate metagenome_profiler
```

## Install
Install this package directly from GitHub with 
``` bash
# # activate conda environment
# conda activate metagenome_profiler
pip install git+https://github.com/canderson318/metagenome_profiler.git@main
```

# Usage
The following example is available at [src/main.py](src/main.py).

1. Imports and file paths

```python
import sys
import time
from Bio import SeqIO
from pathlib import Path

# –– Paths
in_dir = Path("in/")
out_dir = Path("out/")

# path to output plots
fig_path = out_dir / 'figs'

# –– Data files
# Path to reference genome
ref_file = in_dir / "bat_genome/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
# Path to sample reads
samp_file = in_dir / "SRR12464727.fastq"

```

2. QC: Filter fastq reads
```python

# Select which scaffolds to index and reference against. here using the ncbi headers file with each chromosome's headers
genome_headers = [l.split(" ")[0].lstrip(">") 
                  for l in iter(open(ref_file.parent.parent / "headers.txt"))]

#  Filter fastq reads 
from src.lib.qc.filter_reads import filter_reads
filtered_read_path = filter_reads(samp_file, ref_file, 
                                  out_dir = out_dir, 
                                  K = 31,
                                  out_file = "filtered.fasta", 
                                  aln_thresh= .5,
                                  which_host_scaffolds = set(genome_headers)
                                  )

```

3. Classify reads
``` python

# path to kraken database
viral_db_path = in_dir / "kraken_2_viral_db"

# path to output results
result_dir = out_dir / "kraken"

# query reads against kraken viral db and load the formatted result
from src.analysis.classify_reads import classify_reads
result  = classify_reads(filtered_read_path, viral_db_path,result_dir, force = True, ncores = 5)
```

4. Analysis
```python

# –––– Look at Kraken classifications
from src.analysis.result_analysis import result_overview
result_overview(result,fig_path)

# –––– Look at abundances
from src.analysis.result_analysis import result_abundances
top_species_by_prop = result_abundances(result, fig_path)

# –––– Look for pattern of string similarity among same species
from src.analysis.result_analysis import read_str_distance_analysis
read_str_dist_MDS = read_str_distance_analysis(result, top_species_by_prop, out_dir, fig_path)

```


# Specific Libraries

`pandas`

- For data manipulation and csv I/O.

`biopython` 

- For fasta/fastq I/O and handling with `SeqIO` module. 

`skbio`

- For permanova test.

`seaborn` and `matplotlib.pyplot`

- For data visualization. 

`Levenshtein`

- For read to read edit distance calculation.

`sklearn.manifold`

- For read edit distance matrix dimension reduction with `MDS` module.

`kraken2`

- For viral read classification against viral reference database.