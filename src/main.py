from pathlib import Path
from src.lib.qc.filter_reads import filter_reads
import numpy as np

in_dir = Path("in/")
out_dir = Path("out/")

# \\\\
# \\\\
# filter reads
# \\\\
# \\\\

ref_file = in_dir / "bat_genome/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
samp_file = in_dir / "SRR12464727.fastq"

filtered_read_path = filter_reads(samp_file, ref_file,out_dir = out_dir, K = 31, out_file = "host_filtered.fasta", test = True)

# \\\\
# \\\\
# assemble contigs
# \\\\
# \\\\





# \\\\
# \\\\
# reference viral genome
# \\\\
# \\\\


