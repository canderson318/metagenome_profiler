from pathlib import Path

in_dir = Path("in/")
out_dir = Path("out/")


# Data files
ref_file = in_dir / "bat_genome/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
samp_file = in_dir / "SRR12464727.fastq"

# \\\\
# \\\\
# filter reads
# \\\\
# \\\\
from src.lib.qc.filter_reads import filter_reads
filtered_read_path = filter_reads(samp_file, ref_file, out_dir = out_dir, K = 31, out_file = "host_filtered.fasta", test = True)
# out/host_filtered.fasta

# \\\\
# \\\\
# assemble contigs
# \\\\
# \\\\

from src.lib.assembly.DBG import DBG
from Bio import SeqIO as SeqIO

filtered_read_path = out_dir / "host_filtered.fasta"

sum(1 for _ in SeqIO.parse(filtered_read_path, 'fasta'))
G = DBG(SeqIO.parse(filtered_read_path, 'fasta'), k = 31)


# \\\\
# \\\\
# reference viral genome
# \\\\
# \\\\


