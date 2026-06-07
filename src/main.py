from pathlib import Path
from Bio import SeqIO as SeqIO
import pickle as pkl
from src.lib.qc.filter_reads import filter_reads
from src.lib.assembly.DBG import DBG
from src.lib.analysis.classify_reads import classify_reads


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

filtered_read_path = filter_reads(samp_file, ref_file, out_dir = out_dir, K = 31, out_file = "host_filtered.fasta", 
                                  test = False)
# ^^out/host_filtered.fasta

# \\\\
# \\\\
# assemble contigs
# \\\\
# \\\\

# G = DBG(SeqIO.parse(filtered_read_path, 'fasta'), K = 85)
# dbg_dir = out_dir / "DBG"
# dbg_dir.mkdir(exist_ok=True)

# with open(dbg_dir / "dbg.pkl", 'wb') as fout: 
#     pkl.dump(G,fout)
    
# contigs = G.make_contigs(True, join_sub_contigs=True)
# with open(dbg_dir / "dbg_contigs.pkl", 'wb') as fout: 
#     pkl.dump(G,fout)

# with open(out_dir / "contigs.fasta" , 'wt') as fout:
#     for i,seq in enumerate(contigs):
#         fout.write(f">contig_{i}\n{seq}\n")

# covs = np.array([ sum(c)/len(c) for s,n,c in contigs])
# lens = np.array([len(s) for s,n,c in contigs])
# import numpy as np
# import seaborn as sns
# sns.histplot((covs/covs.max()) * 100, kde = True)
# sns.histplot((lens/lens.max()) * 100, kde = True)
# p = 70
# len_perc = np.percentile(lens, p)
# cov_perc = np.percentile(covs, p)
# top_contigs = [s for s,n,c in contigs if (sum(c)/len(c) > cov_perc) and (len(s) > len_perc)]


# \\\\
# \\\\
# reference viral genome
# \\\\
# \\\\


results  = classify_reads(filtered_read_path)