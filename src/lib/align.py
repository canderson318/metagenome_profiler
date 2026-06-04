from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import os
import numpy as np 
import pandas as pd
from itertools import islice
from  multiprocessing import Pool
from functools import partial

wd = Path.home()/ "Documents/school/prelims/day3"
os.chdir(wd)

#\\\\
#\\\\
# Load Data
#\\\\
#\\\\

in_dir = Path("in/")
out_dir = Path("out/")


ref_file = in_dir / "data/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"

samp_file = in_dir / "SRR12464727.fastq"


# #\\\\
# #\\\\
# # 
# #\\\\
# #\\\\

from src.lib.index_host import index_host

K = 31
ref_index_file_path, _ = index_host(K=K, ref_file=ref_file, which_scaffolds = [0], force = False)


#\\\\
#\\\\
# classify reads
#\\\\
#\\\\
    
from src.lib.align_reads import align_reads
    
hits = align_reads(samp_file, ref_index_file_path, K = K, test = True)