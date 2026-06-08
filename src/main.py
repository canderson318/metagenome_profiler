from Bio import SeqIO
from pathlib import Path

# \\
# –– Paths
# \\
in_dir = Path("in/")
out_dir = Path("out/")

fig_path = out_dir / 'figs'
fig_path.mkdir(exist_ok= True)

# Data files
ref_file = in_dir / "bat_genome/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
samp_file = in_dir / "SRR12464727.fastq"

# \\\\
# \\\\
# –––––––––––––––––––––––––– QC: filter reads ––––––––––––––––––––––––––
# \\\\
# \\\\

# select chromosome headers
genome_headers = [l.split(" ")[0].lstrip(">") 
                  for l in iter(open(ref_file.parent.parent / "headers.txt"))]
# genome_headers_set = set(genome_headers)
# seqs = {rec.id: rec.seq for rec in SeqIO.parse(ref_file,'fasta') if rec.id in genome_headers_set}

from src.lib.qc.filter_reads import filter_reads
filtered_read_path = filter_reads(samp_file, ref_file, out_dir = out_dir, K = 31,
                                  out_file = "filtered.fasta", 
                                  aln_thresh= .5,
                                  force_phred = False, # phred.txt does not change depending on other parameters, phred_hits.txt changes based on threshold
                                  force_index=False, # only changes with K
                                  force_aln = True, # changes based on K and thresholds
                                #   which_host_scaffolds = list(genome_headers)[0]
                                  which_host_scaffolds = set(genome_headers[])
                                  )

# \\\\
# \\\\
# ––––––––––––––––––––––––––– reference viral genome –––––––––––––––––––––––––––––––––
# \\\\
# \\\\

# filtered_read_path = out_dir / "filtered.fasta"
viral_db_path = in_dir / "kraken_2_viral_db"
result_dir = out_dir / "kraken"

from src.analysis.classify_reads import classify_reads

# query reads against kraken viral db and load kraken result
result  = classify_reads(filtered_read_path, viral_db_path,result_dir, 
                          force = True,
                          ncores = 5)

# \\\\
# \\\\
# –––––––––––––––––––––––––––––– Analyze result –––––––––––––––––––––––––––––––––––––––
# \\\\
# \\\\

# \\\\
# –––– Look at Kraken classifications
# \\\\
from src.analysis.result_analysis import result_overview
result_overview(result,fig_path)

# \\\\
# –––– Look at abundances
# \\\\
from src.analysis.result_analysis import result_abundances
top_species_by_prop = result_abundances(result, fig_path)

# \\\\
# –––– Look for pattern of string similarity among same species
# \\\\
from src.analysis.result_analysis import read_str_distance_analysis
read_str_dist_MDS = read_str_distance_analysis(result, top_species_by_prop,
                                               out_dir,fig_path, read_sample_N = 500)