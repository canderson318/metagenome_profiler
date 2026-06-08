import subprocess as sp
from pathlib import Path
import pandas as pd
from typing import Union

def classify_reads(fasta_path: Union[str,Path],db_path:Union[str,Path],out_dir:Union[str,Path], force = False, ncores = 1):
    """
    Queries fasta reads against viral genome database using kraken2. 
    Returns a pandas DataFrame of classified reads including their fasta accession number, read index, 
    and the location of alignment with the viral taxon.
    
    :param fasta_path: Path to fasta file.
    :return: pandas DataFrame `result`
    """
    fasta_path = Path(fasta_path)
    
    if not fasta_path.exists():
        raise FileNotFoundError

    if not db_path.exists():
        raise FileNotFoundError

    if not out_dir.exists():
        out_dir.mkdir(exist_ok=True)

    if not (out_dir / "report.txt").exists() or force:
        print("Querying viral reference database...")
        sp.run([
            "kraken2", "--db", str(db_path),
            "--threads", str(ncores),
            "--output", str(out_dir / "results.txt"),
            "--report", str(out_dir / "report.txt"),
            str(fasta_path)
        ], stderr=sp.DEVNULL, stdout = sp.DEVNULL)  
        print("Done.")
    else:
        print(f"Viral taxon classification already exists at {str(out_dir)} so using those.\nUse `force` to query again.")
    
    print("Loading results...")
    
    # columns: pct_reads, reads_covered, reads_direct, rank, taxid, name
    ### rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies
    report = pd.read_csv(out_dir / "report.txt", sep="\t", header=None, 
                    names=["pct","reads_cov","reads_direct","rank","taxid","name"])
    result = pd.read_csv(out_dir / "results.txt", sep = '\t', header = None, 
                        names = ["classified", "readid", "taxid", "length", "loc"])
    print("Done.")
    
    # clean ws
    report = report.apply(lambda col: col.str.strip() if col.dtype == object else col)
    result = result.apply(lambda col: col.str.strip() if col.dtype == object else col)
    
    # set index to readid
    result.set_index("readid",inplace = True, drop =  False)
    
    # make accession column
    result[['acc','idx']] = [x.split('.') for x in result['readid']]
    result['idx'] = result['idx'].astype(int)
    
    # assign taxon name for each read
    result["name"] = result['taxid'].map(
        {taxid: report.loc[report['taxid'] == taxid]["name"].iloc[0] 
         for taxid in report['taxid'].unique()}
    )
    
    result = result.merge(report[['name','rank','taxid']], how = 'left')
    
    return result    
