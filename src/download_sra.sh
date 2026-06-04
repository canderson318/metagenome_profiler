#!/bin/bash

## Run with:
## ./download_sra.sh

ACCESSIONS=("SRR12464727" "SRR12432009")
OUT_DIR="in/"
a_download() {
    ACC=$1
    echo "Downloading $ACC..."
    # prefetch "$ACC"
    # echo "Converting $ACC..."
    # fasterq-dump --split-files --progress --outdir in/ "$ACC"
    fastq-dump \
            --split-3 \
            --skip-technical \
            --clip \
            --qual-filter-1 \
            -O "$OUT_DIR" \
            "$ACC"  
    echo "Done: $ACC"
}

for ACC in "${ACCESSIONS[@]}"; do 
    a_download $ACC || { echo "Failure"; exit 1; }
done

# export -f a_download
# parallel a_download ::: "${ACCESSIONS[@]}"


echo "All downloads complete."


# day3$ ./src/download_sra.sh
# Downloading SRR12464727...
# Rejected 39335609 READS because of filtering out non-biological READS
# Rejected 78117 READS because of Quality-Filtering
# Read 39335609 spots for SRR12464727
# Written 39257492 spots for SRR12464727
# Done: SRR12464727
# Downloading SRR12432009...

