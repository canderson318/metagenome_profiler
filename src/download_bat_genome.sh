#!/bin/bash

echo "Downloading Bat Genome..."
conda run -n metagenome_profiler -- datasets download genome accession GCF_004115265.2 --include genome
echo "Done."

unzip -n ncbi_dataset.zip

rsync -av --update ncbi_dataset/ in/

mv in/data in/bat_genome/

rm -r ncbi_dataset*
rm md5sum.txt
