#!/bin/bash

echo "Downloading Bat Genome..."
conda run -n basic -- datasets download genome accession GCF_004115265.2 --include genome
echo "Done."

unzip -n ncbi_dataset.zip

rsync -av --update ncbi_dataset/ in/

rm -r ncbi_dataset*

rm md5sum.txt
