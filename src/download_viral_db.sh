#!/bin/bash

echo Downloading viral database...
wget --progress=bar:force https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz 

DIR="in/kraken_2_viral_db"
mkdir -p "$DIR"
tar -xzf k2_viral*.tar.gz -C "$DIR"
mv -f k2_viral*.tar.gz in/