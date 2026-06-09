#!/bin/bash

total_reads=$(wc -l < out/phreds.txt)
host_removed=$(wc -l < out/reads_from_host.inds)
# (reads with avg phred < 20)
phred_removed=$(wc -l < out/phred_hits.txt)
clean=$(grep -c '^>' out/filtered.fasta)

conda run -n basic python3 -c "
total = $total_reads
host_removed = $host_removed
after_host = total - host_removed
phred_removed = $phred_removed
clean = $clean       

print(f'total reads         = {total:>12,}  (wc -l out/phreds.txt)')
print(f'host removed        = {host_removed:>12,}  (wc -l out/reads_from_host.inds)')
print(f'after host filter   = {after_host:>12,}  = {total} - {host_removed}')
print(f'phred removed       = {phred_removed:>12,}  (wc -l out/phred_hits.txt)')
print(f'naive clean         = {after_host - phred_removed:>12,}  = {after_host} - {phred_removed}')
print(f'actual clean        = {clean:>12,}  (grep -c \">\" out/filtered.fasta)')
overlap                     = after_host - phred_removed - clean
print(f'overlap (host+phred)= {overlap:>12,}  reads flagged by both filters')
print()
print(f'host %              = {host_removed/total*100:.4f}%')
print(f'phred %             = {phred_removed/total*100:.4f}%')
print(f'total removed %     = {(total-clean)/total*100:.4f}%')
"