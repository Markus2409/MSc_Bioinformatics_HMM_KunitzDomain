#!/usr/bin/bash

# Extract sequences of PF00014 domains from the PDB custom report
cat rcsb_pdb_custom_report_20250410062557.csv | tr -d '"' \
| awk -F ',' '{if (length($2)>0) {name=$2}; print name ,$3,$4,$5}' \
| grep PF00014 \
| awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta

# Cluster the sequences using CD-HIT at 90% identity threshold
cd-hit -i pdb_kunitz.fasta -o pdb_kunitz.clstr -c 0.9

# Extract the most representative ID from each cluster
clstr2txt.pl pdb_kunitz.clstr.clstr | awk '{if ($5==1) print $1}' > pdb_kunitz_rp.ids

# Retrieve the sequences of the representative IDs and store them in a new FASTA file
for i in `cat pdb_kunitz_rp.ids`; do
  grep -A 1 "^>$i" pdb_kunitz.fasta | tail -n 2 >> pdb_kunitz_rp.fasta
done

# This file now contains only the representative sequences selected from CD-HIT

# OPTIONAL: check manually for sequences that are too long and remove them before continuing

# Generate the PDBefold input format (convert IDs to the expected format: PDB:CHAIN)
grep ">" pdb_kunitz_rp.fasta | tr -d ">" | tr "_" ":" > pdb_efold_ids.txt
