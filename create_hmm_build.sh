#!/usr/bin/bash

# Convert the PDBefold alignment to FASTA format (uppercase, one-line), removing 'PDB:' prefix
awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s",toupper($1)}}' pdb_kunitz_rp.ali | sed s/PDB:// | tail -n +2 > pdb_kunitz_rp_clean.fasta

# Build the HMM model from the aligned sequences
hmmbuild pdb_kunitz_rp_clean.hmm pdb_kunitz_rp_clean.fasta

# Create a BLAST database from all Kunitz proteins
makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta

# Perform BLAST search of the 23 representative sequences against the full Kunitz dataset
blastp -query pdb_kunitz_rp.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7

# Extract IDs of sequences with high identity (≥95%) and alignment length ≥50 to remove them from the training/testing pool
grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > pdb_idlist.txt

# Extract all IDs from the full Kunitz dataset
grep ">" all_kunitz.fasta | cut -d "|" -f 2 > all_kunitz.id

echo 'Creating positive test files (set1 and set2)...'

# Remove the overlapping IDs and keep only the remaining Kunitz proteins
comm -23 <(sort all_kunitz.id) <(sort pdb_idlist.txt) > to_keep.ids
sort -r to_keep.ids > random_ok_kunitz.id
head -n 184 random_ok_kunitz.id > pos_1.id
tail -n 184 random_ok_kunitz.id > pos_2.id
python3 get_seq.py pos_1.id uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py pos_2.id uniprot_sprot.fasta > pos_2.fasta

echo 'Creating negative test files (set1 and set2)...'

# Extract all Swiss-Prot IDs
grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 > sp.id

# Remove the Kunitz proteins to get the negative candidates
comm -23 <(sort sp.id) <(sort all_kunitz.id) > sp_negs.ids
sort -r sp_negs.ids > random_sp_negs.id
head -n 286286 random_sp_negs.id > neg_1.id
tail -n 286286 random_sp_negs.id > neg_2.id
python3 get_seq.py neg_1.id uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_2.id uniprot_sprot.fasta > neg_2.fasta

echo 'FASTA files created.'

# Run hmmsearch on all positive and negative FASTA files and create tabular output
hmmsearch -Z 1000 --max --tblout pos_1.out pdb_kunitz_rp_clean.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out pdb_kunitz_rp_clean.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out pdb_kunitz_rp_clean.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out pdb_kunitz_rp_clean.hmm neg_2.fasta

# Check if the number of sequences matches between input and hmmsearch output
echo -n 'Number of sequences in pos_1.out: ' >> hmm_results.txt
grep -v "^#" pos_1.out | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in pos_1.id: ' >> hmm_results.txt
grep -v "^#" pos_1.id | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in pos_2.out: ' >> hmm_results.txt
grep -v "^#" pos_2.out | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in pos_2.id: ' >> hmm_results.txt
grep -v "^#" pos_2.id | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in neg_1.out: ' >> hmm_results.txt
grep -v "^#" neg_1.out | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in neg_1.id: ' >> hmm_results.txt
grep -v "^#" neg_1.id | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in neg_2.out: ' >> hmm_results.txt
grep -v "^#" neg_2.out | cut -d " " -f1 | wc -l >> hmm_results.txt

echo -n 'Number of sequences in neg_2.id: ' >> hmm_results.txt
grep -v "^#" neg_2.id | cut -d " " -f1 | wc -l >> hmm_results.txt


# Extract E-values and build .class files with: ID, label (1/0), full sequence E-value ($5), and best domain E-value ($8)
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_2.class

# Add true negatives not detected by HMM (missing from the .out), with a fake high E-value (10.0)
comm -23 <(sort neg_1.id) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t"0"\t10.0\t10.0}' >> neg_1.class
comm -23 <(sort neg_2.id) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t"0"\t10.0\t10.0}' >> neg_2.class

# Combine positive and negative sets into final training/testing files
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
cat set_1.class set_2.class > temp_overall.class

# Automatically determine the best thresholds (based on highest MCC)
# SET 1
set1_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set1_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

# SET 2
set2_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
set2_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

# OVERALL
overall_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py temp_overall.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
overall_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py temp_overall.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

# Run final performance evaluations using the best thresholds
echo -e "\nPERFORMANCES SET_1 FULL SEQUENCE" > hmm_results.txt
python3 performance.py set_1.class "$set1_best_evalue_full_seq" 1 >> hmm_results.txt
echo -e "\nPERFORMANCES SET_1 ONE DOMAIN" >> hmm_results.txt
python3 performance.py set_1.class "$set1_best_evalue_one_domain" 2 >> hmm_results.txt
echo -e "\nPERFORMANCES SET_2 FULL SEQUENCE" >> hmm_results.txt
python3 performance.py set_2.class "$set2_best_evalue_full_seq" 1 >> hmm_results.txt
echo -e "\nPERFORMANCES SET_2 ONE DOMAIN" >> hmm_results.txt
python3 performance.py set_2.class "$set2_best_evalue_one_domain" 2 >> hmm_results.txt
echo -e "\nOVERALL PERFORMANCES FULL SEQUENCE" >> hmm_results.txt
python3 performance.py temp_overall.class "$overall_best_evalue_full_seq" 1 >> hmm_results.txt
echo -e "\nOVERALL PERFORMANCES ONE DOMAIN" >> hmm_results.txt
python3 performance.py temp_overall.class "$overall_best_evalue_one_domain" 2 >> hmm_results.txt

# Identify false positives: non-Kunitz proteins with E-value below threshold (misclassified as positive) 
echo -e "\nFalse positives for set 1:" >> hmm_results.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num' pos_1.class | sort -grk 3 >> hmm_results.txt
echo -e "\nFalse positives for set 2:" >> hmm_results.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num' pos_2.class | sort -grk 3 >> hmm_results.txt

# Identify false negatives: Kunitz proteins with E-value above threshold (misclassified as negative)
echo -e "\nFalse negatives for set 1:" >> hmm_results.txt
awk -v num="$set1_best_evalue_full_seq" '$3 < num' neg_1.class | sort -grk 3 >> hmm_results.txt
echo -e "\nFalse negatives for set 2:" >> hmm_results.txt
awk -v num="$set2_best_evalue_full_seq" '$3 < num' neg_2.class | sort -grk 3 >> hmm_results.txt
