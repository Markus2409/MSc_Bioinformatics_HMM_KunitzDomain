#!/usr/bin/bash
awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s",toupper($1)}}' pdb_kunitz_rp.ali | sed s/PDB:// | tail -n +2 > pdb_kunitz_rp_clean.fasta #with this command we converted the file aligned from pdbefold in a format that we can use to create hmm model
hmmbuild  pdb_kunitz_rp_clean.hmm pdb_kunitz_rp_clean.fasta #this command creates the hmm model with the sequences that you extracted from pdb!
makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta #questo crea il database di tutte le kunitz protein.
blastp -query pdb_kunitz_rp.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7 #in questo modo abbiamo il blast delle 23 sequenze ovver gli allineamenti con tutte le sequenze all kunitz. 
grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > pdb_idlist.txt #cosi pesco tutti quegli id che hanno identita alta alle proteine che ho usato per la costruzione di hmm model cosi poi potro eliminarle.
grep ">" all_kunitz.fasta | cut -d "|" -f 2 > all_kunitz.id #con questo mi crea un file con tutti gli id delle kunitz esistenti

echo 'creo il test file set1 e 2 dei positivi....' 
comm -23 <(sort all_kunitz.id) <(sort pdb_idlist.txt) > to_keep.ids #effettuiamo la rimozione e salviamo quelli rimanenti in to keep ids. Questi saranno le kunitz protein che useremo poi come train e test model.
sort -r to_keep.ids > random_ok_kunitz.id #randomizzo
head -n 184 random_ok_kunitz.id > pos_1.id
tail -n 184 random_ok_kunitz.id > pos_2.id
python3 get_seq.py pos_1.id uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py pos_2.id uniprot_sprot.fasta > pos_2.fasta
 # in questo modo io prendo le sequenze dei positivi dal fasta file

echo 'creo il test file set1 e 2  dei negativi....' 
grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 > sp.id #questo prende tutti gli id di swiss prot esistenti
comm -23 <(sort sp.id) <(sort all_kunitz.id) > sp_negs.ids # questo cancella tutte le kunitz proteins
sort -r sp_negs.ids >random_sp_negs.id #randomizzo
head -n 286286 random_sp_negs.id > neg_1.id
tail  -n 286286 random_sp_negs.id > neg_2.id
python3 get_seq.py neg_1.id uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_2.id uniprot_sprot.fasta > neg_2.fasta
echo 'fasta files creati'

#ora hai creato i set files è ora di andare a creare dei tabular report usando l'hmmmodel creato utilizzandoi testing sets
#POS_1
hmmsearch -Z 1000  --max --tblout pos_1.out pdb_kunitz_rp_clean.hmm pos_1.fasta # In this way we create the tabular output of our results.


#POS_2
hmmsearch -Z 1000  --max --tblout pos_2.out pdb_kunitz_rp_clean.hmm pos_2.fasta # In this way we create the tabular output of our results.


#NEG_1
hmmsearch -Z 1000  --max --tblout neg_1.out pdb_kunitz_rp_clean.hmm neg_1.fasta # In this way we create the tabular output of our results.


#NEG_2
hmmsearch -Z 1000  --max --tblout neg_2.out pdb_kunitz_rp_clean.hmm neg_2.fasta # In this way we create the tabular output of our results.


echo 'number of sequences in pos_1.out file:' #in this way you can check if all the sequences are present or not in .out file!
grep -v "^#" pos_1.out | cut -d " " -f1 | wc
echo 'number of sequences pos_1.ids file:'
grep -v "^#" pos_1.id | cut -d " " -f1 | wc
echo 'number of sequences in pos_2.out file:' #in this way you can check if all the sequences are present or not in .out file!
grep -v "^#" pos_2.out | cut -d " " -f1 | wc
echo 'number of sequences pos_2.ids file:'
grep -v "^#" pos_2.id | cut -d " " -f1 | wc
echo 'number of sequences in neg_1.out file:' #in this way you can check if all the sequences are present or not in .out file!
grep -v "^#" neg_1.out | cut -d " " -f1 | wc
echo 'number of sequences neg_1.ids file:'
grep -v "^#" neg_1.id | cut -d " " -f1 | wc
echo 'number of sequences in neg_2.out file:' #in this way you can check if all the sequences are present or not in .out file!
grep -v "^#" neg_2.out | cut -d " " -f1 | wc
echo 'number of sequences neg_2.ids file:'
grep -v "^#" neg_2.id | cut -d " " -f1 | wc

#now we recover the e-value for each of sequence  with the respective ids and the 1 for positives and 0 for negatives $5 is the e-value computed for the full sequence, instead the $8 is the e-value computed for the best 1 domain
#the positive ones are the matches if you want more matches you have to edit the threshold of hmmer that is set by default to ten. 
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}'>pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}'>pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}'>neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}'>neg_2.class

#in questo modo aggiungo i veri negativi, ovvero quelli che che non erano presenti in neg_1 e neg_2 class, dandogli un evalue fittizio, dato che non erano stati trovati dall'hmm
comm -23 <(sort neg_1.id) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}'>>neg_1.class
comm -23 <(sort neg_2.id) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}'>>neg_2.class

#adesso mettiamo insieme pos_1 e neg_1 e pos_2 e neg_2 per creare quindi set_1 e set_2, questi ci serviranno per misurare le performances dell'algoritmo
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
#per l'overall
cat set_1.class set_2.class > temp_overall.class

#troviamo la miglior threshold per ciascun set 
#SET 1
#per FULL SEQUENCE
set1_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)


#per BEST 1 DOMAIN
set1_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_1.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
#SET 2
#per FULL SEQUENCE
set2_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
#per BEST 1 DOMAIN
set2_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py set_2.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

#OVERALL
#per FULL SEQUENCE
overall_best_evalue_full_seq=$(
    for i in $(seq 1 9); do
        python3 performance.py temp_overall.class 1e-"$i"
    done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)
#per BEST 1 DOMAIN
overall_best_evalue_one_domain=$(
    for i in $(seq 1 9); do
        python3 performance.py temp_overall.class 1e-"$i"
    done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print $2}'
)

#calcoliamo le performances con il performance.py impostando la threshold a 1e-5
echo -e "\nPERFORMANCES SET_1 FULL SEQUENCE"  > hmm_results.txt
python3 performance.py set_1.class "$set1_best_evalue_full_seq" 1  >> hmm_results.txt
echo -e "\nPERFORMANCES SET_1 ONE DOMAIN" >> hmm_results.txt
python3 performance.py set_1.class "$set1_best_evalue_one_domain" 2 >> hmm_results.txt
echo -e "\nPERFORMANCES SET_2 FULL SEQUENCE" >> hmm_results.txt
python3 performance.py set_2.class $set2_best_evalue_full_seq 1 >> hmm_results.txt
echo -e "\nPERFORMANCES SET_2 ONE DOMAIN" >> hmm_results.txt
python3 performance.py set_2.class $set2_best_evalue_one_domain 2 >> hmm_results.txt
echo -e "\nOVERALL PERFORMANCES FULL SEQUENCE" >> hmm_results.txt
python3 performance.py temp_overall.class $overall_best_evalue_full_seq 1 >> hmm_results.txt
echo -e "\nOVERALL PERFORMANCES ONE DOMAIN" >> hmm_results.txt
python3 performance.py temp_overall.class $overall_best_evalue_one_domain 2 >> hmm_results.txt

#now we want understand which is that false negative, that are values that are positives (kunitz protein) but are upside the threshold and are inserted in negatives ones! in our case we can use the best evalue full seq that we have found for respective set
echo -e "\nFalse positives for set 1:" >>hmm_results.txt
awk -v num="$set1_best_evalue_full_seq" '$3 > num' pos_1.class | sort -grk 3 >>hmm_results.txt
echo -e "\nFalse positives for set 2:" >>hmm_results.txt
awk -v num="$set2_best_evalue_full_seq" '$3 > num' pos_2.class | sort -grk 3 >>hmm_results.txt

#false positives are values that are negatives but are into the positives ones because they have a evalue that is below the threshold! 
echo -e "\nFalse negatives for set 1:" >>hmm_results.txt
awk -v num=$set1_best_evalue_full_seq '$3 < num' neg_1.class | sort -grk 3 >> hmm_results.txt
echo -e "\nFalse positives for set 2:" >>hmm_results.txt
awk -v num=$set2_best_evalue_full_seq '$3 < num' neg_2.class | sort -grk 3 >> hmm_results.txt

