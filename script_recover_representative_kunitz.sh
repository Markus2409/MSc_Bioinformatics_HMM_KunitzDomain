cat rcsb_pdb_custom_report_20250410062557.csv |tr -d '"' | awk -F ',' '{if (length($2)>0) {name=$2}; print name ,$3,$4,$5}' | grep PF00014| awk '{print">"$1"_"$3;print$2}' > pdb_kunitz.fasta #a partire dal file pdb report ti scarica le sequenze dei domini pf00014 di tutte le proteine. 
cd-hit -i pdb_kunitz.fasta -o pdb_kunitz.clstr -c 0.9 #questo ti clusterizza in base alla similarità. -c 0.9 hanno una similarita di 0.9
clstr2txt.pl pdb_kunitz.clstr.clstr | awk '{if ($5==1) print $1}' > pdb_kunitz_rp.ids #questo ti prende le piu rappresentative per ciascun cluster e te le mette in un file a parte

for i in `cat pdb_kunitz_rp.ids`; do
  grep -A 1 "^>$i" pdb_kunitz.fasta | tail -n 2 >> pdb_kunitz_rp.fasta
done 
# with this command i create a file that contain the sequences of ids representative that we got from cd-hit
#at this point you have to check if there are proteins that has sequences too long. In this case delete it. otherwise launch pdbefold!
grep ">" pdb_kunitz_rp.fasta |tr -d ">" | tr "_" ":" > pdb_efold_ids.txt #generate the ids that you have to insert in pdbefold to do the msa of structure