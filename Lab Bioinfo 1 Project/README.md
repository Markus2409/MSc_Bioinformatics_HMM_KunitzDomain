# Building a Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain
### Package that are needed:
cd-hit - a cosa serve
comando install via conda
hmmer - a cosa serve
comando install via conda
blast - a cosa serve
comando install via conda


### What is contened into the repository?
The repository contains:
rcsb_pdb_custom_report_20250410062557.csv a custom reports downloaded from pdb using this query: 
 magari mettila in un quote: QUERY: Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
 pdb_kunitz_rp.ali  an empty fail that will contains the pdbefold multi structural alignments.
all_kunitz.fasta, a file that contains all the kunitz proteins (human and not human).
uniprot_sprot.fasta , a file that contains all sequences of the swissprot database, and it is useful to get from the all_kunitz, all the non_kunitz proteins sequences
the script script_recover_representative_kunitz.sh extract the ids of most rapresentative kunitz proteins and inserts them in a properly formatted pdb_efold_ids.txt that is useful to make the MultiAlignment structural via PDBe fold 
the script create_hmm_build.bash that creates the hmm model using the all_kunitz.fasta file, and the uniprot_sprot.fasta 

##What to do?
After installed the Package (insert the link to the paragraph packlage that are needed) you have to run the script script_recover_representative_kunitz.sh. After that go to the website fo pdb efold https://www.ebi.ac.uk/msd-srv/ssm/
from this select multiple and source: List of PDB codes, and load the file that the script created (pdb_efold_ids.txt).
After that you have to click download fasta alignment and copy all the contents in the empty file called pdb_kunitz_rp.ali 
Now you can run the creat_hmm_build.bash that will create for you the hmm model, and will test two random sets of positives and negatives taken from the positives and negatives found thanks to the hmm model. In this way you will have a file callend hmm_results.txt that will contains the best performances for both sets and the overall performance gained combining the two sets.Also will individuate the false positives and negatives.
