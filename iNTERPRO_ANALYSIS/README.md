## iNTERPRO_ANALYSIS Folder

This folder contains all files and datasets used to repeat the main analysis after identifying a crucial issue related to missing annotations in the original Pfam-based domain filtering.

During the first analysis (based solely on Pfam ID `PF00014`), a false positive predicted by the model was revealed to actually be a true Kunitz domain protein — it was simply missing the `PF00014` annotation. Further inspection showed that it was correctly annotated with the corresponding InterPro domain (`IPR036880`), indicating a **cross-reference inconsistency** between Pfam and InterPro.

To address this, a second query was issued on the PDB:

> **Updated Query Used**:  
> `Data Collection Resolution <= 3.5 AND ((Identifier = "PF00014" AND Annotation Type = "Pfam") OR (Lineage Identifier = "IPR036880" AND Annotation Type = "InterPro")) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45`

This allowed the inclusion of all relevant Kunitz domain proteins, whether annotated through Pfam or InterPro, leading to a broader and more reliable training and testing dataset.

### Contents

This folder includes:

- The updated `rcsb_pdb_custom_report_20250410062557.csv` PDB custom report file containing both Pfam and InterPro annotated Kunitz entries  
- An extended version of `all_kunitz.fasta` including proteins annotated **only** via InterPro  
- A revised version of `create_testing_sets.sh`, which redefines the positive and negative sets based on the expanded pool of Kunitz proteins  

> **Note**: The rest of the pipeline and associated scripts remain unchanged. Representative selection, alignments, and HMM construction still follow the procedure described in the main workflow.  
> The results obtained from this extended analysis — including false positives, false negatives, and overall performance — are documented in the  `OUTPUT_EXAMPLES` folder.
