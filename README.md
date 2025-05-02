# Building a Profile Hidden Markov Model for the Kunitz-Type Protease Inhibitor Domain

## Required Packages

The following packages are required:

- **CD-HIT**  
  **Purpose**: Clustering and redundancy reduction of protein sequences.  
  **Install via conda**:
  ```bash
  conda install -c bioconda cd-hit
  ```

- **HMMER**  
  **Purpose**: Building and searching Hidden Markov Models (HMMs) for protein domain detection.  
  **Install via conda**:
  ```bash
  conda install -c bioconda hmmer
  ```

- **BLAST+**  
  **Purpose**: Performing protein sequence similarity searches using `blastp`.  
  **Install via conda**:
  ```bash
  conda install -c bioconda blast
  ```

---

## Repository Contents

This repository contains:

- `rcsb_pdb_custom_report_20250410062557.csv` [see Additional Notes](#-contact-and-additional-notes)  
  A custom report downloaded from the PDB using the following query:

  > **QUERY:**  
  > Data Collection Resolution ≤ 3.5  
  > AND (Identifier = "PF00014" AND Annotation Type = "Pfam")  
  > AND Polymer Entity Sequence Length between 45 and 80

- `pdb_kunitz_rp.ali`  
  An initially empty file. You will paste into this file the multi-structure alignment obtained from PDBeFold.

- `all_kunitz.fasta` [see Additional Notes](#-contact-and-additional-notes)  
  A FASTA file containing all Kunitz proteins (human and non-human).

- `script_recover_representative_kunitz.sh`  
  A script that extracts representative Kunitz domain PDB IDs and writes them to `pdb_efold_ids.txt`.

- `create_hmm_build.sh`  
  A script that builds the HMM model and evaluates its performance using the datasets.

---

## What to Do

### 0. Download the Swiss-Prot FASTA file

Visit [UniProt](https://www.uniprot.org/) and download all protein sequences annotated in the Swiss-Prot database. Insert the file into the folder with the scripts and recall it `uniprot_sprot.fasta`.  
This file will be used to extract non-Kunitz protein sequences.

---

### 1. Install the required packages

Install the packages listed in the [Required Packages](#-required-packages) section using `conda`.

---

### 2. Run the representative ID extraction script

```bash
bash script_recover_representative_kunitz.sh
```

This will create the file `pdb_efold_ids.txt`.

---

### 3. Perform structure-based multiple alignment

1. Go to the [PDBeFold Multi Alignment Tool](https://www.ebi.ac.uk/msd-srv/ssm/)
2. Choose:
   - **Mode**: Multiple
   - **Source**: List of PDB codes
3. Upload the file `pdb_efold_ids.txt`
4. Click **Download FASTA Alignment**
5. Paste the downloaded content into the file `pdb_kunitz_rp.ali`

---

### 4. Build the HMM and evaluate performance

Run the following script:

```bash
bash create_hmm_build.sh
```

This will:

- Build an HMM model from the structural alignment
- Generate random subsets of positive and negative sequences, which will be used to create test sets  
- Run `hmmsearch` on those sets
- Automatically identify the best E-value thresholds
- Calculate performance metrics (q2, MCC, etc.)
- Identify false positives and false negatives
- Save all results into `hmm_results.txt`

---

## Output

- **`hmm_results.txt`**  
  This file contains:
  - The best E-value thresholds selected by maximizing the Matthews Correlation Coefficient (MCC)
  - Performance metrics for each test set and overall, calculated using the E-value that yielded the highest MCC, based on either full sequence or best single domain evaluations
  - Lists of false positives and false negatives

---

## Contact and Additional Notes

You are not required to use the provided `rcsb_pdb_custom_report_20250410062557.csv` or `all_kunitz.fasta` files.

- To generate your own **PDB custom report**, visit the [RCSB PDB Advanced Search](https://www.rcsb.org/search/advanced) and customize the filters to suit your needs (e.g., domain = PF00014, resolution limits, sequence length).
- To build your own **all_kunitz.fasta**, go to [UniProt](https://www.uniprot.org/), search for all proteins annotated with the Pfam domain `PF00014`, and download the resulting sequences in FASTA format.

---

For any questions, suggestions, or contributions, feel free to open an issue or contact the maintainer: **Marco Cuscunà**.
