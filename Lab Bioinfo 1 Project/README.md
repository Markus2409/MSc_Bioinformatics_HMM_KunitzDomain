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

- `rcsb_pdb_custom_report_20250410062557.csv`  
  A custom report downloaded from the PDB using the following query:

  > **QUERY:**  
  > Data Collection Resolution ≤ 3.5  
  > AND (Identifier = "PF00014" AND Annotation Type = "Pfam")  
  > AND Polymer Entity Sequence Length between 45 and 80

- `pdb_kunitz_rp.ali`  
  An initially empty file. You will paste into this file the multi-structure alignment obtained from PDBeFold.

- `all_kunitz.fasta`  
  A FASTA file containing all Kunitz proteins (human and non-human).

- `uniprot_sprot.fasta`  
  A FASTA file containing the SwissProt database, used to extract non-Kunitz proteins.

- `script_recover_representative_kunitz.sh`  
  A script that extracts representative Kunitz domain PDB IDs and writes them to `pdb_efold_ids.txt`.

- `create_hmm_build.bash`  
  A script that builds the HMM model and evaluates its performance using the datasets.

---

## What to Do

### 1. Install the required packages

Install the packages listed in the [Required Packages](#-required-packages) section using conda.

---

### 2. Run the representative ID extraction script

```bash
bash script_recover_representative_kunitz.sh
```

This creates the file `pdb_efold_ids.txt`.

---

### 3. Perform structure-based multiple alignment

1. Go to [PDBeFold Multi Alignment Tool](https://www.ebi.ac.uk/msd-srv/ssm/)
2. Choose:
   - **Mode**: Multiple
   - **Source**: List of PDB codes
3. Upload the `pdb_efold_ids.txt` file.
4. Click "Download FASTA Alignment"
5. Paste the content into the file `pdb_kunitz_rp.ali`.

---

### 4. Build the HMM and run evaluation

Execute the script:

```bash
bash create_hmm_build.bash
```

This will:

- Build an HMM model from the structural alignment
- Generate training and test datasets from positive and negative sets
- Run `hmmsearch` on those sets
- Automatically identify the best E-value thresholds
- Calculate performance metrics (accuracy, q2, MCC, etc.)
- Detect false positives and false negatives
- Save the output to `hmm_results.txt`

---

## Output

- `hmm_results.txt`:  
  Contains:
  - Best E-value thresholds for full sequence and single domain predictions
  - Performance metrics for each test set and overall
  - Lists of false positives and false negatives

---

## Contact

For any questions, suggestions or contributions, feel free to open an issue or contact the maintainer Marco Cuscunà.
