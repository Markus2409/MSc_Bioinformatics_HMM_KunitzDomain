## `graphics_and_tables` folder

This folder contains a collection of plots and tables summarizing key results from the HMM performance evaluations. Many of these have been used in the final report, while others are provided here as additional material for completeness, even if not explicitly included in the report.

### Contents:
- `Best thresholds plots/`  
  A subfolder containing plots of **MCC vs. E-value threshold** curves. These plots help visualize how the MCC varies across different thresholds for each test set, HMM model (sequence-based or structure-based), and E-value type (full sequence or best domain).

- `best_threshold_seqali.pdf`  
  A LaTeX table listing the best E-value thresholds (those that maximize MCC) for each test set using the **sequence-based HMM**, along with the corresponding MCC values.

- `best_threshold_strali.pdf`  
  Same as above, but using the **structure-based HMM** derived from PDBeFold alignment.

- `confusion_matrices_seqali.pdf`  
  Confusion matrices for all test conditions related to the **sequence-based model**, including cross-validation and overall evaluations.

- `confusion_matrices_strali.pdf`  
  Confusion matrices for the **structure-based model**, in the same format.

- `overall_performances.pdf`  
  A comprehensive table listing all **false positives** and **false negatives** for each model (sequence- and structure-based) across both E-value types (full sequence and best domain), for both the **standard Pfam-based** and **InterPro-augmented** analyses. Also includes the final confusion matrices and average MCC values obtained on the merged dataset (`temp_overall.class`), using thresholds selected from both Set 1 and Set 2.

---

> **Note**:  
> The plots and tables in `best_threshold_seqali.pdf`, `best_threshold_strali.pdf`, and the folder `Best thresholds plots` all refer exclusively to the initial analysis based on the Pfam ID `PF00014`.  
> For simplicity, results from the InterPro-based extension were not reformatted into polished LaTeX tables or plots.  
> Instead, all key insights and false positive/negative counts for the InterPro analysis are already summarized in `overall_performances.pdf`.
