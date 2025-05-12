# Best Thresholds Plots

This folder contains a series of plots showing the relationship between the **E-value threshold** and the corresponding **Matthews Correlation Coefficient (MCC)** for different test configurations.  
The goal of these plots is to visualize how the threshold selection affects model performance.

Each file name encodes the following information:

- `s1` or `s2`: Set 1 or Set 2
- `seq` or `str`: Sequence-based or Structure-based HMM
- `full` or `sing`: : Using full sequence E-value or best single domain E-value
- The **red dot** in each plot marks the threshold that gives the highest MCC for that combination of dataset, alignment type, and E-value type.
---

### File Descriptions

- `gr-s1-seq-full.png`  
  MCC vs. E-value threshold for **Set 1**, **sequence-based HMM**, using **full sequence** E-values.

- `gr-s1-seq-sing.png`  
  Same as above, but using **best domain** E-values.

- `gr-s1-str-full.png`  
  MCC vs. E-value threshold for **Set 1**, **structure-based HMM**, using **full sequence** E-values.

- `gr-s1-str-sing.png`  
  Same as above, but using **best domain** E-values.

- `gr-s2-seq-full.png`  
  MCC vs. E-value threshold for **Set 2**, **sequence-based HMM**, using **full sequence** E-values.

- `gr-s2-seq-sing.png`  
  Same as above, but using **best domain** E-values.

- `gr-s2-str-full.png`  
  MCC vs. E-value threshold for **Set 2**, **structure-based HMM**, using **full sequence** E-values.

- `gr-s2-str-sing.png`  
  Same as above, but using **best domain** E-values.

---

Each plot is useful to identify the **optimal threshold** that maximizes MCC, guiding the choice of E-value cutoff during classification.

> **Note**: These plots are generated from the first phase of analysis using **Pfam ID PF00014**, not the InterPro-based evaluation.
