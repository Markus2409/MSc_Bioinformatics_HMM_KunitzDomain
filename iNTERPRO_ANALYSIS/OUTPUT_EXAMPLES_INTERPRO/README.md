# OUTPUT_EXAMPLES_INTERPRO

This folder contains the full set of output files generated from the **second analysis**, which integrates Kunitz proteins annotated **only via InterPro (IPR036880)** in addition to those annotated via Pfam (PF00014), and which are used to perform the comparison in the final report.

It mirrors the structure and content of the main [`OUTPUT_EXAMPLES`](../OUTPUT_EXAMPLES) folder, but all files here reflect the **updated dataset** obtained using the broader, combined annotation query described in the [`INTERPRO_ANALYSIS`](../INTERPRO_ANALYSIS) section.

These outputs can be used independently to assess how the model performs when the training and testing sets are expanded to include proteins missed due to the Pfam-InterPro cross-reference inconsistency.

All output types — including `.fasta`, `.class`, `.hmm`, and performance summaries — are organized identically to those in the main pipeline.

> **Note**: The alignment and model-building procedures are unchanged. The only difference is in the composition of the input datasets, which now include previously excluded InterPro-only annotated Kunitz entries.
