# Acute Myeloid Leukemia (AML) Variant Analysis Project

> **Course**: BINF/STAT 5354  
> **Semester**: Fall 2022  
> **Project Type**: Group Final Project  

---

## Overview
This project explores genetic sequence variants from 57 patients diagnosed with Acute Myeloid Leukemia (AML). The goal was to identify protein-coding genes that may be associated with AML by analyzing tumor and normal samples using bioinformatics and statistical methods.

We processed and annotated VCF files, conducted comparative mutation analysis, and applied tools like FATHMM to identify potentially pathogenic variants.

---

## Objectives
- Parse and process genetic variant data from VCF files
- Construct a binary matrix of SNV occurrences across samples
- Analyze mutation frequency patterns
- Annotate SNVs to identify coding and non-coding regions
- Identify nonsynonymous SNVs (nsSNVs)
- Use FATHMM to score potential pathogenicity
- Score genes for potential association with AML

---

## Key Findings
- Out of 7238 SNVs in coding sequences, 4117 were nonsynonymous
- Tumor samples showed a higher frequency of mutations than normal samples
- FATHMM scores helped filter out likely pathogenic nsSNVs
- Several genes showed recurrent mutation patterns and high impact scores

---

## Visualizations

### Top Scoring AML-Associated Genes
<p align="center">
  <img src="./top_genes_by_score.png" alt="Top 10 Genes by Score" width="400"/>
</p>

---

## How to Run the Analysis
1. **Environment Setup**
   - Python 3
   - Required packages: `pandas`, `PyVCF`

2. **Prepare VCF Files**
   - Place VCFs into the appropriate input folder (see paths inside scripts)

3. **Run Scripts in Order**
   ```bash
   python VCF2CSV.py                  # Converts raw VCFs to CSV
   python parsevcf.py                 # Adds gene annotations and change type
   python FATHHM.py                   # Adds FATHMM scores
   python SNPnexus.py                 # Adds SNPnexus scores
   python tumor_vs_normal.py          # Merges data and analyzes variant stats
   
   ```

4. **Review Outputs**
   - Outputs are saved as CSV or Excel files in the `results/` folder

---

## Tools Used
- Python (pandas, os, vcf)
- FATHMM (Functional Analysis through Hidden Markov Models)
- SNPnexus
- Reference genome annotations from UCSC

---

## Notes
- This project was completed in Fall 2022 as part of coursework in BINF/STAT 5354.
- Some file paths and directories inside the scripts are specific to local environments and may need to be updated to run correctly.

---

## References
- [Cancer.org - AML Overview](https://www.cancer.org/cancer/acute-myeloid-leukemia/about/what-is-aml.html)
- [Mayo Clinic - AML](https://www.mayoclinic.org/diseases-conditions/acute-myelogenous-leukemia/symptoms-causes/syc-20369109)
- [OMIM - AML Genes](https://www.omim.org)

---

## Contact
For questions or collaboration, feel free to open an issue or fork the repo.
