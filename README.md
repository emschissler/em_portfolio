# Em Schissler - Bioinformatics Portfolio

Welcome to my bioinformatics portfolio! This repository contains various projects showcasing my expertise in data analysis, bioinformatics workflows, and RNA-seq analysis. Below, you will find descriptions and files for each project, detailing the objectives, methodologies, and key outcomes.

---

## 1. **Differential Gene Expression (DGE) Analysis Workflow (Spring 2024)**

This repository contains a complete RNA-seq pipeline for identifying differentially expressed genes (DEGs) between conditions. Key workflow steps include pre-processing RNA-seq data, conducting statistical analysis using DESeq2, and generating visualizations to interpret DEG results.

### Key Steps:
- **Data Pre-processing:** Trimming adapters and pseudo-aligning reads.
- **Statistical Analysis:** Conducting differential expression analysis using DESeq2.
- **Data Visualization:** Generating PCA plots and MA plots to visualize gene expression changes.

**Files Include:** 
- `RNA-seq Analysis.pdf`
- `RNA-seq Analysis.md`

---

## 2. **Gene Read Matching: Parsing SAM/BAM and GTF Files (Fall 2023)**

This repository contains a bioinformatics workflow for parsing SAM/BAM files and comparing them with GTF gene annotation files to determine the number of reads aligned to specific genes. The workflow involves:
- Reading and parsing SAM/BAM files (containing aligned sequencing reads).
- Comparing aligned reads to gene annotations (from GTF files), sorted by chromosome, exon start, and end positions.
- Counting how many reads align to each gene based on the gene's position.

This project demonstrates methods for efficient parsing and comparison of large genomic datasets using example data.

**Files Include:** `Gene Read Matching.ipynb`

---

## 3. **Analysis of Biomarkers in Acute Kidney Injury within COVID Patients (Summer 2023)**

This project focuses on identifying potential biomarkers that may predict acute kidney injury in COVID-19 patients. Using clinical and imaging data from The Clinical Proteomic Tumor Analysis Consortium (CPTAC) and The Cancer Imaging Archive (TCIA), the analysis explores relationships between biomarker levels and clinical outcomes.

### Objectives:
1. **Identify Predictive Biomarkers:** Explore whether clinical variables like potassium and creatinine can predict acute kidney injury.
2. **Examine Biomarker Patterns:** Investigate biomarker level patterns to predict kidney injury in hospitalized patients.

**Files Include:** 
- `Analysis of Biomarkers in Acute Kidney Injury.html`
- `RAWCOVID.csv`
- `COVIDMETA.csv`

---

## 4. **Analysis of Factors Affecting Mammal Longevity (Spring 2021)**

This project investigates how body mass and diet influence mammalian lifespan. Traditional assumptions about the relationship between these factors and longevity are challenged, and this analysis provides insights into mammalian lifespans using datasets from the Animal Diversity Web (ADW) and MammalDIET.

### Objectives:
1. **Analyze the Impact of Diet and Body Mass:** Evaluate whether diet and body mass significantly affect mammalian lifespan.
2. **Validate Data Representativeness:** Compare ADW and MammalDIET datasets to ensure accurate representation of mammalian diets.

**Files Include:** 
- `Analysis of Factors Affecting Mammal Longevity.html`
- `FinalMammalDIET.csv`
- `MammalMetadata.xlsx`

---

Feel free to explore each repository and view the projects in detail. If you have any questions or would like to discuss the projects further, don't hesitate to reach out!

