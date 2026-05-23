# SKO-1

### Contributors

* **Elena Kokinos**: management, BCR analysis
* **Olga Loginova**: IFN status classification
* **Alexey Senkovenko**: nextflow pipeline, HLA analysis
* **Elizaveta Perepelitsa**: Gene expression analysis 
* **Sergei Losev**: TCR analysis
* **Anastasia Patrusheva**: HLA analysis

Supervisor: **Pavel Skoptsov**

## Introduction

Multiple sclerosis (MS) and systemic lupus erythematosus (SLE) are chronic autoimmune diseases that differ in their clinical features and underlying mechanisms, but both involve serious immune system dysfunction. Type I interferons (IFNs), particularly interferon-beta (IFN-β), play opposite roles in these conditions. In MS, IFN-β is commonly used as a first-line immunomodulatory treatment that decreases disease activity and reduces the frequency of relapses [[Jakimovski]](https://doi.org/10.1101/cshperspect.a032003). In SLE persistent overactivation of the type I IFN pathway known as the "interferon signature" or interferonopathy is a major factor driving inflammation, the breakdown of immune tolerance, and the production of autoantibodies [[Ishihara]](https://doi.org/10.3390/biom15111586).
Interestingly, both MS and SLE can involve the production of neutralizing autoantibodies against IFN-β. In MS, these antibodies often appear as a result of IFN-β therapy and may reduce how well the treatment works [[Gilli]](https://doi.org/10.1093/brain/awh028). In SLE, they may arise spontaneously as part of the natural autoimmune response [[Grenmyr]](https://doi.org/10.1177/09612033261432154). Still, only a subset of patients in each disease develops these antibodies. Understanding why this heterogeneity exists could help identify shared molecular mechanisms that explain differences in response to interferon and the loss of immune tolerance.


## Nextflow pipeline

### Installation

Check Java installation (must be above Java 11).

```
java -version
```

If Java is not installed.
```
curl -s https://get.sdkman.io | bash
sdk install java 17.0.10-tem
java -version 
```
Install Nextflow in separate conda environment using [official guide](https://docs.seqera.io/nextflow/install#conda)

The pipeline was adapted to the latest Nextflow version [v. 26.04.0](https://doi.org/10.1038/nbt.3820).

### Tree of pipeline:
	
```
├── bin
│   └── download_fastq_from_csv.sh
├── envs
│   ├── arcashla.yaml
│   ├── download_aria.yaml
│   ├── fastp.yaml
│   ├── kallisto.yaml
│   ├── mixcr.yaml
│   ├── star.yaml
│   └── trust4.yaml
├── main.nf
├── modules
│   ├── MIXCR.nf
│   ├── TRUST4.nf
│   ├── arcasHLA.nf
│   ├── download_aria.nf
│   ├── fastp.nf
│   ├── kallisto.nf
│   └── star.nf
└── nextflow.config
```

* `main.nf` – the core file containing workflow and orchestrating the use of modules.
* `nextflow.config` – the core file with settings for working profiles and executors, contains parameters for running different modules.
* `modules` – separate nf-files with process logic.
* `envs` – .yaml files to build separate conda environments for each tool in workflow.
* `bin` – directory to contain other executable scripts.

Processes include:

* `DOWNLOAD_ARIA` – process to download fastq files by SRR IDs. Uses [aria2](https://aria2.github.io/) as high-speed multiprotocol utility for downloading files.
* `FASTP` – process using [fastp](https://github.com/opengene/fastp) for preprocessing and QC.
* `STAR` – process using [STAR](https://github.com/alexdobin/STAR) splice-aware alignment tool to generate .bam files.
* `MIXCR` – process running [MIXCR](), commercial standard tool for NGS-data analysis. In this project MIXCR is used for clonotype repertoire reconstruction. If a MIXCR license is available, it is recommended to use this module.
* `ARCAS_HLA` – process for HLA typing, uses [arcasHLA](https://github.com/RabadanLab/arcasHLA) tool.

Additional processes:
* `KALLISTO` – process to replace STAR, to save computing resources and retain high-speed. Uses [kallisto](https://github.com/pachterlab/kallisto) pseudoalignment tool.
* `TRUST4` – process that uses [TRUST4](https://github.com/liulab-dfci/TRUST4) for clonotype repertoire reconstruction, can be used as an alternative for MIXCR when license is unavailable.

### Prepare to work


### Future adjustments
* Make separate docker profile and docker images to operate with.

## Data Collection & Cohort


Transcriptomic data for **MS** and **SLE** were retrieved from public repositories (GEO/ArrayExpress). From an initial pool of 2,300+ samples, we curated a high-quality cohort based on strict technical and clinical criteria.


**Inclusion Criteria:**
*  **Sequencing:** Paired-end reads, length ≥100 nt.
*  **Source:** Peripheral blood (to capture diverse immune populations).
*  **Clinical:** Baseline samples only (treatment-naïve) to minimize confounding.


**Final Cohort (n = 433):**


| Group | Count |
|-------|-------|
| **SLE** | 249 |
| **CLE** | 62 |
| **MS**  | 63 |
| **Healthy Controls** | 59 |

![](images/data_cohorts.png)

## Classifier

Supervised Random Forest classifier for stratifying samples into **IFN-High** / **IFN-Low** categories based on the expression of 9 canonical type I interferon response genes (*IFI27, IFIT1, IFIT3, MX1, OAS1, RSAD2, STAT1, IFNAR1, IFNAR2*).

**Key Details:**
*   **Training Data:** SLE cohort (GSE116006, n=152) with RT-PCR validated interferon status.
*   **Implementation:** Python `scikit-learn` (v.1.6.1), stratified 5-fold CV, optimized threshold (0.67).
*   **Performance:** AUC = 0.996, Test Accuracy = 96.8%.

>  *Validated on SLE; application to other conditions (e.g., MS) requires further validation.*

![](images/classifier.png)

## Gene expression analysis

A key focus of the analysis was evaluating the biological validity of our classification by assessing an Interferon type I (IFN) score. Using a curated set of IFN-inducible genes (IFI27, IFIT1, IFIT3, MX1, OAS1, RSAD2, STAT1, IFNAR1, IFNAR2), we calculated an average expression score, which was Z-scored separately within each dataset to eliminate batch effect. 
Mann—Whitney U test was applied to compare the IFN scores between the samples predicted as "IFN-High" and "IFN-Low" by our classifier, grouped by diagnosis, as well as by potential technical (RNA type, object) and biological (sex) confounders. The significant difference observed in these scores between "IFN-High" and "IFN-Low" samples supports the biological relevance of the classification. 

![](images/ifn_score.png)

Additional gene expression analysis involved calculating similar scores for various immune, matrix, and angiogenesis gene signatures to confirm the integrity of the blood-derived data.

To explore the transcriptomic differences between specific groups, we performed differential gene expression analysis using [DESeq2](https://github.com/thelovelab/DESeq2). In particular, we assessed the differences between healthy controls, as the presence of healthy samples classified as "IFN-High" by our model was unexpected. Notably, these "IFN-High" healthy samples were found in different data sets, eliminating the possibility of a confounding batch effect. Interestingly, only IFN-inducible antiviral genes were differentially expressed between "IFN-High" and "IFN-Low" healthy controls.

Results are summarized in the [notebook](notebooks/expression_analysis.ipynb).

## TCR
T-cell receptor repertoire analysis was based on [MiXCR](https://github.com/milaboratory/mixcr/) clonotype tables generated from .fastq files. Frequency matrices were generated from those tables and we used perMANOVA and pairwise Mann-Whitney criterion with FDR in order to compare V, D and J-segment gene usage profiles between high and low interferon response groups. 

For V-segments of TCR beta chains only 2 genes have shown higher frequencies in SLE samples with Low interferon beta response. This might imply their involvement in anti-interferon beta autoimmune response in SLE, but it should be noted that the effect size is small so these 2 changes are to be regarded as potential candidates with caution. No significant differences have been observed for D- and J-gene usage in the TCR beta-chain.

You can find all the results of this analysis in these 2 notebooks:
[TCR_analysis_diversity.ipynb](notebooks/TCR_analysis_diversity.ipynb)
[TCR_analysis_gene_usage.ipynb](notebooks/TCR_analysis_gene_usage.ipynb)

MiXCR clonotype tables for both TCR and BCR are available at this [link](https://drive.google.com/file/d/1F0lG7-75RYOUBOmnxsbJyyBxJdB44i3N/view?usp=drive_link).

## BCR

BCR repertoire analysis was performed using MiXCR, followed by downstream processing of clonotype tables in a Jupyter notebook. The analysis focused on heavy-chain BCR V(D)J segments and their usage frequencies across samples.

In particular, the frequencies of IGHV, IGHD, and IGHJ segments were compared between groups with high and low predicted interferon status in patients with multiple sclerosis (MS), systemic lupus erythematosus (SLE), and cutaneous lupus erythematosus (CLE). This analysis was designed to identify repertoire features potentially associated with reduced IFN status, which was hypothesized to reflect the presence of neutralizing antibodies against type I interferons.

For each sample, V, D, and J segment usage frequencies were derived from MiXCR clonotype tables. Segment usage was summarized at the diagnostic-group level using median within-sample frequencies, and the top 10 segments were identified separately for each segment class (IGHV, IGHD, IGHJ). These top segments were then compared between Low- and High-IFN groups. Segment usage frequencies were compared between groups using a two-sided Mann-Whitney U test with Benjamini-Hochberg FDR correction for multiple testing.

Significant differences in several BCR segments were detected between Low- and High-IFN groups in SLE, whereas no such differences were observed in the other diagnostic groups.
We next examined whether these effects remained consistent across individual GSE datasets or whether they were primarily driven by dataset-specific signals.
For each selected segment, IFN-Low and IFN-High samples were compared separately within each GSE using a two-sided Mann-Whitney U test. This produced a dataset-specific estimate of effect size and direction for every segment.

![](images/IGHD3-10.png)

Among the analyzed BCR segments, only IGHD3-10 demonstrated a significant and reproducible difference between IFN-High and IFN-Low groups in SLE across independent datasets. No analogous segment-level effects were observed in MS or CLE.

Results of the analysis are contained in [BCR_analysis](notebooks/BCR_analysis.ipynb) notebook.

## HLA typing

To analyze the results of HLA genotyping obtained using arcasHLA, the data were merged into a tsv-table. In the  [notebook](notebooks/hla_analysis.ipynb) for HLA-analysis, there this tsv-table is concatenated with the metadata.

To analyze the results of HLA genotyping, data were divided into ethnic groups, since different ethnicities are characterized by a different set of HLA alleles. Within each ethnic group, the observations were grouped by disease.

All statistical processing methods are combined in the HLA class:

`count alleles` - calculates the occurrence of alleles for each group of *diagnosis-IFN status* for the uploaded batch; the output is a pandas.Dataframe (which is saved as an attribute `allele_counts`) containing the number of carriers of this allele, the number of copies of it, the number of its occurrences in homozygotes and heterozygotes, the number of occurrences of this allele in the *diagnosis-IFN status* group (`n_samples`), as well as the frequency of occurrence in fractions by carriers and the number of copies;
based on `allele_counts`, the Fisher’s exact test and Chi-Square tests can be applied, at the same time, if allele has low expected counts (<5) or zero expected frequencies when forming the contingency table, Fisher’s exact test instead of Chi-Square is applied;
barplots and pie charts are also drawn based on `allele_counts`.

Forest plot is used for meta-analysis, and a separate module is imported to create it.

![](images/HLA_forest_plot.png)

Specific usage examples are given in [notebook](notebooks/hla_analysis.ipynb) for HLA-analysis.


