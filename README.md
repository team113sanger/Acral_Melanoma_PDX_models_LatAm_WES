# Alignment, filtering and somatic variant calling of PDX from Latin - Whole Exome Sequencing
[![DOI]()]()

## Overview 
This repository contains the code for the alignment, filtering of mouse reads and variant calling performed in relation to the results presented in the manuscript:

_Modelling Acral Melanoma in Admixed Brazilians Uncovers Genomic Drivers and Targetable Pathways_ paper.

## Experimental design

The table containing the metadata of all the samples and the experimental group and conditions that each sample belongs to can be found in the (`Group` column) on the [metadata table here](./metadata/7688_3365_sample_exp_metadata.tsv) inside the `metadata` directory. 

## Analysis workflow overview
The analyses were performed using a combination of R, bash and Perl scripts. The workflow depicted below shows the steps followed in the analysis of the sequencing data and the differences required due to the experimental designs.
![image](./documentation/diagrams/WES_analysis_final_workflow.drawio.png "Analysis workflow")

The overall analysis workflow after the sequencing data was generated was divided into two stages:

2. **Data Pre-processing**: This section includes the alignment to the required reference genome GRCh38 for Human and [NOD_ShiLtJ_V1_PDX](./reference/NOD_ShiLtJ_V1_PDX_ref/README.md) customised mouse reference, and filter of mouse reads from the Human mapped BAM files using **XenofilteR**. Alignment against the Human reference genome was performed in the same way for both datasets. Instructions on how the custom mouse reference and analyses at these stage were performed can be found in the [README here](./documentation/Alignment_and_Filtering_of_mouse_reads_wXenofilteR.md). 

3. **Coverage depth QC analysis**: This section includes the scripts used for the analysis of the coverage depths across the pulled regions and the filtering of samples based on coverage depth. See [README](./documentation/Coverage_depth_check.md) file in the `documentation` folder for more details on how to reproduce the analysis.

3. **Variant calling and annotation**: This section includes the SNV, MNV and INDEL variant calling, filtering, and variant effect prediction of the variants using negative control lines as normal samples. Resulting variants were collated and VCFs were transformed into a Mutation Annotation Format(MAF) file. Instructions on how the analyses at these stage were performed can be found in the [README here](./documentation/Off-pipe_calling_of_Xenofiltered_WES_data.md)

4. **Variant calling and annotation**: This section includes the scripts used for variant calling from whole exome sequencing data. It includes the generation of VCF files, the filtering of variants based on various criteria and finally the generation of summary MAF files and plots. See [README](./documentation/Somatic_Variant_calling.md) file in the `documentation` folder for more details on how to reproduce the analysis.

5. **Variant calling and annotation**: This section includes the SNV, MNV and INDEL variant calling, filtering, and variant effect prediction of the variants using negative control lines as normal samples. Resulting variants were collated and VCFs were transformed into a Mutation Annotation Format(MAF) file. Instructions on how the analyses at these stage were performed can be found in the [README here](./documentation/Off-pipe_calling_of_Xenofiltered_WES_data.md)


## Results 
 
All the variants that passed standard filtering criteria and had a Variant Allele Frequency >0.04 were retained. The results  and plots used in the paper are located within the `analysis/somatic_variant_plots` directory. 

To reproduce the results files and plots, please follow the instructions in the [README here](./documentation/Somatic_variant_plotting.md).

## Software dependencies
Analyses were performed using a combination of R, Perl and shell scripts.

- All R scripts were run using `R v4.2.2` and the packaged dependencies for each analysis are detailed within the `renv.lock` file of each analysis folder in the `scripts` directory.
- Perl scripts were run using Perl version `v5.38.0`
- The following software was used:
  - `samtools` version`v1.14` [**here**](https://github.com/samtools/samtools)
  - `bwa-mem` version `0.7.17` [**here**](https://github.com/lh3/bwa)
  - `XenofilteR` version `1.6` [**here**](https://github.com/NKI-GCF/XenofilteR/releases/tag/v1.6)
  - `bcftools` version `1.9` [**here**](https://github.com/samtools/bcftools/)
  - `tabix` version `1.9` [**here**](https://github.com/samtools/tabix/)
  - `CaVEMan` version `1.18.2` [**here**](https://github.com/cancerit/CaVEMan)
  - `cgpCaVEManwrapper` version `1.18.2` [**here**](https://github.com/cancerit/cgpCaVEManWrapper)
  - `Smart-Phase` version `0.1.8`[casm docker image here](https://github.com/cancerit/CASM-Smart-Phase/tree/main)
  - `cgpPindel` version `3.11.0` [**here**](https://github.com/cancerit/cgpPindel)
  - ENSEMBL VEP version `103`[**here**](http://feb2021.archive.ensembl.org/info/docs/tools/vep/index.html)
- The repositories with the scripts used for variant QC and VCF to MAF conversion can be found in the following links:
    - [**QC**](https://github.com/team113sanger/dermatlas_analysis_qc) `v0.5.2`
    - [**MAF**](https://github.com/team113sanger/dermatlas_analysis_maf) `v0.6.4` 


## Contact

- Martin Del Castillo Velasco-Herrera (mdc1@sanger.ac.uk)
- Jacqueline Boccacino (jb62@sanger.ac.uk)