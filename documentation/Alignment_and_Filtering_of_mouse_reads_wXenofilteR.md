# Alignment and Filtering of mouse reads from Human BAM file from Xenografted samples

## Overview

This document describes the steps taken to align and filter out mouse reads from Human BAM files from the Patient Derived Xenograft (PDX) samples. Mapping parameters used to align reads against the Human GRCh38 reference genome and the mouse reference used [NOD_ShiLtJ_V1_PDX](../reference/NOD_ShiLtJ_V1_PDX_ref/README.md) were the same.

All the scripts and code mentioned below can be found in the `scripts` directory.

## Alignment to the Human GRCh38 reference genome 

The WES sequencing data was aligned to the GRCh38 Human reference genome using `bwa mem`. PCR duplicates were marked using `samtools markdup` function. The same process was applied for all the samples from both targeting experiments in this project. This process was preformed through an internal pipeline.

## Filtering of mouse reads with Xenofilter

This process was used to remove mouse reads from the Human BAM files with the whole exome sequencing data generated for the PDX samples.

### Required Environment variables and software

The following environment variables are required to be set before running the scripts:
- **PROJECTDIR**: The path to the project directory where this repo got cloned into
- **STUDY**: The study ID,  6633 for this analysis
- **PROJECTID**: The project ID, 2729 for this analysis

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.2.2`
- **samtools**: samtools `v1.14`
- **bwa mem**: bwa mem `v0.7.17`
- **XenofilteR**: XenofilteR `v1.6`

- Load the following variables and software
```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729

# Example of a source me file containing the required environment variables and software to load
# iRODs is used internally to access the sequencing data but this should be replaced for the human unfiltered BAM files 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh
```
#### R environment
If you're interested in reproducing the R environment, for a the code used in R 4.2.2 run the following commands within `R v4.2.2`, change the path on `projectdir` to the path where the repository was cloned into:

```R
projectdir<-"/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES"
pdx_processing_dir<- file.path(projectdir,"scripts/pdx_processing")

setwd(pdx_processing_dir)
install.packages("renv")
library(renv)
renv::restore(lockfile=file.path(pdx_processing_dir, "renv.lock")) # To rebuild an environment from the renv.lockfile

```

### CRAM to FASTQ files

#### Manifest generation and import of sequencing data

To import the sequencing metadata from iRODS to generate manifests with the sequencing statistics and information
we ran the following: 

**IMPORTANT**: All the manifest generated can be found within the `metadata/manifests` directory.

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing


#load iRODS module 
module load IRODS/1.0
# initiate session 
iinit

#To call the script to build the manifests with project ID Name
mkdir -p ${PROJECTDIR:?unset}/metadata/manifests

source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/Build_manifest_from_irods_cram_information.R --seqscape_proj_id ${STUDY} --outdir ${PROJECTDIR:?unset}/metadata/manifests
```

**OUTPUTS**:
- **6633_cram_manifest_INFO_from_iRODS.txt** : Contains the information of the samples.

- After generating the manifest, we split the information to only contain the Patient Derived Xenografted (PDX) samples

```bash 
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

# Reformat the manifest and filter 
mv ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS.txt ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all.txt
# Filter only the information for the PDX samples but keep the header of the file
head -n 1 ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all.txt >${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_PDXs.txt
grep -f ${PROJECTDIR:?unset}/metadata/${STUDY}_${PROJECTID}_unfilt_PDX_sample_names.tsv ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all.txt >>${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_PDXs.txt

```
**OUTPUTS**:
- **6633_2729_cram_manifest_INFO_from_iRODS_all.txt** : Contains the information of all the samples for both tumours and PDX samples.
- **6633_2729_cram_manifest_INFO_from_iRODS_PDXs.txt** : Contains the information of the PDX samples only.

- To generate the list of jobs to transform the CRAM files to fastq files the files we ran

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

cd ${PDXSCRIPTS_DIR:?unset}

# Load environment with requiring 
source ${PDXSCRIPTS_DIR:?unset}/source_me.sh

#This script takes the cram manifest and generates the SH file with the jobs to import and transform to fastq all of the cram files from iRODs
Rscript ${PDXSCRIPTS_DIR:?unset}/cramtofastq_from_iRODs_based_cram_manifest.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS_PDXs.txt --projectdir ${PROJECTDIR:?unset} --studyID ${STUDY} --mem 16000

```
Output: 
- `scripts/pdx_processing/6633_cramtofastq_from_iRODs_jobs.sh` : Contains the list of jobs to transform the CRAM to fastq files

- Then we proceed to execute the jobs import the BAM files and transform them into fastqs using samtools

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

cd ${PROJECTDIR:?unset}
# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh
# Login into iRODs
#iinit 
sh ${PDXSCRIPTS_DIR:?unset}/${STUDY}_cramtofastq_from_iRODs_jobs.sh

```

#### Generate the mouse genome reference files and bwa index

To be able to generate the mouse referenced use the steps mentioned in the [**NOD_ShiLtJ_V1_PDX_ref**](../reference/NOD_ShiLtJ_V1_PDX_ref/README.md) README file. 


#### Mapping against the NOD_ShiLtJ_V1_PDX mouse reference genome with bwa-mem 

To generate the jobs to map the fastq files against the mouse reference genome, we used the `PDX_bwa_mem_mapping_jobs_from_master_manif.R` script:

- **INPUT**: Use the file : `metadata/manifests/6633_cram_manifest_INFO_from_iRODS_wbam_counts_qc.txt`

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

NOD_PDXV1_REFDIR=${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref

# Load environment with requiring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

# SET the directory with the UNFILTERED BAM files by symlinking bam and bai files to the data bams unfiltered directory
WES_UNFILT_BAMDIR=${PROJECTDIR:?unset}/data/bams/WES_UNFILT
STAGEDIR=/lustre/scratch125/casm/staging/team113/${PROJECTID}

mkdir -p ${PROJECTDIR:?unset}/data/bams/WES_UNFILT
for sample in $( cat ${PROJECTDIR:?unset}/metadata/${STUDY}_${PROJECTID}_unfilt_PDX_sample_names.tsv); do
	mkdir -p ${WES_UNFILT_BAMDIR:?unset}/${sample}
	ln -s ${STAGEDIR}/${sample}/mapped_sample/${sample}.sample.dupmarked.bam ${WES_UNFILT_BAMDIR:?unset}/${sample}/
	ln -s ${STAGEDIR}/${sample}/mapped_sample/${sample}.sample.dupmarked.bai ${WES_UNFILT_BAMDIR:?unset}/${sample}/
	ln -s ${STAGEDIR}/${sample}/mapped_sample/${sample}.sample.dupmarked.bas ${WES_UNFILT_BAMDIR:?unset}/${sample}/
done

# Add column to the tab delimited file at the end entitled Proc_as_PDX to the manifest with the PDX samples and Y on every row,
#Add PDX to the header as a last column in the manifest ${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc.txt  keepeing UTF-8 encoding
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{$(NF+1)="Proc_as_PDX"} NR>1{$(NF+1)="Y"} 1' ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc.txt > ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot.txt 


#This scrip
Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/PDX_bwa_mem_mapping_jobs_from_master_manif.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot.txt --projectdir ${PROJECTDIR:?unset} --referencedir ${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref/bwa_mem

```
This will generate two outputs:
 1. `bwamem_mapping_perlanrun_to_NOD_PDXV1_tum_only_jobs.sh`: which contains with the list of jobs to perform the mapping against the mouse reference genome with **bwa-mem**.
 2. `samtools_psample_merge_nodv1_tum_only_jobs.sh`: which contains the list of jobs to merge and index the BAM files per sample.

Submit the remapping jobs with the mouse reference using **bwa-mem**

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

NOD_PDXV1_REFDIR=${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref

# Load environment with requiring 
cd ${PROJECTDIR:?unset}/scripts/pdx_processing

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

sh bwamem_mapping_perlanrun_to_NOD_PDXV1_tum_only_jobs.sh

```

Submit the merging per sample

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing

# Load environment with requiring 
cd ${PDXSCRIPTS_DIR:?unset}

# Load environment with requring 
source ${PDXSCRIPTS_DIR:?unset}/source_me.sh

sh samtools_psample_merge_nodv1_tum_only_jobs.sh
```

### Filter the mouse reads from the Xenofilter jobs for the samples

First, we generated the manifest with the final filtered file name and example of the Xenofilter jobs for each sample if a single per samples file could have been used using the `run_Xenofilter_from_WES_master_manif.R` script.  

However, given the amount of sequencing data generated for the grafted samples, we required to split both Human and mouse BAM files with the same matching read numbers and read names to be able to run XenofilteR successfully.  This is due to a known issue with XenofilteR when and one of it's dependencies, `Rsamtools`, see issue [here](https://github.com/NKI-GCF/XenofilteR/issues/7) 


We generate the manifest:

**INPUT**: `metadata/manifests/6633_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot_psamp_mouse.txt`

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing
NOD_PDXV1_REFDIR=${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref

# Load environment with requiring 
cd ${PROJECTDIR:?unset}/scripts/pdx_processing

# Load environment with requring 
source ${PROJECTDIR}/scripts/pdx_processing/source_me.sh

# Set the job to create the .sh XenofilteR job submissions 
Rscript ${PROJECTDIR}/scripts/pdx_processing/run_Xenofilter_from_WES_master_manif.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot_psamp_mouse.txt --projectdir ${PROJECTDIR} --outdir ${PROJECTDIR}/data/bams/WES_xfilt
```
**OUTPUTS**:
 `metadata/manifests/6633_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot_psamp_mouse_xfb.txt`


#### Split input BAM files, split by read names from Human BAM files and run XenofilteR for mouse read filtering

To split the reads we used the script: **split_bam_files_and_get_xenofilter_jobs.R**. 

The script takes the a manifest with BAM file information and generates the jobs to take the unfiltered human BAM files, sort by read name, obtain a plain text file with the read names for all the reads present in the file, then it will split this file in the number of files required that have a maximum of `NREADS_SPLIT` per file. **This approach was followed as the samples were sequenced across a sinlge land and only a single readgroup was present.** Subsequently it will split the Human and Mouse BAM files by read names and then generate the jobs to run XenofilteR to filter out the mouse reads for the matching Human & Mouse BAM files.

**INPUT**: `metadata/manifests/6633_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot_psamp_mouse_xfb.txt`

```bash
PROJECTDIR=/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES
STUDY=6633
PROJECTID=2729
SCRIPTS_DIR=${PROJECTDIR:?unset}/scripts
LOGDIR=${PROJECTDIR:?unset}/logs
PDXSCRIPTS_DIR=${SCRIPTS_DIR:?unset}/pdx_processing
NOD_PDXV1_REFDIR=${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref

cd ${PROJECTDIR:?unset}/tmp


# Load environment with requiring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

# Variables:
NCORES=8
NSPLIT_CORES=22
NREADS_SPLIT=50000000

bsub -q basement -n 24 -M64000 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]" -J ${STUDY}_split_bam_files_and_get_xenofilter_jobs_WES -o ${LOGDIR:?unset}/${STUDY}_split_bam_files_and_get_xenofilter_jobs_WES_b2_b3.o -e ${LOGDIR:?unset}/${STUDY}_split_bam_files_and_get_xenofilter_jobs_WES_b2_b3.e \
"Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/split_bam_files_and_get_xenofilter_jobs.R --study_id ${STUDY:?unset} \
--manifest ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_PDXs_wbam_counts_qc_PDX_annot_psamp_mouse_xfb_b2_3.txt \
--projectdir ${PROJECTDIR:?unset} \
--xfilter_outdir ${PROJECTDIR:?unset}/bams/WES_xfilt/NOD_PDXV1 \
--split_nreads ${NREADS_SPLIT:?unset} \
--nsplit_cores ${NSPLIT_CORES:?unset} \
--nthreads ${NCORES} 
"

/lustre/scratch125/casm/teams/team113/projects/6633_2729_3248_PDX_models_from_Latin_America_WES/metadata/6633_2729_unfilt_PDX_sample_names_batch2_3.txt

```

**OUTPUTS**:
The scripts will take the unfiltered human BAM files, sort by read name, split by read name  following output files:
- [`6633_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt`](../metadata/manifests/6633_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt): Manifest that contains the metadata and file name information of the split files for all the samples.

Runner files:
- **6633_bamsplittin_by_read_names_jobs.sh** : Contains the list of jobs to split the BAM files by read names using `samtools` 
- **6633_Xenofilter_by_read_names_parts_jobs.sh** : Contains the list of jobs used to run XenofilteR for the split Human and Mouse BAM files using the script [`scripts/pdx_processing/bam_Xenofilter_rg.R`](../scripts/pdx_processing/bam_Xenofilter_rg.R)
- **6633_Xenofilter_merged_filtered_parts_jobs.sh**: Contains the list of jobs to merge the filtered BAM files per sample and index them using `samtools`

To submit the BAM file splitting by read names:

```bash
PROJECTDIR=/lustre/6633_PDX_models_Latin_America_WES
STUDY=6633

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_bamsplittin_by_read_names_jobs.sh

```

Submit Xenofilter filtering of split bamfiles using `6633_Xenofilter_by_read_names_parts_jobs.sh`

```bash
PROJECTDIR=/lustre/6633_PDX_models_Latin_America_WES
STUDY=6633

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_Xenofilter_by_read_names_parts_jobs.sh
```

Submit filtering of BAM files with `XenofilteR`

```bash
PROJECTDIR=/lustre/6633_PDX_models_Latin_America_WES
STUDY=6633

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_Xenofilter_merged_filtered_parts_jobs.sh
```

### Collate filtered read stats and plots

Finally, once all jobs are complete we collate information on the proportion of filtered reads per file per sample and plot it using the script `collate_xfilter_stats_from_manif_and_plots.R `

INPUT: `metadata/manifests/6633_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt`

```bash
PROJECTDIR=/lustre/6633_PDX_models_Latin_America_WES
STUDY=6633

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

# Run the collation of the xenofiltered stats information and plots
Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/collate_xfilter_stats_from_manif_and_plots.R --study_id ${STUDY:?unset} \
--manifest ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt \
--projectdir ${PROJECTDIR:?unset} \
--xfilter_outdir ${PROJECTDIR:?unset}/analysis/bam_xfilterstats 
```

**OUTPUTS**:
- [`6633_cohort_per_part.xfiltstats`](../analysis/bam_xfilterstats/6633_cohort_per_part.xfiltstats): Contains the information of the proportion of filtered reads per file per sample.
- [`6633_cohort_per_sample.xfiltstats`](../analysis/bam_xfilterstats/6633_cohort_per_sample.xfiltstats): Contains the information of total of filtered reads per sample.
- [`6633_per_sample_pdx_xfiltstats_bpl.pdf`](../analysis/bam_xfilterstats/6633_per_sample_pdx_xfiltstats_bpl.pdf): Bar plot with the percentage of reads filtered per sample.
- [`6633_per_sample_pdx_xfiltstats_vpl.pdf`](../analysis/bam_xfilterstats/6633_per_sample_pdx_xfiltstats_vpl.pdf): Bar plot with the percentage of reads filtered per sample.

