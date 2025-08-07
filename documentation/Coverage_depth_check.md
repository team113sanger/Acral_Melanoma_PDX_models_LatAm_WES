# Coverage depth check for somatic variant calling	

## Overview

To analyse the evenness of coverage depth of the samples, we used a file that contained the coordinates of the baits used for exome capture. Coordinates from hg19 were converted to GRCh38 coordinates. Following this, coverage observed was performed using `samtools depth` command to calculate the coverage depth of each sample. The output was then processed to generate a summary of the coverage depth across all samples.

Samples with a coverage depth of **>20X across at least 80% of the baits** used were considered for somatic variant calling. 

## Dependencies

## Required software

The scripts used for this analysis are located in the `scripts/coverage_qc` directory. The analysis was performed using a combination of R, Perl, and shell scripts. The following software versions were used:

- Perl scripts were run using Perl version `v5.38.0`[**here**](https://www.perl.org/)
- `samtools` version `1.14` [**here**](https://www.htslib.org/)
- `R` version `4.3.3` [**here**](https://www.r-project.org/)
- `bwa-mem` version `0.7.17` [**here**](https://github.com/lh3/bwa)
- `biobambam2` version `2.0.146` [**here**](https://github.com/gt1/biobambam2)

### Required datasets

- The bait set used for the exome capture was the Agilent SureSelect Human All Exon V6+UTR, the file was obtained from Agilent (file S07604624_hs_hg19.zip). Coordinates for GRCh38 were obtained by performing a liftover from hg19 to GRCh38 coordinates using [`LifOver`](https://genome.ucsc.edu/cgi-bin/hgLiftOver). See file [`resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed`](resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed) 

## Data alignment and processing

The WES sequencing data was aligned to the GRCh38 Human reference genome using `bwa-mem` version `0.7.17` [**here**](https://github.com/lh3/bwa). PCR duplicates were marked using `bammarkduplicates2` from biobambam2 version `2.0.146`[**here**](https://github.com/gt1/biobambam2). The same process was applied for all samples and done through an internal pipeline.


## Coverage Depth analysis

The Coverage Depth analysis was run under a cluster environment with **Ubuntu 22.04** and **LSF job execution**. It is required to replace the `PROJECTDIR` shell environment variable with the path of where the current repository was downloaded. 

### Download BAM files 

BAM files are available for Download from EGA repository with Study Accession ID **EGAS00001008212**, EGA website [*here**](https://ega-archive.org/).  The BAM files should be named as `<sample_name>.sample.dupmarked.bam`, where `<sample_name>` is the name of the sample. Follow EGA procedure to get access to the data.  Sample BAM files should be downloaded and placed in the `BAMS` directory under the `PROJECTDIR`.

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
BAMDIR=${PROJECTDIR}/BAMS
mkdir -p ${BAMDIR}
# Download BAM files from EGA repository here
```

### Calculate depth of coverage

To make a job array to calculate the depth of coverage for each sample across the bait set coordinates we ran the following commands:

```bash
export PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
export RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
export BAITSET=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed
export LOGDIR=${PROJECTDIR}/logs/depth

mkdir -p ${RESULTSDIR}
mkdir -p ${LOGDIR}

cd ${RESULTSDIR}
#To create the sample list file 
ls -1 ${PROJECTDIR}/BAMS/*.bam >sample_bams.fofn
ls -1 ${PROJECTDIR}/BAMS/*.bam | cut -f 9 -d "/" | sed 's/\.sample\.dupmarked\.bam//' >sample.list
#To create the jobs for submission of sample depth calculation
for f in `cat sample.list`; do bsub -e ${LOGDIR}/depth.${f}.e -o ${LOGDIR}/depth.${f}.o -n 2 -M2000 -R"select[mem>2000] rusage[mem=2000]" "samtools depth -a -b ${BAITSET} -o ${f}.depth.tsv -J -s -@ 2 ${PROJECTDIR}/BAMS/${f}.sample.dupmarked.bam"; done
```
OUTPUT : `${f}.depth.tsv` files for each sample containing 

Subsequently to generate the summaries of the proportion of bases covered with at least a given coverage threshold in steps of 10X, the following commands were run:

```bash
export PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
export RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
export DPSCRIPTDIR=${PROJECTDIR}/scripts/coverage_qc
export LOGDIR=${PROJECTDIR}/logs/depth

cd ${RESULTSDIR}

for f in `cat sample.list`; do bsub -e ${LOGDIR}/count.${f}.e -o ${LOGDIR}/count.${f}.o -M2000 -R"select[mem>2000] rusage[mem=2000]" -q small "${DPSCRIPTDIR}/count_region_coverage.pl ${f}.depth.tsv > ${f}.covstats.tsv"; done
```
OUTPUT : `${f}.covstats.tsv` files for each sample containing the coverage depth statistics for each sample at 10X intervals from 11+(>10X) to 121+.

- Then to summarise the coverage stats per cohort

```bash
export PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
export RESULTSDIR=${PROJECTDIR}/results/qc_plots/depth
export DPSCRIPTDIR=${PROJECTDIR}/scripts/coverage_qc
#Togenerate the cov_stats summary for the cohort
cd ${RESULTSDIR}
#Get the header from the first file
for f in `ls -1 *.covstats.tsv | head -n1`; do head -n1 ${f} | sed 's/^11/Sample\t11/' >cov_stats_summary.tsv; done
grep -v Mean *stats.tsv | sed 's/.covstats.tsv:/\t/' >> cov_stats_summary.tsv

# Generate the tables and plots for the coverage stats summary
cd ${PROJECTDIR}
Rscript ${DPSCRIPTDIR}/02_plot_cov_stats_summary_sortByMeanCov_mod.R
```
The outputs of the above script will be located in the `results/qc_plots/depth` directory. The script will generate the following files:

- `cov_stats_summary.tsv`: file containing the coverage depth statistics for each sample at 10X intervals from 11+(>10X) to 121+
- `summary_cov_stats_ordered.png`: Plot showing the coverage depth statistics for each sample with samples on the Y-axis and proportion of bases within the bait regions with X coverage at given coverage depth intervals, from 11+ to 121+ on the X-axis. Text on the plot indicates the highest interval of coverage depth at which the sample has at least 80% of the bases covered. 
- `Final_cov_stats_summary.tsv`: file containing the list of samples that passed the coverage depth QC, i.e., samples with at least 80% of the bases covered at 20X (21+ or above).

Samples with at least 80% of the bases covered at 20X (21+ or above) passed the coverage depth QC and were considered for somatic variant calling.