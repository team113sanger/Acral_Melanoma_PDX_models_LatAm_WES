# Somatic variant calling parameters, annotation filtering and results

## Overview
This document describes the steps taken to identify somatic single nucleotide variants (SNVs) and short INDEL mutations present on the whole exome of the biopsied samples. This was done using CaVEMAN and Pindel, followed by variant consequence prediction with VEP.  All samples were checked for evenness of coverage before the variant calling was performed (see [documentation/Coverage_depth_check.md](./Coverage_depth_check.md)).   

The somatic calling, flagging and annotation was performed with an internal pipeline. The steps, and parameters used for the variant calling are described below including some example scripts. However, **Filtering, MAF conversion and plotting** were performed using the commands and scripts mentioned below on their respective section

## Results

MAF file containing the summary of the results of the variant calling can be found in FigShare Project [here](https://figshare.com/projects/A_post-zygotic_disruptive_germline_NFATC1_variant_in_a_patient_with_segmental_cherry_angiomas/254267). 

Raw VCF files containing the flagged and VEP annotated variant calls can be downloaded from EGA study Accession **EGAS00001008212**

Information on how and where to copy the  result MAF file can be found in the [`analysis/variants_combined/version1/README.md`](../analysis/variants_combined/version1/README.md) file.

## Sample Pairs used for somatic calling 

The list of sample pairs used for variant calling using the somatic callers can be found in the file [`metadata/8117-biosample_manifest-completed.tsv`](../metadata/8117-biosample_manifest-completed.tsv). 
 
## Required Environment variables and software

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.3.3`[**here**](https://cran.r-project.org/)
- **samtools**: `v1.14` [**here**](https://github.com/samtools/samtools)
- **bcftools**: `1.9` [**here**](https://github.com/samtools/bcftools/)
- **CaVEMan**: `1.17.2` [**here**](https://github.com/cancerit/CaVEMan)
- **cgpCaVEManwrapper**:`1.17.2` [**here**](https://github.com/cancerit/cgpCaVEManWrapper)
- **cgpPindel**: `v3.5.0` [**here**](https://github.com/cancerit/cgpPindel)
- ENSEMBL VEP version `103`[**here**](http://feb2021.archive.ensembl.org/info/docs/tools/vep/index.html)
- The repositories with the scripts used for variant QC and VCF to MAF conversion have been added as **submodules** to this repository:
    - [**QC**](https://github.com/team113sanger/dermatlas_analysis_qc) `v0.4.3`
    - [**MAF**](https://github.com/team113sanger/dermatlas_analysis_maf) `v0.5.4` 
    - To clone the submodules used if not done when cloning the repository, run the following command:
```bash 
git submodule update --init --recursive
```

:warning: **IMPORTANT NOTES** :warning:
- The scripts below are written to be run in our internal servers, which have Ubuntu 22.04.5 and submit jobs using `lsf bsub`. They are written to show an example of the commands used, this is because the calling was performed with an internal pipeline that uses CaVEMan and Pindel which utilise the same singularity images can be made from the repositories mentioned above.  Paths and environment variables would need to be adjusted to run in a different environment.
- `cgpCaVEManwrapper` and `cgpPindel` were downloaded and used as singularity images within our internal pipelines using `bpipe`[https://docs.bpipe.org/](https://docs.bpipe.org/). 
- bpipe running parameters used for CaVEMan and Pindel can be found in the [bpipe_CaVEMan_Pindel_running_params.md](./bpipe_CaVEMan_Pindel_running_params.md) file.
- Path, or scripts module call modifications may need to be made to run in a different environment.

### Somatic variant calling

The variant calling was performed using `CaVEMan` and `cgpPindel` somatic callers. To see the samples pairs used as tumour-normal for the calling see [metadata/8117-biosample_manifest-completed.tsv](../metadata/8117-biosample_manifest-completed.tsv). Metadata for the samples can be located in the[metadata/8117_2744_metadata.tsv](../metadata/8117_2744_metadata.tsv) table. Functional annotation was done with ENSEMBL v103 Variant Effect Predictor (VEP). 

####  **STEP 1- CAVEMAN v1.15.1 SNV calling**

Running parameters for the Caveman calling can be found inside the [run_Cavemanwrapper_1.17.2.sh](../scripts/offpipe_calling/run_Cavemanwrapper_1.17.2.sh) script.

**NOTE**
Set the `PROJECTDIR` variable to the path where the repository was cloned into and run the following commands in the terminal:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
PREFIX="${STUDY}_${PROJECT}"

cd ${PROJECTDIR}/scripts

# BAM Directory for WES data
BMDIR=${PROJECTDIR}/bams

mkdir -p ${PROJECTDIR}/analysis/CAVEMAN
cd ${PROJECTDIR}/analysis/CAVEMAN
mkdir -p logs;
num=0; 
for Tumour_PDID in `cut -f 1 ${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv |cut -f 2`; 
bsub -e ./logs/caveman.${Tumour_PDID}-vs-${Normal_PDID}.e -o ./logs/caveman.${Tumour_PDID}-vs-${Normal_PDID}.o \
-J"cavemanrun[$num]" -n 7 -M40000 -R"select[mem>40000] rusage[mem=40000] span[hosts=1]"\
 -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_Cavemanwrapper_1.17.2.sh ${Tumour_PDID:?unset} ${Normal_PDID:?unset} 6 ${BMDIR:?unset} "; 
done

```
#### STEP 2-  CAVEMAN v1.17.2 flagging - cgpFlagCaVEMan.pl

For this step the flagging .ini files are required, these contain the parameters for the flagging of the variants. A copy of these are shared within this repository in the resources folder. The files are:
- [**flag.to.vcf.convert.ini**](../resources/caveman/flag.to.vcf.convert.ini) -  contains conversions of flags to FLAG ID
- [**flag.vcf.config.ini**](../resources/caveman/flag.vcf.config.ini) - specifies which flags are cutoffs to use for WXS, WGS, etc. 

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

# BAMDIR
BMDIR=${PROJECTDIR}/bams
# Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
# Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
CAVECONFIG=/lustre/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/caveman/flag.vcf.config.ini

cd ${PROJECTDIR}/analysis/CAVEMAN
num=0; 
for Tumour_PDID in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
Normal_PDID=`grep ${Tumour_PDID} ${SAMPLETSV} |cut -f 2`;
PAIR=${Tumour_PDID}"_vs_"${Normal_PDID}; 
bsub -e ./logs/caveman.cgpFlag.${Tumour_PDID}-vs-${Normal_PDID}.e -o ./logs/caveman.cgpFlag.${Tumour_PDID}-vs-${Normal_PDID}.o \
-J"cavemanflag[$num]" -n 4 -M8000 -R"select[mem>8000] rusage[mem=8000] span[hosts=1]"\
 -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_cgpFlagCaVEManvv1.17.2.sh ${CAVECONFIG:?unset} v1 ${Tumour_PDID} ${Normal_PDID} ${PAIR:?unset} ${BMDIR:?unset} ${CAVEMANDIR:?unset} "; 
done
```

Then the VCFs were compressed and indexed using the following commands:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

# Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
# Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
cd ${CAVEMANDIR}
cat ${SAMPLETSV:?unset} | cut -f 1,2 | sed 's/\t/_vs_/g' >samples.list

module load bcftools-1.9/python-3.11.6
find ${PROJECTDIR}/analysis/CAVEMAN/*/* -name "*_DNA.flag.vcf" -type f -exec bgzip {} \;
find ${PROJECTDIR}/analysis/CAVEMAN/*/* -name "*_DNA.flag.vcf.gz" -exec tabix -p vcf {} \;

```
#### STEP 3- CGP CAVEMAN v1.17.2  VEP103 annotation - run_vep_caveman_vdermEns103.target.custom.grch38.sh

To predict the effect of the identified variants, ENSEMBLs' VEP v103 was used with additional custom annotations regarding their presence in external datasets such as:
    - gnomAD [`**v3.1**`](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/)
    - dbSNP [`**v155**`](https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz) Instructions on how the file was parsed can be found [`resources/dbsnp`](../resources/dbsnp)
    - COSMIC [`**v97**`](https://cancer.sanger.ac.uk/cosmic/download/cosmic)
    - ClinVar [`**20220115**`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220115.vcf.gz)

To run the VEP annotation the following commands were used :

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
#Caveman directory
CAVEMANDIR=${PROJECTDIR}/analysis/CAVEMAN
cd ${CAVEMANDIR}
# List of caveman vcf files
ls -1 ${CAVEMANDIR}/*/*.caveman_c.flag.vcf.gz >samples.list
#Run the VEP annotation
bash ${PROJECTDIR}/scripts/VEP/run_vep_pindel_Ens103.grch38.sh ./samples.list

#Then index the VCFs
module load bcftools-1.9/python-3.11.6
find ./*/* -name "*.caveman_c.flag.vep.vcf.gz" -exec tabix -p vcf {} \;
```

#### STEP 4- Pindel calling   -CGP Pindel v3.5 calling 

The list of FLAGS for cgpPindel v3.5 is required to run the calling. A copy of the rules file is shared with this repository in the resources folder. The file is:
- [**WXS_Rules.lst**](../resources/pindel/WXS_Rules.lst)
- [**empty_germline_bed_for_caveman.bed**](../resources/caveman/example/empty_germline_bed_for_caveman.bed)

To run the calling 
```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

# BAMDIR
BMDIR=${PROJECTDIR}/bams
#Sample pairs file
SAMPLETSV=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL
# path to germline indels bedfile
export GERM_INDEL=${PROJECTDIR}/resources/caveman/example/empty_germline_bed_for_caveman.bed
#Pindel Rules file - copy in ${PROJECTDIR}/resources/pindel
export CGPINDRULES=/lustre/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/pindel/WXS_Rules.lst

mkdir -p ${PINDELDIR}
cd ${PINDELDIR}
mkdir -p logs;
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} |cut -f 2`; 
bsub -e logs/pindel.$f-vs-$g.e -o logs/pindel.$f-vs-$g.o \
-J"pindelrun[$num]" -n 6 -M40000 -R"select[mem>40000] rusage[mem=40000] \
span[hosts=1]" -q long "bash ${PROJECTDIR}/scripts/offpipe_calling/run_cgpPindel_3.5.sh ${f:?unset} ${g:?unset} 6 ${BMDIR:?unset} ${CGPINDRULES:?unset} "; 
done

```
#### STEP 5-  CGP Pindel v3.5  -  VEP annotation

To predict the effect of the identified variants, ENSEMBLs' VEP v103 was used with additional custom annotations regarding their presence in external datasets such as:
- gnomAD [`v3.1`](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/)
- dbSNP [`v155`](https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz) Instructions on how the file was parsed can be found [`resources/dbsnp`](../resources/dbsnp)
- COSMIC [`v97`](https://cancer.sanger.ac.uk/cosmic/download/cosmic)
- ClinVar [`20220115`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220115.vcf.gz)

To run the VEP annotation the following commands were used:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

# BAMDIR
BMDIR=${PROJECTDIR}/bams
SAMPLETSV=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
PINDELDIR=${PROJECTDIR}/analysis/PINDEL
# path to germline indels bedfile
export GERM_INDEL=${PROJECTDIR}/resources/caveman/example/empty_germline_bed_for_caveman.bed
#Pindel Rules file - copy in resources/pindel
export CGPINDRULES=${PROJECTDIR}/resources/pindel/WXS_Rules.lst

cd ${PINDELDIR}

#to create the file with the sample list for analysis (taking the sample pairs names
ls -1 ${PINDELDIR:?unset}/*/*.flagged.vcf.gz >flag_samples_files.list

bash ${PROJECTDIR}/scripts/VEP/run_vep_pindel_Ens103.grch38.sh ./flag_samples_files.list

#Then sort the naming Tumour_vs_Normal to Tumour only output to allow easy conversion to MAF
num=0; 
for f in `cut -f 1 ${SAMPLETSV}`;  do 
echo $num; 
let num=num+1; 
g=`grep $f ${SAMPLETSV} | cut -f 2`;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vcf.gz ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz ; 
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vcf.gz.tbi ${PINDELDIR:?unset}/${f}/${f}.pindel.flagged.vcf.gz.tbi ;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vep.vcf.gz ${PINDELDIR:?unset}/${f}/${f}.pindel.vep.vcf.gz ;
mv ${PINDELDIR:?unset}/${f}/${f}.flagged.vep.vcf.gz.tbi ${PINDELDIR:?unset}/${f}/${f}.pindel.vep.vcf.gz.tbi ;
cp -R ${PINDELDIR:?unset}/${f} ${PROJECTDIR}/analysis/pindel_files/
done

```

### Post processing filter, MAF file generation and plotting

To reproduce this steps download the raw annotated VCFs from EGA study Accession **EGAS00001008212** and place them in the `analysis/caveman_files/PD*` and `analysis/pindel_files/PD*` directories respectively. The VCF files should be named as `<SAMPLE>.caveman_c.vep.vcf.gz` and `<SAMPLE>.pindel.vep.vcf.gz` respectively, where `<SAMPLE>` is the sample name.

#### Adding annotations for dbSNP155 common variants  - CavEMan

We retain variants that passed the Flagging criteria from CaVEMAn and kept variants that were located within the baitset regions.  To do this we ran 
the following commands, using the script `select_vcf_pass_ontarget.sh` from the [`QC`](../scripts/QC) repository. This script filters the VCF files to retain only those variants that are on-target and have the PASS flag.

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

cd ${PROJECTDIR}/analysis

BEDFILE=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed

cd ${PROJECTDIR}/analysis/caveman_files
# First get the on-target PASS variants from the CaVEMan  VCF file 

for f in */*.caveman_c.vep.vcf.gz; do nice -n 5 bash ${PROJECTDIR:?unset}/scripts/QC/select_vcf_pass_ontarget.sh ${f} ${BEDFILE:?unset}; done
```
OUTPUT:
- **`<SAMPLE>.caveman_c.vep.filt.vcf.gz`** files are created in the `analysis/caveman_files/` directories. These files contain the variants that passed the filtering criteria.

Subsequntly, to add the information regarding the presence of variants in dbSNP155, the following steps were performed, using script `add_commonSNPs2vcf.sh` from the [`QC`](../scripts/QC) repository. 

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

BEDFILE=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed

cd ${PROJECTDIR}/analysis/caveman_files
# Now add dbSNP common annotations to the filtered VCF from above
for f in `dir -1 |grep PD`; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p $PROJECTDIR -o ./$f -v $f/$f.caveman_c.vep.filt.vcf.gz; done
 
# This creates *snpflagged.vcf.gz files.
```
OUTPUT:
- **`<SAMPLE>.caveman_c.vep.filt.snpflagged.vcf.gz`** files are created in the `analysis/caveman_files/PD*` directories. These files contain the variants that passed the filtering criteria and have been annotated with dbSNP155 common variants information.


#### Adding annotations for dbSNP155 common variants  - Pindel

Similar to what was done for CaVEMAn, we retained variants that passed the Flagging criteria from Pindel and kept variants that were located within the baitset regions. To do this we ran the following commands, using the script `select_vcf_pass_ontarget.sh`

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

cd ${PROJECTDIR}/analysis/pindel_files

#  First get the on-target PASS variants from the Pindel files 
# The output is *pindel.vep.filt.vcf.gz
BEDFILE=${PROJECTDIR}/resources/baits/SureSelect_Human_All_Exon_V6_plusUTR_GRCh38_liftover.bed

for f in */*pindel.vep.vcf.gz; do bash ${PROJECTDIR}/scripts/QC/select_vcf_pass_ontarget.sh ${f} ${BEDFILE}; done
 
```
OUTPUT:
- **`<SAMPLE>.pindel.vep.filt.vcf.gz`** files are created in the `analysis/pindel_files/PD*` directories. These files contain the variants that passed the filtering criteria.

Subsequently, to add the information regarding the presence of variants in dbSNP155, the following steps were performed, using script `add_commonSNPs2vcf.sh` from the [`QC`](../scripts/QC) repository.

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744

cd ${PROJECTDIR}/analysis/pindel_files

# Now add dbSNP common annotations to the filtered VCF from above
# This creates *snpflagged.vcf.gz files.
 for f in `dir -1 |grep PD`; do echo $f; bash ${PROJECTDIR}/scripts/QC/add_commonSNPs2vcf.sh -p ${PROJECTDIR} -o ./${f} -v ${f}/${f}.pindel.vep.filt.vcf.gz; done

```
OUTPUT:
- **`<SAMPLE>.pindel.vep.filt.snpflagged.vcf.gz`** files are created in the `analysis/pindel_files/PD*` directories. These files contain the variants that passed the filtering criteria and have been annotated with dbSNP155 common variants information. These files are used to create the MAF files in the next step.

#### Get the list of CaVEMan and Pindel files 

We get the list of all the VCFs and their paths into a single file `caveman_pindel_vcfs_all.list` to be used for the MAF generation. The following commands were used to create the directory and the lists of VCF files to be used for the MAF generation. 

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
i=1
mkdir -p ${PROJECTDIR}/analysis/variants_combined/version${i}

cd ${PROJECTDIR}/analysis/variants_combined/version${i}

# Generate the file witha ll the paths of the VCF files
dir -1  ${PROJECTDIR}/analysis/caveman_files/*/*.snpflagged.vcf.gz ${PROJECTDIR}/analysis/pindel_files/*/*snpflagged.vcf.gz | grep -f ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-analysed_all_tum.txt > caveman_pindel_vcfs_all.list
 
```

#### Get the MAF files

To make MAF files and plots, the VCF files containing the SNV, MNV, INDELs per samples were combined and converted into MAF file. In this step all variants in coding and splice sites regions (including synonymous variants) were included. This was done using the `reformat_vcf2maf.pl` scripts form the MAF repository. The following commands were used :

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
i=1

# Variables
VCFLIST=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_vcfs_all.list
BUILD="GRCh38"
SAMPLELIST=${PROJECTDIR}/metadata/8117_2744-analysed_all.tsv
TEMP_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/caveman_pindel_pass_keep.maf
FINAL_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/${STUDY}_${PROJECT}-filtered_mutations_matched_allTum_keep.maf
SCRIPTDIR=${PROJECTDIR}/scripts/MAF

cd ${PROJECTDIR}/analysis/variants_combined/version${i}

# Transform the filtered VCFs into Combined MAFs
# Convert VCFs to maf
${SCRIPTDIR}/reformat_vcf2maf.pl --vcflist ${VCFLIST} --build ${BUILD} --pass --af_col "gnomAD_AF" --sample_list ${SAMPLELIST} --dbsnp_filter > ${TEMP_MAF}

#Get the keep variants - all variants in coding and splice sites regions (including synonymous variants)
head -n1 ${TEMP_MAF} > ${FINAL_MAF}
awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' ${TEMP_MAF} >> ${FINAL_MAF}
rm ${TEMP_MAF}

# Transform to Excel 
Rscript ${SCRIPTDIR}/maf2xlsx.R ${FINAL_MAF}

```

OUTPUT:
The script outputs a MAF and an excel file, both contain all variants in coding and splice sites regions (including synonymous variants)
- [`analysis/variants_combined/version1/8117_2744-filtered_mutations_matched_allTum_keep.maf`](../analysis/variants_combined/version1/8117_2744-filtered_mutations_matched_allTum_keep.maf)
- [`analysis/variants_combined/version1/8117_2744-filtered_mutations_matched_allTum_keep.xlsx`](../analysis/variants_combined/version1/8117_2744-filtered_mutations_matched_allTum_keep.xlsx)


#### Summary plots

Finally to generate the summary plots of with all of the variants found, we ran the `03_somatic_var_plots.R` script.  This script takes the MAF file generated above and creates a summar tileplot and Figure 1B. The script is located in the [`scripts/var_plots`](../scripts/var_plots) directory. The following commands were used to run the script:

```bash 
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
STUDY=8117
PROJECT=2744
i=1
# Variables
export FINAL_MAF=${PROJECTDIR}/analysis/variants_combined/version${i}/${STUDY}_${PROJECT}-filtered_mutations_matched_allTum_keep.maf
export SCRIPTDIR=${PROJECTDIR}/scripts/var_plots

# Results dir
cd ${PROJECTDIR}/analysis/variants_combined/version${i}

Rscript ${SCRIPTDIR}/03_somatic_var_plots.R

```

OUTPUT:
- [`results/somatic_var_plots/Cherry_angioma_HS_som_mut_oncoplot.pdf`](../analysis/variants_combined/version1/results/somatic_var_plots/Cherry_angioma_HS_som_mut_oncoplot.pdf) : oncoplot with the summary of the somatic variants found in the samples.
- [`results/somatic_var_plots/NFATC1_loliplot.pdf`](../analysis/variants_combined/version1/results/somatic_var_plots/NFATC1_loliplot.pdf) : loliplot with the summary of the NFATC1 variants found in the samples. This figure is used in the manuscript as **Figure 1B**. See [`results/figures/Fig1B_NFATC1_lollipop_Maftools.pdf`](../results/figures/Fig1B_NFATC1_lollipop_Maftools.pdf)

