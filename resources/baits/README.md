# Bait set

This project uses Agilent SureSelect V5 bait coordinates which were obtained from Agilent and converted to GRCh38 using `liftOver` from referred to as  `SureSelect_Whole_Human_Exome_v5_GRCh38_liftover_160` below.

## Bait set with just the canonical chromosomes

For coverage depth evenness calculation, we used a sorted BED file that contains only canonical chromosomes. To obtain this file we ran the following commands:

```bash
module load bedtools/2.29.0

awk '$0!~/_/' SureSelect_Whole_Human_Exome_v5_GRCh38_liftover_160.bed | bedtools sort -g canonical-chrs.txt | bedtools merge -i - > GRCh38_WES5_canonical.bed

```

⚠️ Note: We use the fact that the canonical chromosomes do not contain the `_` character, whilst the non-canonical chromosomes do.

## Padded Bait set 

Most analyses require a 100-nt padded bait set so that reads that overlap baits will be captured. This process requires the chromosome lengths.  These are pulled from the genome DNA FASTA index file located `

This is generated from the canonical BED file as:

```bash
bedtools slop -b 100 -g GRCh38_full_analysis_set_plus_decoy_hla/genome.fa.fai GRCh38_full.fasta.fai -i ./GRCh38_WES5_canonical.bed > GRCh38_WES5_canonical_pad100.bed
```

## Generate file with regions merged

The above file has 7MB of overlap, for analysis, use a BED file with non-overlapping regions a bedtools merge step is needed:

```bash 

bedtools slop -b 100 -g /nfs/cancer_ref02/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa.fai -i ./GRCh38_WES5_canonical.bed | bedtools merge -i | sort -Vk 1,2 > GRCh38_WES5_canonical_pad100.merged.bed

```

To share the files they were compressed and available on the `resources/baits` directory.

```bash
gzip GRCh38_WES5_canonical_pad100.merged.bed
gzip GRCh38_WES5_canonical.bed

```