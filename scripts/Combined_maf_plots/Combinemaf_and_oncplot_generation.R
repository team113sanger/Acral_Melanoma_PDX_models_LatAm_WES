#!/usr/bin/env Rscript
# Start by reading the libraries that are required
# 0. Set the variables for the Project -------
# 0.1 Set the working directory
here::i_am(file.path("scripts", "Combined_maf_plots", "Combinemaf_and_oncplot_generation.R"))
# Load the libraries
library(renv)
library(here)
library(RColorBrewer)
library(writexl)
#library(data.table)
library(maftools)


# 0 Set the directroies to work ----
wdir<-here()
wdir<-file.path(wdir, "analysis","somatic_variants","release_v4","Combined_2729_3248")
pdx_dir<-file.path(wdir, "PDX")
tum_dir<-file.path(wdir, "Tumour")
# Read the COSMIC cancer gene census file
#cosmic_genesv97<-read.csv(file =file.path(wdir, "COSMIC","cancer_gene_census.v97.genes.tsv" ), header = T, sep = "\t", stringsAsFactors = F )
# Write the list of genes that are in cosmic 
cosmic_genesv97<-read.csv(file =file.path(wdir, "COSMIC","AM_cancer_gene_census.v97.genes_overlap.tsv" ), header = T, sep = "\t", stringsAsFactors = F )


# 1 Merge and plot the data for Combined all_samples  matched & unmatched ----
prefix<-"all_samples"
suffix<-"all_samples"
comball_dir<-file.path(wdir, prefix)

#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, "all_samples", "6633_3248-filtered_mutations_all_allTum_keep.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, "all_samples", "6633_2729-filtered_mutations_all_allTum_keep.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)

#Read the PDX and Tumour TMBs mafs to merge them
pdx_tmb<- as.data.frame(data.table::fread(input = file.path(pdx_dir, "all_samples", "mutations_per_Mb.tsv"), 
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_tmb<- as.data.frame(data.table::fread(input = file.path(tum_dir, "all_samples", "mutations_per_Mb.tsv"),
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_tmb<- rbind(pdx_tmb, tum_tmb)
colnames(comb_tmb)<-c("PDID", "mutations_per_MB")
#Create the directory for the files 
dir.create(comball_dir, recursive = T)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_all_allTum_keep")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )
#Save the TMB 
writexl::write_xlsx(comb_tmb, paste0(comb_maf_fname, "mutations_per_MB.xlsx"),col_names = T, format_headers = T )
#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, "all_samples", "6633_3248-filtered_mutations_all_allTum_keepPA.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, "all_samples", "6633_2729-filtered_mutations_all_allTum_keepPA.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_all_allTum_keepPA")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )


# 2 Read the Combined all_samples matched & unmatched maf file for oncoplots ----
#Try reading the maf
comb_maf<-read.maf(maf=paste0(comb_maf_fname, ".maf"), clinicalData = file.path(wdir, "6633_METADATA_table.tsv") )

#Set the colours for the Tumour or PDX
pcol<-c("#295884", "#FCCDE5")
#pcol<-c("#295884", "#D20C60")
names(pcol)<- unique(comb_maf@clinical.data$sample_type)
#Set the colours for the  PDX  PAssage
passagecol<- c("#CCEBC5", "#FFED6F", "#8DD3C7","#BEBADA","#A6761D")  
names(passagecol)<- unique(comb_maf@clinical.data$PDX_passage )
#Set the colour for tumour of origin 
thiscol<- c("#BC80BD", "#B3DE69",  "#D20C60", "#1B9E77") 
names(thiscol)<- unique(comb_maf@clinical.data$tissue_histology )

clin_col<-list(sample_type=pcol,
               PDX_passage=passagecol,
               tissue_histology=thiscol
               )

#Mutation colours - standard used by oncoplots
#Number of Top Mutated genes by frequency
topn<-40 # Change this to print top 30
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,   #Show title or not 
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
#        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()
#Oncoplot with all the susbstitution types proportions and sorted by Annotation
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS", suffix,"_TiTv.pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not 
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), 
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByMutation = T,
         sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = TRUE #Draw the SNVs by substitution type 
)
dev.off()
# Oncoplot showing COSMIC genes from the CANCER gene census only
#Top mutated genes that are in cosmic 
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,   #Show title or not 
         genes = cosmic_genesv97$Gene.Symbol[cosmic_genesv97$Gene.Symbol %in% comb_maf@gene.summary$Hugo_Symbol[1:topn]], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()

# Oncoplot showing ALLS COSMIC genes from the CANCER gene census only
topn<-40
observed_mut_canc_genes<-cosmic_genesv97$Gene.Symbol[cosmic_genesv97$Gene.Symbol %in% comb_maf@gene.summary$Hugo_Symbol]

##  PLOT the list of genes that present in cosmic and in the topN mutated genes

fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_ALL_top_", topn, "_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not
         genes = observed_mut_canc_genes[1:topn], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1,
         legendFontSize = 1.6, 
         barcode_mar = 6.5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()



# 3 Merge and plot the data for Combined all_samples  matched & unmatched ----
prefix<-"matched_samples"
suffix<-"matched_samples"
comball_dir<-file.path(wdir, prefix)

#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "6633_3248-filtered_mutations_matched_allTum_keep.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "6633_2729-filtered_mutations_matched_allTum_keep.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)

#Read the PDX and Tumour TMBs mafs to merge them
pdx_tmb<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "mutations_per_Mb.tsv"), 
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_tmb<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "mutations_per_Mb.tsv"),
                                          header = F, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_tmb<- rbind(pdx_tmb, tum_tmb)
colnames(comb_tmb)<-c("PDID", "mutations_per_MB")
#Create the directory for the files 
dir.create(comball_dir, recursive = T)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_matched_allTum_keep")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )
#Save the TMB 
writexl::write_xlsx(comb_tmb, paste0(comb_maf_fname, "mutations_per_MB.xlsx"),col_names = T, format_headers = T )
#Read the PDX and Tumour mafs to merge them
pdx_maf<- as.data.frame(data.table::fread(input = file.path(pdx_dir, prefix, "6633_3248-filtered_mutations_matched_allTum_keepPA.maf"), 
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F) 
tum_maf<- as.data.frame(data.table::fread(input = file.path(tum_dir, prefix, "6633_2729-filtered_mutations_matched_allTum_keepPA.maf"),
                                          header = T, sep="\t", stringsAsFactors = F),stringsAsFactors = F)
comb_maf<- rbind(pdx_maf, tum_maf)
#Save the MAF files
comb_maf_fname<-file.path(comball_dir, "6633_2729_3248-filtered_mutations_matched_allTum_keepPA")
data.table::fwrite(comb_maf, paste0(comb_maf_fname, ".maf"), sep = "\t", col.names = T, row.names = F, quote = F)
writexl::write_xlsx(comb_maf, paste0(comb_maf_fname, ".xlsx"),col_names = T, format_headers = T )


# 4 Read the Combined all_samples matched & unmatched maf file for oncoplots ----
#Try reading the maf
comb_maf<-read.maf(maf=paste0(comb_maf_fname, ".maf"), clinicalData = file.path(wdir, "6633_METADATA_table.tsv") )

#Set the colours for the Tumour or PDX
pcol<-c("#295884", "#FCCDE5")
#pcol<-c("#295884", "#D20C60")
names(pcol)<- unique(comb_maf@clinical.data$sample_type)
#Set the colours for the  PDX  PAssage
passagecol<- c("#CCEBC5", "#FFED6F", "#8DD3C7","#BEBADA","#A6761D")  
names(passagecol)<- unique(comb_maf@clinical.data$PDX_passage )
#Set the colour for tumour of origin 
thiscol<- c("#BC80BD", "#B3DE69",  "#D20C60", "#1B9E77") 
names(thiscol)<- unique(comb_maf@clinical.data$tissue_histology )

clin_col<-list(sample_type=pcol,
               PDX_passage=passagecol,
               tissue_histology=thiscol
)

#Mutation colours - standard used by oncoplots
#Number of Top Mutated genes by frequency
topn<-40 # Change this to print top 30
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not 
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()
#Oncoplot with substitution types 
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS", suffix,"_TiTv.pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not 
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'),
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByMutation = T,
         sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = TRUE #Draw the SNVs by substitution type 
)
dev.off()

# Oncoplot showing COSMIC genes from the CANCER gene census only
#Top mutated genes that are in cosmic 
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_", "Top", topn, "_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,   #Show title or not 
         genes = cosmic_genesv97$Gene.Symbol[cosmic_genesv97$Gene.Symbol %in% comb_maf@gene.summary$Hugo_Symbol[1:topn]], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1.6,
         legendFontSize = 1.6, 
         barcode_mar = 5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()

# Oncoplot showing ALLS COSMIC genes from the CANCER gene census only
topn<-40
observed_mut_canc_genes<-cosmic_genesv97$Gene.Symbol[cosmic_genesv97$Gene.Symbol %in% comb_maf@gene.summary$Hugo_Symbol]
fig_fnamepdf<-file.path(comball_dir, paste0("Combined_6633_2729_3248_mutated_gene_Keep_SNVsand_VAF_filt_INDELS_mutsort_ALL_top_", topn, "_COSMICv97cgenes_", suffix,".pdf"))
pdf(file =fig_fnamepdf ,width = 18, height = 16)
oncoplot(comb_maf,
         top = topn,
         altered=F,
         showTitle = T,  #Show title or not
         genes = observed_mut_canc_genes[1:topn], # This is to look only the genes out of the topN most frequently mutated genes that are known cancer genes Cancer gene sensus COSMICv97
         showTumorSampleBarcodes = T, # show sample names or not
         SampleNamefontSize = 1, # sample name font size
         additionalFeatureCex = 0.9,
         genesToIgnore =  maftools:::flags()[1:20], # To ignore specific genes e.g. top 20 FLAG genes 
         #        sampleOrder = c("PD52540c", "PD52540a", "PD52540d"), # This is to ensure the order as the tissue section 
         annotationFontSize = 1,
         legendFontSize = 1.6, 
         barcode_mar = 6.5,
         gene_mar=6,
         clinicalFeatures = c('sample_type','PDX_passage', 'tissue_histology'), #Clinical data being labelled
         sortByMutation = T,
         removeNonMutated = F, # This is to ensur the plotting of all samples
         #sortByAnnotation = TRUE, # Sort by any clinical feature or not
         # colors=varc_cols, # To set the colours of the mutations 
         annotationColor = clin_col, # This is to ensure the colours match between slide and the samples
         draw_titv = FALSE #Draw the SNVs by substitution type 
)
dev.off()





