#!/software/R-4.1.0/bin/Rscript 


#.libPaths(c("/nfs/casm/team113da/users/kw10/lib/R-4.0.3/"))
#.libPaths()
# install.packages(c("viridis",  "cowplot"))
library(ggplot2)
library(viridis)
library(dplyr)
library(reshape2)
library(cowplot)

#Get wdir to rpoject
projdir<-"/lustre/scratch119/casm/team113da/projects/6351_Melanoma_patients_lacking_established_disease_risk_factors_WES/results/qc_plots/depth"
projdir<-"/Users/mdc1/Desktop/team113sc119/projects/6351_Melanoma_patients_lacking_established_disease_risk_factors_WES/results/qc_plots/depth"
manif_dir<-"/Users/mdc1/Desktop/team113sc119/projects/6351_Melanoma_patients_lacking_established_disease_risk_factors_WES/manifests/" # Location of manifest folder
manif<-read.csv(file=file.path(manif_dir,"casm",  "Submitted_cgp_contact_manifest_Mcompilation_info.txt"), header = T, stringsAsFactors = F, sep="\t")

setwd(projdir)

slist<-read.csv(file = file.path(projdir, "sample.list"), header = F, stringsAsFactors = F, sep="\t")
slist<-slist$V1

header<-read.csv(file = file.path(projdir, "cov_stats_summary.tsv"), header = F, stringsAsFactors = F, sep="\t")
header<-header[1,]
data<-NULL
i<-NULL
for (i in 1:length(slist)){
    temp<-read.csv(file = file.path(projdir, paste(slist[i], ".covstats.tsv", sep="")), header = F, stringsAsFactors = F, sep="\t")
    data<-rbind(data, c(as.vector(t(temp[2,]))))  
}
data<-cbind(slist,data)
data<-as.data.frame(data)
colnames(data)<-as.vector(t(header[1,]))
data$`11+`<-as.numeric(data$`11+`)
data$`21+`<-as.numeric(data$`21+`)
data$`31+`<-as.numeric(data$`31+`)
data$`41+`<-as.numeric(data$`41+`)
data$`51+`<-as.numeric(data$`51+`)
data$`61+`<-as.numeric(data$`61+`)
data$`71+`<-as.numeric(data$`71+`)
data$`81+`<-as.numeric(data$`81+`)
data$`91+`<-as.numeric(data$`91+`)
data$`101+`<-as.numeric(data$`101+`)
data$`111+`<-as.numeric(data$`111+`)
data$`121+`<-as.numeric(data$`121+`)
data$`131+`<-as.numeric(data$`131+`)
data$`141+`<-as.numeric(data$`141+`)
data$`151+`<-as.numeric(data$`151+`)
data$Cov_Mean<-as.numeric(data$Cov_Mean)
#file="cov_stats_summary.tsv"
#data <- read.table(file, sep = "\t", header = T, check.names = F)
head(data)
mode(data)
class(data)
mode(data$Cov_Mean)
class(data)


data1 <- data %>% select(-Cov_Mean)
data1.melt <- melt(data1)

colnames(data1.melt) <- c("Sample", "Coverage", "Percent")

labels <- data1.melt %>% group_by(Sample) %>% filter(Percent > 80)  %>% top_n(n=-1, wt=Percent)

data2 <- data %>% select(Sample, Cov_Mean)
order <- data2 %>% arrange(Cov_Mean) %>% select(Sample)
order$Sample

data2.melt<-melt(data2)
colnames(data2.melt)<-c("Sample","Variable","Coverage")
data2.melt$Sample<-factor(data2.melt$Sample, levels=order$Sample)
data1.melt$Sample<-factor(data1.melt$Sample, levels=order$Sample)
str(data2.melt)
str(data1.melt)

plot<-ggplot(data1.melt, aes(Coverage, Sample, fill = Percent))+
  geom_tile()+
  scale_fill_viridis(option="turbo", direction=-1, begin=0.15, end=0.85)+
  theme(axis.text=element_text(size=11), axis.title=element_blank(), axis.ticks=element_blank())+
  scale_x_discrete(position = "top", expand=c(0,0))+
  geom_text(data=labels, aes(label=Percent))

quants <- as.numeric(quantile(data2$Cov_Mean ,c(.50, 0.25, .05, .01)))
quants1 <- which.max(sort(data2$Cov_Mean) > quants[1]) + 0.5
quants2 <- which.max(sort(data2$Cov_Mean) > quants[2]) + 0.5
quants3 <- which.max(sort(data2$Cov_Mean) > quants[3]) + 0.5
quants4 <- which.max(sort(data2$Cov_Mean) > quants[4]) + 0.5
quants.pos <- c(quants1, quants2, quants3, quants4)
quants
quants.pos

mean_cov_plot<-ggplot(data2.melt, aes(Variable, Sample, fill=Coverage)) + geom_tile() + 
  scale_fill_viridis(option="turbo", direction=-1, begin=0.15, end=0.85) + 
  scale_x_discrete(position = "top", expand=c(0,0)) + 
  geom_text(aes(label=Coverage)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=12))

png("summary_cov_stats_ordered.png", width=1000, height=2000)
plot_grid(plot, mean_cov_plot, align = "h", rel_widths = c(5,1))
dev.off()


#Make the final plot with the addition of sample batch and place of origin
data3<- manif %>% select(Supplier_Name, seq_batch)
source_data<- manif %>% select(Supplier_Name, sample_source)
#Keep only the ones that were sequenced
data3<- data3[data3$Supplier_Name %in% data$Sample,]
source_data<- source_data[source_data$Supplier_Name %in% data$Sample,]

colnames(data3)[1]<-"Sample"
colnames(source_data)[1]<-"Sample"
#Then re-order the tables in the same order as the matrix with the coverage data
data3<- data3[match( as.character(data2.melt$Sample), data3$Sample),]
data3$Sample<- factor(data3$Sample, levels = levels(data2.melt$Sample)) #change the variable type to factor
source_data<- source_data[match( as.character(data2.melt$Sample), source_data$Sample),]
source_data$Sample<- factor(source_data$Sample, levels = levels(data2.melt$Sample))
#Add the variable column so it can be plotted all in a single tile column with geom_tile
data3$Variable<-rep("seq_batch", n=dim(data3)[1])
source_data$Variable<-rep("sample_origin", n=dim(source_data)[1])

#Generate the plot for the batch
batch_plot<-ggplot(data3, aes(Variable, Sample, fill=seq_batch)) +
  geom_tile() + 
  scale_fill_manual(values = c("#CAB2D6","#6A3D9A")) + 
  scale_x_discrete(position = "top", expand=c(0,0)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=12))
batch_plot
#Generate the plot for the batch
orig_plot<-ggplot(source_data, aes(Variable, Sample, fill=sample_source)) +
  geom_tile() + 
  scale_fill_brewer(type="Seq", palette="Paired") + 
  scale_x_discrete(position = "top", expand=c(0,0)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=12))
orig_plot

#Plot the four plots containg the coverage matrtix, mean coverage and the sequencing batch and sample origin information 
png("summary_cov_stats_ordered_wbatch.png", width=1100, height=2000)
plot_grid(plot,batch_plot,orig_plot, mean_cov_plot, align = "h", rel_widths = c(5,1,1,1),nrow = 1,ncol = 4)
dev.off()


#Get the samples which had a coverage across 80% of the baits positions 
samp_20cov<-data1.melt[data1.melt$Coverage=="21+", ]
samp_20cov<-samp_20cov<-samp_20cov[ samp_20cov$Percent>79.99, ]
write.table(samp_20cov$Sample, file=file.path(projdir, "cov_pass_sample.list"),quote = F, row.names = F, col.names = F, sep = "\t")


#Add the important information to the manifest about the seuqenced samples successfully and the samples that passed Coverage QC
manif<- manif  %>% mutate(sequenced_status= Supplier_Name %in% data$Sample )
manif<- manif  %>% mutate(cov_qc_status= ifelse(Supplier_Name %in% samp_20cov$Sample, "Pass", "Failed") )

#Get only the pass samples 
fin_manif<- manif[ manif$cov_qc_status %in% "Pass",  ]
table(fin_manif$seq_batch)
table(fin_manif$sample_source)

#Get only the pass samples 
fail_manif<- manif[ manif$cov_qc_status %in% "Failed",  ]
table(fail_manif$seq_batch)
table(fail_manif$sample_source)


#Write the Final Manifest
write.table(manif, file = file.path(manif_dir,paste0("Submitted_cgp_contact_manifest_Mcompilation_with_COV_QC_passed_info.txt")), quote = F, col.names = T, row.names = F, sep="\t")


