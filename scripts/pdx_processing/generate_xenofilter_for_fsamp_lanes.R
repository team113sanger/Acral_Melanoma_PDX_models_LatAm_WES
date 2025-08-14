#!/software/R-4.1.0/bin/Rscript 
# generate_xenofilter_for_fsamp_lanes.R
#Script to generate the jobs to filter the failed 3 samples due to read numbers and samtools issue 

manifdir<-"/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/manifests"
outdirscripts<-"/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/scripts"
outdir<-"/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/BAMS/WES_xfilt"
xfilt_logdir<-"/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/logs/xfilt_logdir"


#Read the table with the information for the samples that need filtering by Lane 
table<-read.csv(file = file.path(manifdir, "Failed_xfilter_samples_WES_RG_info.txt"), header = T, sep = "\t", stringsAsFactors = F)
# table<-read.csv(file = file.path("/Users/mdc1/Documents/Projects/PDX_models_Latin_America/6633_PDX_models_Latin_America_WES/manifests", "Failed_xfilter_samples_WES_RG_info.txt"), header = T, sep = "\t", stringsAsFactors = F)
#

#PAHTS to the programs and things that are required
run_xenofilter<-"/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/scripts/bam_Xenofilter_rg.R"
refname_nod1<-"NOD_PDXV1"

# Get the commands to submit the jobs
cmds_nodv1<-NULL
for(i in 1:dim(table)[1]){
  print(paste0("Processing sample ", table$PDID[i]," ", table$Hum_PU[i], " NodV1"))
  tempsamp<-NULL
  cmd<-NULL
  sdoutf<- NULL
  serrf<- NULL
  thbam<-NULL #Temp human bam
  tmbam<-NULL #Temp mouse bam
  tempsamp<- table$PDID[i]
  #Assing job output names 
  sdoutf<-file.path(xfilt_logdir, paste("xenofilter_log_nodv1_",tempsamp,"_", i,".o", sep = ""))
  serrf<-file.path(xfilt_logdir, paste("xenofilter_log_nodv1_",tempsamp,"_", i,".e", sep = ""))
  #Sample analysing
  temp<-table[i, ]
  thbam<- temp$unfilt_run_lane_tid_bam_path[1]
  tmbam<- temp$bam_run_lane_tid_path_nodv1[1]
  #Command for Nodv3 per sample file and index Min 64GBs RAM
  cmd<-as.character(paste("bsub -q normal -M 64000 -R",shQuote("select[mem>64000] rusage[mem=64000] span[hosts=1]") , " -n 2 -o ", sdoutf," -e ",serrf,
                          " ' module load R/4.1.0  ; Rscript ", run_xenofilter, " --sample_name ", tempsamp,
                          " --human_bam ", thbam,
                          " --mouse_bam ", tmbam,
                          " --rg_id ", temp$Hum_PU[1],
                          " --outdir ", file.path(outdir,refname_nod1),
                          " --ncpu ", "1",
                          " '",
                          sep=""))  
  #Append the command to the list of NODv1 merging
  cmds_nodv1<- c(cmds_nodv1, cmd)
}

#Generate the .sh files with the submissions 
write.table(c("#!/bin/sh", cmds_nodv1), file=file.path(outdirscripts, paste("xenofilter_nodv1_pdx_xfilt_fsamp_jobs.sh", sep="")), quote = F, col.names = F, row.names = F, sep = "\n")
