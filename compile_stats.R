# compile_stats

#Read in a Roary gene_presence_absence.Rtab file and compile some stats
#How many core, soft core, shell, cloud genes? Which genes belong to each category?


#Thresholds
#Core threshold (99% <= strains <= 100%)
core_thresh <- 0.99

#Soft core threshold (95% <= strains < 99%)
score_thresh <- 0.95

#Shell threshold (15% <= strains < 95%)
shell_thresh <- 0.15

#Cloud threshold (0% <= strains < 15%)


#97+108 merge (99):
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/gene_pres_abs_97_108.Rtab

#97+108 Roary:
#/home/nadege/Desktop/roary/97_108_refined_meta/gene_presence_absence.Rtab

#97+108 merge (100):
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/gene_pres_abs_97_108.Rtab
#97+108 merge (99):
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/gene_pres_abs_97_108.Rtab


inFile <- "/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/gene_presence_absence_97_108.Rtab"

df <- read.csv(file = inFile, header = T, row.names = 1, sep=",")

#Store the gene names (row names)
gene_names <- rownames(df)

#How many isolates is each gene in?
num_isolates <- apply(df,1,sum)

#How many total isolates?
total_isolates <- length(colnames(df))

#Gene lists
core_genes <- c()
score_genes <- c()
shell_genes <- c()
cloud_genes <- c()

#Gene counts
num_core <- 0
num_score <- 0
num_shell <- 0
num_cloud <- 0

#Go through each gene and sort it into a category + update counts
for (i in (1:length(num_isolates)))
{
  spread <- num_isolates[i]/total_isolates
  
  if (spread < shell_thresh)
  {
    num_cloud <- num_cloud + 1
    cloud_genes[num_cloud] <- gene_names[i]
  }
  
  else if (spread < score_thresh)
  {
    num_shell <- num_shell + 1
    shell_genes[num_shell] <- gene_names[i]
  }
  
  else if (spread < core_thresh)
  {
    num_score <- num_score + 1
    score_genes[num_score] <- gene_names[i]
  }
  
  else
  {
    num_core <- num_core + 1
    core_genes[num_core] <- gene_names[i]
  }
}


#Print the stats file
stats_str <- 
sprintf(
"\tCore genes (%d%% <= strains <= 100%%)\t%d
\tSoft core genes (%d%% <= strains < %d%%)\t%d
\tShell genes (%d%% <= strains < %d%%)\t%d
\tCloud genes (0%% <= strains < %d%%)\t%d
\tTotal genes (0%% <= strains <= 100%%)\t%d", 
        core_thresh*100, num_core,
        score_thresh*100, core_thresh*100, num_score,
        shell_thresh*100, score_thresh*100, num_shell,
        shell_thresh*100, num_cloud,
        length(gene_names))

#97+108 merge (99):
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/summary_stats.txt

#97+108 Roary:
#/home/nadege/Desktop/metagenome/comparison_to_cd-hit-est_attempt/summary_stats.txt

#97+108 merge (100):
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/summary_stats.txt

write(stats_str,file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/summary_stats.txt")


#97+108 merge (99):
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/core_genes.txt
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/soft_core_genes.txt
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/shell_genes.txt
#/home/nadege/Desktop/metagenome/cd-hit-est_attempt/cloud_genes.txt

#97+108 Roary:
#/home/nadege/Desktop/metagenome/comparison_to_cd-hit-est_attemp/core_genes.txt
#/home/nadege/Desktop/metagenome/comparison_to_cd-hit-est_attemp/soft_core_genes.txt
#/home/nadege/Desktop/metagenome/comparison_to_cd-hit-est_attemp/shell_genes.txt
#/home/nadege/Desktop/metagenome/comparison_to_cd-hit-est_attemp/cloud_genes.txt

#97+108 merge (100):
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/core_genes.txt
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/soft_core_genes.txt
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/shell_genes.txt
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/cloud_genes.txt

#Write our gene lists to separate files
write(core_genes, file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/core_genes.txt")
write(score_genes, file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/soft_core_genes.txt")
write(shell_genes, file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/shell_genes.txt")
write(cloud_genes, file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/cloud_genes.txt")













