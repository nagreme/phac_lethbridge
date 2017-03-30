#extract the accessory from the huge .Rtab file
#and subset genes?

# Nadège Pulgar-Vidal
# 1 mars 2017

library("magrittr")

#/home/nadege/Desktop/roary/big_set/gene_presence_absence.Rtab
#/home/nadege/Desktop/mist/big_set/gene_presence_absence_100.Rtab

filename <- "/home/nadege/Desktop/mist/big_set/gene_presence_absence_100.Rtab"

df <- read.table(file=filename, sep=",", header=T, row.names = 1)

dft <- as.data.frame(t(df))

#Filter the genomes we want:

#read in the file Dillon gave me with the genomes we want as row names (by heights)
#/home/nadege/Desktop/acc_clustering/goeburst_clusters.tsv

heights <- read.csv(file="/home/nadege/Desktop/acc_clustering/goeburst_clusters.tsv", sep="\t", header=T, row.names = 1)

#retrieve rownames from that
wanted_genomes <- rownames(heights)

#filter new_df => if any(curr_genome_name == rownames)
filtered_df <- dft[which(rownames(dft) %in% wanted_genomes),]

#we can use mean in this case because all the values are 0 or 1
carriage <- sapply(filtered_df,mean)
#carriage <- sapply(dft,mean)

acc_df <- filtered_df[,which(carriage < 0.999 & carriage > 0)] %>% t() 
#acc_df <- dft[,which(carriage < 0.999 & carriage > 0)] %>% t() #transpose for jaccard distance

write.csv(acc_df,file="/home/nadege/Desktop/acc_clustering/acc_presence_absence_big_set.Rtab")

num_isolates <- lapply(acc_df,sum)











