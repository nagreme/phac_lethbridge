#intrafile merging script

# Nadège Pulgar-Vidal
# 27 fév 2017

#Merge rows within a .Rtab file according to the parsed .clstr file (.csv) using an OR function

library("dplyr")
library("magrittr")


#97 
#.Ratb: /home/nadege/Desktop/roary/old_97_refined_meta/gene_presence_absence.Rtab
#.csv: /home/nadege/Desktop/metagenome/cd-hit-est_attempt/pan_gen_97_clstrd.csv (99)
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/pan_gen_97_clstrd.csv (100)

#108 
#.Rtab: /home/nadege/Desktop/roary/corrected_prokka_db_full_set_exclude_CP017855/gene_presence_absence.Rtab
#.csv: /home/nadege/Desktop/metagenome/cd-hit-est_attempt/pan_gen_108_clstrd.csv (99)
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/pan_gen_108_clstrd.csv (100)

#108+97 (100) control:
#.Rtab: /home/nadege/Desktop/roary/97_108_refined_meta/gene_presence_absence.Rtab
#.csv: /home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/control/pan_gen_100_clstrd.csv

#108+97 (99) control:
#.Rtab: /home/nadege/Desktop/roary/97_108_refined_meta/gene_presence_absence.Rtab
#.csv: /home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_99/control/pan_gen_100_clstrd.csv

#5695 (100):
#.Rtab: /home/nadege/Desktop/roary/big_set/gene_presence_absence.Rtab
#.csv; "/home/nadege/Desktop/mist/big_set/pan_gen_clstrd.csv"

rtab_file <- "/home/nadege/Desktop/roary/big_set/gene_presence_absence.Rtab"
clstr_file <- "/home/nadege/Desktop/mist/big_set/pan_gen_clstrd.csv"

df <- read.table(file=rtab_file, header=T,  sep="\t", row.names = 1)
clusts <-  read.csv(file=clstr_file, header=F, sep="\t", col.names = paste0("V",seq_len(18)), fill = T, stringsAsFactors = F, row.names = 1) 
clusts[] <- lapply(clusts, function(x) {x[x==""] <- NA; x})

#I can't get this to work how I want it to on one row so I'll just use this to precompute it and store it
row_lens <- apply(clusts,1,function(x) length(na.omit(x)))

#Go through each cluster
for (i in (1:nrow(clusts)))
{
  #If there is more than one member in the cluster (+row/cluster name)
  if (row_lens[i] > 1)
  {
    #We merge them with the representative
    for (j in (2:row_lens[i]))
    {
      #Assumes .1 and .2 were already added to the gene annotations in/before the clstr file but not in .Rtab
      #If not remove -2 after nchar()
      gene1 <- substr(clusts[i,1],1,nchar(clusts[i,1])) 
      gene2 <- substr(clusts[i,j],1,nchar(clusts[i,j]))
      df[gene1,] <- df[gene1,] | df[gene2,]

      #then remove that row
      df <- subset(df,rownames(df)!=gene2)
    } #inner for (members to merge)
  }#if merge needed
}#outer for (each cluster)

write.csv(df1,file="/home/nadege/Desktop/mist/big_set/gene_presence_absence_100.Rtab")

  
