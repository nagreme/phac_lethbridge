#interfile merging script

# Nadège Pulgar-Vidal
# 28 fév 2017

# Read in both .Rtabs do a full join then merge the rows using the parsed .clstr file (.csv) 
# The row merging will be the same as the row merging in the intrfile merge script

library("dplyr")
library("magrittr")


#97+108
#.Rtab 97: /home/nadege/Desktop/metagenome/cd-hit-est_attempt/gene_presence_absence_97.Rtab
#.Rtab 108: /home/nadege/Desktop/metagenome/cd-hit-est_attempt/gene_presence_absence_108.Rtab
#.csv: /home/nadege/Desktop/metagenome/cd-hit-est_attempt/pan_gen_97_108_clstrd.csv

file1 <- "/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/gene_presence_absence_97.Rtab"
file2 <- "/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/gene_presence_absence_108.Rtab"
clstr_file <- "/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/pan_gen_97_108_clstrd.csv"

#read the .Rtab  files
df1 <- read.table(file=file1, header=T,  sep=",")#, row.names = 1)
df2 <- read.table(file=file2, header=T,  sep=",")#, row.names = 1)

#Add suffixes to avoid same names merging together
df1[,1] <- paste0(df1[,1],".1",collapse=NULL)
df2[,1] <- paste0(df2[,1],".2",collapse=NULL)

#get the cluster data ready
clusts <-  read.csv(file=clstr_file, header=F, sep="\t", col.names = paste0("V",seq_len(3)), fill = T, stringsAsFactors = F, row.names = 1) 
clusts[] <- lapply(clusts, function(x) {x[x==""] <- NA; x})
row_lens <- apply(clusts,1,function(x) length(na.omit(x)))

#Do a full join on the .Rtab files
#the full join needs the tbl objects to work right
#df1_tbl <- df1 %>% tbl_df()
#df2_tbl <- df2 %>% tbl_df()

#Unfortunately this will merge rows with the same name that should not be merged
#I need a join/merge without matching: just add rows and columns and prolong the existing ones, no matching
#df1_and_2 <- full_join(df1,df2, copy=F) 
df1_and_2 <- merge(df1,df2,all=T)
df1_and_2[] <- lapply(df1_and_2, function(x){x[is.na(x)] <- 0; x})

#set row names to genes then delete those from the table
rownames(df1_and_2) <- df1_and_2[,1] 
df1_and_2 <- df1_and_2[,-1]


#then merge rows (intrafile merge)

#Go through each cluster
for (i in (1:nrow(clusts)))
{
  #If there is more than one member in the cluster (+row/cluster name)
  if (row_lens[i] > 1)
  {
    gene1 <- clusts[i,1]
    gene2 <- clusts[i,2]
    #We merge it with the representative (there should be at most 2, one from each file)
    df1_and_2[gene1,] <- df1_and_2[gene1,] | df1_and_2[gene2,]
      
    #then remove that row
    df1_and_2 <- subset(df1_and_2,rownames(df1_and_2)!=gene2)
  }#if merge needed
}#outer for (each cluster)

#get rid of the suffixes (I'm not doing this anymore I think)
#rownames(df1_and_2) <- sub("\.[01]$","",rownames(df1_and_2))
#rownames(df1_and_2) <- substr(rownames(df1_and_2),1,nchar(rownames(df1_and_2))-2)

write.csv(df1_and_2,file="/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/gene_pres_abs_97_108.Rtab")


