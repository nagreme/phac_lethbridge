library("dplyr")
library("magrittr")

# setwd("~/Desktop/scripts")
df1 <- read.table(file="sample_file1.csv", header=T,  sep=",")
df2 <- read.csv(file="sample_file2.csv", header=T,  sep=",")

#this one is matching by columns... I'll have to transpose it I guess NOPE
#df1_and_2 <- full_join(df1,df2, by=NULL, copy=F)

#test renaming genes with suffixes for the interfile merge 
df1[,1] <- paste0(df1[,1],".1",collapse=NULL)
df2[,1] <- paste0(df2[,1],".2",collapse=NULL)


df1_tbl <- df1 %>% tbl_df()
df2_tbl <- df2 %>% tbl_df()

df1_and_2 <- full_join(df1_tbl,df2_tbl, copy=F) 
df1_and_2[] <- lapply(df1_and_2, function(x){x[is.na(x)] <- 0; x})


#only matches same row names, doesn't work if row.names = 1)
merge(df1,df2,all = T) 

merge(df1,df2, by.x=NULL,by.y=NULL)



#intrafile merge
df2["new",] <- df2["fourth",] | df2["fifth",]
#df2 <- df2[-"new",] #doesn't work?

subset(df2,rownames(df2)!="new") #this is cleaner for more complex logicals

df2[rownames(df2)!="new",]#simple way for what I'm doing





