#!/usr/bin/env Rscript
#SBATCH --get-user-env

## Configured to run for for each chromosome
#### To run sets of 25kb sites per chromosome, we would need to:
#### minimum <- c((MAXVAL-24999))
#### sites2=c(minimum:MAXVAL)
#### for (i in sites2) {GETINFO}
#### write.table(counts[minimum:MAXVAL,........
####

library(data.table)
data=fread('all_mratios_CHROMNAME_v2.txt',header=F)
data[,V9:=paste(data[,V1],data[,V2],sep="_")]

setkey(data,V9,V8)
temp=data
data=data.table(unique(data))
names=sort(unique(data[,V8]))
sites=unique(data[,V9])

#generate counts table
counts=as.data.table(unique(data[,V9]))
setnames(counts,"V1","site")
for(i in 1:length(names))
{
print(i)
counts[,names[i] := data[.(as.factor(counts[,site]),names[i]),V5]]
}

#generate mcounts table
mcounts=as.data.table(unique(data[,V9]))
setnames(mcounts,"V1","site")
for(i in 1:length(names))
{
print(i)
mcounts[,names[i] := data[.(as.factor(mcounts[,site]),names[i]),V6]]
}

counts[is.na(counts)] <- 0
mcounts[is.na(mcounts)] <- 0
mratios<-counts
mratios[,2:(length(names)+1)]<-mcounts[,2:(length(names)+1)]/counts[,2:(length(names)+1)]

#Calculate mean mratio, sd, and avg coverage
mratios$avg=apply(mratios[,2:(length(names)+1)],1,mean,na.rm=TRUE)
mratios$sd=apply(mratios[,2:(length(names)+1)],1,sd,na.rm=TRUE)
mratios$depth<-apply(counts[,2:(length(names)+1)],1,mean,na.rm=TRUE)

# these values are sample + 1 for the rownames column
mratios$cov_anu <- apply(counts[,2:10], MARGIN=1, mean); mratios$cov_gui <- apply(counts[,11:16], MARGIN=1, mean); mratios$cov_ham <- apply(counts[,17:30], MARGIN=1, mean); mratios$cov_kin <- apply(counts[,31:34], MARGIN=1, mean); mratios$cov_mac <- apply(counts[,41:45], MARGIN=1, mean); mratios$cov_yel <- apply(counts[,35:40], MARGIN=1, mean)


subset(counts, mratios$cov_anu > 5  & mratios$cov_gui > 5 & mratios$cov_ham > 5 & mratios$cov_kin > 5 & mratios$cov_mac > 5 & mratios$cov_yel > 5) -> counts2
subset(mcounts, mratios$cov_anu > 5 & mratios$cov_gui > 5 & mratios$cov_ham > 5 & mratios$cov_kin > 5 & mratios$cov_mac > 5 & mratios$cov_yel > 5) -> mcounts2
subset(mratios, mratios$cov_anu > 5 & mratios$cov_gui > 5 & mratios$cov_ham > 5 & mratios$cov_kin > 5 & mratios$cov_mac > 5 & mratios$cov_yel > 5) -> mratios2


# save raw methylation ratios
write.table(mratios2[,2:(length(names)+1)],"methratio_CHROMNAME.txt",row.names=F,col.names=F,sep="\t")

# write methylated and total counts
write.table(mcounts2,file="mcounts_table_CHROMNAME.txt",row.names=F,sep="\t")
write.table((counts2[,2:(length(names)+1)]),file="counts_table_CHROMNAME.txt",col.names=F,row.names=F,sep="\t")

# write info
info=as.data.frame(cbind(counts2$site,mratios2$avg,mratios2$sd,mratios2$depth, mratios2$cov_anu, mratios2$cov_ham, mratios2$cov_gui, mratios2$cov_kin, mratios2$cov_yel, mratios2$cov_mac))
colnames(info) <- c("site", "avg", "sd", "depth", "anu", "ham", "gui", "kin", "yel", "mac")
# order = site ,mean, sd, depth
write.table(info,"info_CHROMNAME.txt",row.names=F,sep="\t")

# add the bed file header to the info file
library(plyr)
strsplit(as.character(info$site), "_") -> info2
df <- ldply(info2)
#paste(df$V1, df$V2, sep="_") -> df$V1
df$V3 <- df$V2
colnames(df) <- c("chrom", "site1", "site2")
cbind(df, info) -> info
write.table(info,"info_CHROMNAME.txt",row.names=F,sep="\t")

