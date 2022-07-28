rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# setwd('~/Desktop/athilafinder_in_progress_23092021/')

args=commandArgs(trailingOnly = TRUE)

input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]
input6 = args[6]
input7 = args[7]
input8 = args[8]

# STEP ARGUMENT FOR BASH SCRIPT STEP 3

############ CHECK SECOND WINDOW ANALYSIS OUTPUT

outside_internal=read.table(input1, sep="\t")

# outside_internal=read.table("t2t-col.20210610.fasta.F2B.ALL_PBSPPT_DP_OUTSIDE_INTERNAL.bed", sep="\t")

# fghjkl=read.table("t2t-col.20210610.fasta.F2B.ALL_PBSPPT_DP_WINDOW2_SORTED.bed", sep="\t")
# dim(fghjkl)
# 
# unique(fghjkl$V4)
# unique(fghjkl$V9)
# 
# rty1=fghjkl%>%
#   select(V4)
# 
# colnames(rty1)[1]="oligomer"
# 
# rty2=fghjkl%>%
#   select(V9)
# 
# colnames(rty2)[1]="oligomer"
# 
# rty3=rbind(rty1, rty2)
# dim(rty3)
# unique(rty3$oligomer)

# BEDTOOLS window command also finds parts (PBS and PPT) that are parts of these internal elements,
# so we have to calculate the following distances and keep the non-zero ones

# outside_internal=outside_internal %>%
#   mutate(middle_of_oligomer=(as.numeric(as.character(V17)) + as.numeric(as.character(V18)))/2)

outside_internal=outside_internal %>%
  mutate(before_start=as.numeric(as.character(V2)) - as.numeric(as.character(V17)),
         after_end=as.numeric(as.character(V3)) - as.numeric(as.character(V18)))

without_zero=outside_internal[which(outside_internal$before_start!=0 & outside_internal$after_end!=0),]

new_without_zero=without_zero %>%
  select(V1:V15)

write.table(new_without_zero, input2, sep = "\t", row.names = F, quote= F, col.names = F)

# GET the internal elements without a PBS or a PPT close to them [2,5kb] (even if other internal element)   

input_window=read.table(input3, sep="\t")

# input_window=read.table("t2t-col.20210610.fasta.F2B.ALL_PBSPPT_DP_WINDOW2.bed", sep="\t")

rrr=anti_join(input_window, new_without_zero)

write.table(rrr, input4, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(rrr, "t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_OUTSIDE_INTERNAL_CLEAN.bed", sep = "\t", row.names = F, quote= F, col.names = F)

# separate elements based on the strand in order to getfasta

strd=rrr[which(rrr$V15=="+"),]

write.table(strd, input5, sep = "\t", row.names = F, quote= F, col.names = F)


strp=rrr[which(rrr$V15=="-"),]

write.table(strp, input6, sep = "\t", row.names = F, quote= F, col.names = F)

# Proximal elements

try_smt2=intersect(without_zero$V4, without_zero$V19)
try_smt2=as.data.frame(try_smt2)
colnames(try_smt2)=c("try_something")

try_smt3=intersect(without_zero$V9, without_zero$V19)
try_smt3=as.data.frame(try_smt3)
colnames(try_smt3)=c("try_something")

overall_try_smt=rbind(try_smt2, try_smt3)

proximal_elements=NULL;
for (i in overall_try_smt$try_something)
{
  tmp=without_zero[which(without_zero$V19==i),]
  proximal_elements=rbind(proximal_elements, tmp)
}

write.table(proximal_elements, input7, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(proximal_elements, "t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_OUTSIDE_INTERNAL_PROXIMAL_ELEMENTS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

# Orphan oligomers 

try_smt4=setdiff(without_zero$V19, without_zero$V4)

try_smt5=setdiff(try_smt4, without_zero$V9)
try_smt5=as.data.frame(try_smt5)
colnames(try_smt5)=c("something_try")

orphan_oligomers=NULL;
for (i in try_smt5$something_try)
{
  tmp=without_zero[which(without_zero$V19==i),]
  orphan_oligomers=rbind(orphan_oligomers, tmp)
}

write.table(orphan_oligomers, input8, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(orphan_oligomers, "t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_OUTSIDE_INTERNAL_ORPHAN_OLIGOMERS.bed", sep = "\t", row.names = F, quote= F, col.names = F)