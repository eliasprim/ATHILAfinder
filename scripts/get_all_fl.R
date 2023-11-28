rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))


# setwd('~/Desktop')


args=commandArgs(trailingOnly = TRUE)

input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]


nonover=read.table(input1, sep="\t")

# nonover=read.table("Col-CEN_v1.2.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING.bed", sep="\t")

spends=read.table(input2, sep="\t")

# spends=read.table("Col-CEN_v1.2.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep="\t")

result=nonover[!(nonover$V4 %in% spends$V7),]

result_new=result %>%
  select(-V5, -V6) %>%
  mutate(start="non_specific_start", end="non_specific_end", origin=V4) %>%
  select(V1, V2, V3, V4, start, end, origin, V7)

colnames(result_new)=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")

allfl=rbind(spends, result_new)

write.table(allfl, input3, sep = "\t", row.names = F, quote= F, col.names = F)

allfl_d=allfl[which(allfl$V8=="+"),]

write.table(allfl_d, input4, sep = "\t", row.names = F, quote= F, col.names = F)

allfl_p=allfl[which(allfl$V8=="-"),]

write.table(allfl_p, input5, sep = "\t", row.names = F, quote= F, col.names = F)


