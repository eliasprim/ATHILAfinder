rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# setwd('~/Desktop/athilafinder_in_progress_23092021/')

args=commandArgs(trailingOnly = TRUE)

input1 = args[1]
input2 = args[2]

bed_file=read.table(input1, sep="\t")

mychr=unique(bed_file$V1)

for (chr in mychr)
{
  tmp1=bed_file[!duplicated(bed_file[,c('V2')]),]
  tmp2=tmp1[!duplicated(tmp1[,c('V3')]),]
}

write.table(tmp2, input2, sep = "\t", row.names = F, quote= F, col.names = F)
