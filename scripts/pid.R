rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(pwalign))
suppressMessages(library(seqinr))
suppressMessages(library(dplyr))

# setwd('~/Desktop')


args=commandArgs(trailingOnly = TRUE)

input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]
input6 = args[6]
input7 = args[7]
input8 = args[8]
input9 = args[9]
input10 = args[10]



ltr5=read.table(input1, sep="\t")

# ltr5=read.table("t2t-col.20210610.fasta.F2B.DP_5LTR_BLAST_TOGETHER.fasta_table.txt", sep="\t")

ltr5$V1=gsub("_5prime", "", ltr5$V1)

colnames(ltr5)=c("FL_ID", "LTR5_SEQ")


ltr3=read.table(input2, sep="\t")

# ltr3=read.table("t2t-col.20210610.fasta.F2B.DP_3LTR_BLAST_TOGETHER.fasta_table.txt", sep="\t")

ltr3$V1=gsub("_3prime", "", ltr3$V1)

colnames(ltr3)=c("FL_ID", "LTR3_SEQ")


ltrs=merge(ltr5, ltr3, by.x = c("FL_ID"), by.y = c("FL_ID"), all.x = TRUE)


ltrs_pid1=ltrs %>%
  mutate(perc_ident1=pid(pairwiseAlignment(LTR5_SEQ, LTR3_SEQ, type="global"), type = "PID1"))

ltrs_pid2=ltrs %>%
  mutate(perc_ident2=pid(pairwiseAlignment(LTR5_SEQ, LTR3_SEQ, type="global"), type = "PID2"))

ltrs_pid3=ltrs %>%
  mutate(perc_ident3=pid(pairwiseAlignment(LTR5_SEQ, LTR3_SEQ, type="global"), type = "PID3"))

ltrs_pid4=ltrs %>%
  mutate(perc_ident4=pid(pairwiseAlignment(LTR5_SEQ, LTR3_SEQ, type="global"), type = "PID4"))


ltr53_pid1=ltrs_pid1 %>%
  select(FL_ID, perc_ident1)

ltr53_pid1$perc_ident1=round(ltr53_pid1$perc_ident1, digits = 5)

ltr53_pid2=ltrs_pid2 %>%
  select(FL_ID, perc_ident2)

ltr53_pid2$perc_ident2=round(ltr53_pid2$perc_ident2, digits = 5)

ltr53_pid3=ltrs_pid3 %>%
  select(FL_ID, perc_ident3)

ltr53_pid3$perc_ident3=round(ltr53_pid3$perc_ident3, digits = 5)

ltr53_pid4=ltrs_pid4 %>%
  select(FL_ID, perc_ident4)

ltr53_pid4$perc_ident4=round(ltr53_pid4$perc_ident4, digits = 5)



write.table(ltr53_pid1, input3, sep = "\t", row.names = F, quote= F, col.names = T)

write.table(ltr53_pid2, input4, sep = "\t", row.names = F, quote= F, col.names = T)

write.table(ltr53_pid3, input5, sep = "\t", row.names = F, quote= F, col.names = T)

write.table(ltr53_pid4, input6, sep = "\t", row.names = F, quote= F, col.names = T)


# write.table(ltr53_pid1, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID1_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)

# write.table(ltr53_pid2, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID2_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)

# write.table(ltr53_pid3, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID3_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)

# write.table(ltr53_pid4, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID4_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)


png(as.character(input7))
hist(ltr53_pid1$perc_ident1, breaks=100, xlim=c(0,100), main='Percentage Identity (PID1) of Full-length Elements',
     xlab="Percentage Identity", ylab="Frequency", col=1)
abline(v=mean(ltr53_pid1$perc_ident1), col='red', lty=2, lwd=2)
dev.off()

png(as.character(input8))
hist(ltr53_pid2$perc_ident2, breaks=100, xlim=c(0,100), main='Percentage Identity (PID2) of Full-length Elements',
     xlab="Percentage Identity", ylab="Frequency", col=1)
abline(v=mean(ltr53_pid2$perc_ident2), col='red', lty=2, lwd=2)
dev.off()

png(as.character(input9))
hist(ltr53_pid3$perc_ident3, breaks=100, xlim=c(0,100), main='Percentage Identity (PID3) of Full-length Elements',
     xlab="Percentage Identity", ylab="Frequency", col=1)
abline(v=mean(ltr53_pid3$perc_ident3), col='red', lty=2, lwd=2)
dev.off()

png(as.character(input10))
hist(ltr53_pid4$perc_ident4, breaks=100, xlim=c(0,100), main='Percentage Identity (PID4) of Full-length Elements',
     xlab="Percentage Identity", ylab="Frequency", col=1)
abline(v=mean(ltr53_pid4$perc_ident4), col='red', lty=2, lwd=2)
dev.off()



