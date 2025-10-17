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
input9 = args[9]

# STEP ARGUMENT FOR BASH SCRIPT STEP 4 

############ ADD 2Kb IN THE PREVIOUS BED FILE ---> MAKE IT READY FOR START/END OLIGOMERS SEARCHING WITH VMATCH

extended=read.table(input1, sep = "\t")

# extended=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_OUTSIDE_INTERNAL_CLEAN.bed", sep="\t")

extended_D=extended[which(extended$V17=="+"),]
extended_P=extended[which(extended$V17=="-"),]

extended_D_4_SE=extended_D
extended_P_4_SE=extended_P

extended_D$V2=extended_D$V2-as.numeric(input2)
extended_D$V3=extended_D$V3+as.numeric(input2)

extended_P$V2=extended_P$V2-as.numeric(input2)
extended_P$V3=extended_P$V3+as.numeric(input2)

# extended_D$V2=extended_D$V2-2000
# extended_D$V3=extended_D$V3+2000
# 
# extended_P$V2=extended_P$V2-2000
# extended_P$V3=extended_P$V3+2000

extended_D$V6=substr(extended_D$V6, start = 1, stop = 5)
extended_D$V11=substr(extended_D$V11, start = 16, stop = 21)

extended_P$V11=substr(extended_P$V11, start = 1, stop = 5)
extended_P$V6=substr(extended_P$V6, start = 16, stop = 21)

write.table(extended_D, input3, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(extended_P, input4, sep = "\t", row.names = F, quote= F, col.names = F)


############ TAKE 15/20/25MERS FROM THE PPT/PBS SIGNATURES IN ORDER TO FIND THE CORRECT STARTS/ENDS OF THE ELEMENTS

# head(extended_D_4_SE)
# head(extended_D)

oligomer=as.numeric(input5) - 5

extended_D_4_SE=extended_D_4_SE %>%
  mutate(V2_minus_10_equals_ES=as.numeric(as.character(V2)) - oligomer, V2_plus_5_equals_EE=as.numeric(as.character(V2)) + 5,
         V3_minus_5_equals_SS=as.numeric(as.character(V3)) - 5, V3_plus_10_equals_SE=as.numeric(as.character(V3)) + oligomer)

extended_D_4_SE$V2=extended_D_4_SE$V2-as.numeric(input2)
extended_D_4_SE$V3=extended_D_4_SE$V3+as.numeric(input2)

extended_D_4_SE=unite(extended_D_4_SE, coord, c(V2, V3), remove=FALSE, sep="-")
extended_D_4_SE=unite(extended_D_4_SE, ext_id, c(V1, coord), remove=FALSE, sep=":")

extended_D_START=extended_D_4_SE %>%
  select(V1, V3_minus_5_equals_SS, V3_plus_10_equals_SE, ext_id, V17)

extended_D_END=extended_D_4_SE %>%
  select(V1, V2_minus_10_equals_ES, V2_plus_5_equals_EE, ext_id, V17)


write.table(extended_D_START, input6, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(extended_D_END, input7, sep = "\t", row.names = F, quote= F, col.names = F)


# head(extended_P_4_SE)
# head(extended_P)

extended_P_4_SE=extended_P_4_SE %>%
  mutate(V3_minus_5_equals_EE=as.numeric(as.character(V3)) - 5, V3_plus_10_equals_ES=as.numeric(as.character(V3)) + oligomer,
         V2_minus_10_equals_SE=as.numeric(as.character(V2)) - oligomer, V2_plus_5_equals_SS=as.numeric(as.character(V2)) + 5)

extended_P_4_SE$V2=extended_P_4_SE$V2-as.numeric(input2)
extended_P_4_SE$V3=extended_P_4_SE$V3+as.numeric(input2)

extended_P_4_SE=unite(extended_P_4_SE, coord, c(V2, V3), remove=FALSE, sep="-")
extended_P_4_SE=unite(extended_P_4_SE, ext_id, c(V1, coord), remove=FALSE, sep=":")

extended_P_START=extended_P_4_SE %>%
  select(V1, V2_minus_10_equals_SE, V2_plus_5_equals_SS, ext_id, V17)

extended_P_END=extended_P_4_SE %>%
  select(V1, V3_minus_5_equals_EE, V3_plus_10_equals_ES, ext_id, V17)


write.table(extended_P_START, input8, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(extended_P_END, input9, sep = "\t", row.names = F, quote= F, col.names = F)


### FOR D

# mers=extended_D %>%
#   select(V1, V2, V3, V6, V11)
# 
# mers=unite(mers, half_id, c(V2, V3), remove=TRUE, sep="-")
# mers=unite(mers, full_id, c(V1, half_id), remove=TRUE, sep=":")

# mers$V11=substr(mers$V11, start = 12, stop = 16)
# mers$V6=substr(mers$V6, start = 1, stop = 8)

# mers=mers %>%
#   select(full_id, V11, V6)
# 
# colnames(mers)=c("TE_ID", "LTR_start", "LTR_end")
# 
# 
# start_ltr_D=as.data.frame(table(mers$LTR_start))
# start_ltr_D=start_ltr_D[order(start_ltr_D$Freq, decreasing = TRUE),]
# 
# start_ltr_D=start_ltr_D %>%
#   mutate(Oligo_ID="famous_start_D", nr=row_number())
# 
# start_ltr_D=unite(start_ltr_D, Fasta_ID, c(Oligo_ID, Freq, nr), remove=TRUE, sep="_")
# 
# write.fasta(as.list(start_ltr_D$Var1), as.list(start_ltr_D$Fasta_ID), input4)
# 
# end_ltr_D=as.data.frame(table(mers$LTR_end))
# end_ltr_D=end_ltr_D[order(end_ltr_D$Freq, decreasing = TRUE),]
# 
# end_ltr_D=end_ltr_D %>%
#   mutate(Oligo_ID="famous_end_D", nr=row_number())
# 
# end_ltr_D=unite(end_ltr_D, Fasta_ID, c(Oligo_ID, Freq, nr), remove=TRUE, sep="_")
# 
# write.fasta(as.list(end_ltr_D$Var1), as.list(end_ltr_D$Fasta_ID), input5)


### FOR P

# mers1=extended_P %>%
#   select(V1, V2, V3, V6, V11)
# 
# mers1=unite(mers1, half_id, c(V2, V3), remove=TRUE, sep="-")
# mers1=unite(mers1, full_id, c(V1, half_id), remove=TRUE, sep=":")

# mers1$V11=substr(mers1$V11, start = 1, stop = 8)
# mers1$V6=substr(mers1$V6, start = 12, stop = 16)

# mers1=mers1 %>%
#   select(full_id, V6, V11)
# 
# colnames(mers1)=c("TE_ID", "LTR_start", "LTR_end")
# 
# start_ltr_P=as.data.frame(table(mers1$LTR_start))
# start_ltr_P=start_ltr_P[order(start_ltr_P$Freq, decreasing = TRUE),]
# 
# start_ltr_P=start_ltr_P %>%
#   mutate(Oligo_ID="famous_start_P", nr=row_number())
# 
# start_ltr_P=unite(start_ltr_P, Fasta_ID, c(Oligo_ID, Freq, nr), remove=TRUE, sep="_")
# 
# write.fasta(as.list(start_ltr_P$Var1), as.list(start_ltr_P$Fasta_ID), input6)
# 
# end_ltr_P=as.data.frame(table(mers1$LTR_end))
# end_ltr_P=end_ltr_P[order(end_ltr_P$Freq, decreasing = TRUE),]
# 
# end_ltr_P=end_ltr_P %>%
#   mutate(Oligo_ID="famous_end_P", nr=row_number())
# 
# end_ltr_P=unite(end_ltr_P, Fasta_ID, c(Oligo_ID, Freq, nr), remove=TRUE, sep="_")
# 
# write.fasta(as.list(end_ltr_P$Var1), as.list(end_ltr_P$Fasta_ID), input7)

