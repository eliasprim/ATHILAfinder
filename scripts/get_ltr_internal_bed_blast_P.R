rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# setwd('~/Desktop')

args=commandArgs(trailingOnly = TRUE)

input1 = args[1]
input2 = args[2]
input3 = args[3]
input4 = args[4]
input5 = args[5]


# # STEP ARGUMENT FOR BASH SCRIPT STEP 6


############### TAKE LTRs AND INTERNAL REGIONS

##### P STRAND

int_tab_p=read.table(input1, sep="\t")

# int_tab_p=read.table("Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_INTERNAL.bed", sep="\t")

colnames(int_tab_p)=c("Chr", "Start_INT", "End_INT", "BLAST_ID", "PPT_Type", "PPT_Seq", "PPT_Mismatches", "PBS_Type", "PBS_Seq", "PBS_Mismatches",
                      "Strand")

full_tab_p=read.table(input2, sep="\t")

# full_tab_p=read.table("Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep="\t")

colnames(full_tab_p)=c("Chr", "Start_FL", "End_FL", "FL_ID", "Start_Seq", "End_Seq", "BLAST_ID", "Strand")


to_ltrs_p=merge(full_tab_p, int_tab_p, by.x = c("BLAST_ID"), by.y =  c("BLAST_ID"), all.x = TRUE)
to_ltrs_p=na.omit(to_ltrs_p)

to_ltrs_new_p=to_ltrs_p %>%
  mutate(End_5LTR=as.numeric(as.character(End_INT)) + 1,
         Start_3LTR=as.numeric(as.character(Start_INT)) - 1,
         LTR5_ID="5prime",
         LTR3_ID="3prime",
         INT_ID="internal")

to_ltrs_five_p=to_ltrs_new_p %>%
  select(Chr.x, End_5LTR, End_FL, FL_ID, LTR5_ID, Strand.x)

to_ltrs_five_p=unite(to_ltrs_five_p, ID_5LTR, c(FL_ID, LTR5_ID), remove=TRUE, sep="_")

to_ltrs_three_p=to_ltrs_new_p %>%
  select(Chr.x, Start_FL, Start_3LTR, FL_ID, LTR3_ID, Strand.x)

to_ltrs_three_p=unite(to_ltrs_three_p, ID_3LTR, c(FL_ID, LTR3_ID), remove=TRUE, sep="_")

to_internal_p=to_ltrs_new_p %>%
  select(Chr.x, Start_INT, End_INT, FL_ID, INT_ID, Strand.x)

to_internal_p=unite(to_internal_p, ID_INT, c(FL_ID, INT_ID), remove=TRUE, sep="_")


# write.table(to_ltrs_five_p, "Col-CEN_v1.2.fasta.P_5LTR_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_ltrs_five_p, input3, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(to_ltrs_three_p, "Col-CEN_v1.2.fasta.P_3LTR_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_ltrs_three_p, input4, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(to_internal_p, "Col-CEN_v1.2.fasta.P_INTERNAL_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_internal_p, input5, sep = "\t", row.names = F, quote= F, col.names = F)