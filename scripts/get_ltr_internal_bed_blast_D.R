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

##### D STRAND

int_tab_d=read.table(input1, sep="\t")

# int_tab_d=read.table("Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_INTERNAL.bed", sep="\t")

colnames(int_tab_d)=c("Chr", "Start_INT", "End_INT", "BLAST_ID", "PBS_Type", "PBS_Seq", "PBS_Mismatches", "PPT_Type", "PPT_Seq", "PPT_Mismatches",
                      "Strand")


full_tab_d=read.table(input2, sep="\t")

# full_tab_d=read.table("Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep="\t")

colnames(full_tab_d)=c("Chr", "Start_FL", "End_FL", "FL_ID", "Start_Seq", "End_Seq", "BLAST_ID", "Strand")

full_tab_d$BLAST_ID=sub("::.*", "", full_tab_d$BLAST_ID)

to_ltrs_d=merge(full_tab_d, int_tab_d, by.x = c("BLAST_ID"), by.y =  c("BLAST_ID"), all.x = TRUE)
to_ltrs_d=na.omit(to_ltrs_d)

to_ltrs_new_d=to_ltrs_d %>%
  mutate(End_5LTR=as.numeric(as.character(Start_INT)) - 1,
         Start_3LTR=as.numeric(as.character(End_INT)) + 1,
         LTR5_ID="5prime",
         LTR3_ID="3prime",
         INT_ID="internal")

to_ltrs_five_d=to_ltrs_new_d %>%
  select(Chr.x, Start_FL, End_5LTR, FL_ID, LTR5_ID, Strand.x)

to_ltrs_five_d=unite(to_ltrs_five_d, ID_5LTR, c(FL_ID, LTR5_ID), remove=TRUE, sep="_")

to_ltrs_three_d=to_ltrs_new_d %>%
  select(Chr.x, Start_3LTR, End_FL, FL_ID, LTR3_ID, Strand.x)

to_ltrs_three_d=unite(to_ltrs_three_d, ID_3LTR, c(FL_ID, LTR3_ID), remove=TRUE, sep="_")

to_internal_d=to_ltrs_new_d %>%
  select(Chr.x, Start_INT, End_INT, FL_ID, INT_ID, Strand.x)

to_internal_d=unite(to_internal_d, ID_INT, c(FL_ID, INT_ID), remove=TRUE, sep="_")


# write.table(to_ltrs_five_d, "Col-CEN_v1.2.fasta.D_5LTR_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_ltrs_five_d, input3, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(to_ltrs_three_d, "Col-CEN_v1.2.fasta.D_3LTR_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_ltrs_three_d, input4, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(to_internal_d, "Col-CEN_v1.2.fasta.D_INTERNAL_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(to_internal_d, input5, sep = "\t", row.names = F, quote= F, col.names = F)