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
input10 = args[10]
# input11 = args[11]
# input12 = args[12]
# input13 = args[13]
# input14 = args[14]
# input15 = args[15]
# input16 = args[16]
# input17 = args[17]
# input18 = args[18]
# input19 = args[19]


# # STEP ARGUMENT FOR BASH SCRIPT STEP 6


############### TAKE LTRs AND INTERNAL REGIONS

##### D STRAND

ext_tab=read.table(input1, sep="\t")

# ext_tab=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_D_EXTENDED_INTERNAL_CLEAN.bed", sep="\t")

ext_tab=ext_tab %>%
  select(V1, V2, V3, V4, V9, V14)

ext_tab=unite(ext_tab, half_id, c(V2, V3), remove=TRUE, sep="-")
ext_tab=unite(ext_tab, full_id, c(V1, half_id), remove=TRUE, sep=":")

full_tab=read.table(input2, sep="\t")

# full_tab=read.table("t2t-col.20210610.fasta.F2B.D_FULLLENGTH.bed", sep="\t")


whole_tab=merge(full_tab, ext_tab, by.x = c("V4"), by.y = c("full_id"), all.x = TRUE)

### change for P
colnames(whole_tab)=c("Ext_Inter_ID", "Chr", "Full_Start", "Full_End", "Start_Mer",
                      "End_Mer", "Strand", "PBS_ID", "PPT_ID", "Distance_PBS_PPT")

whole_tab=whole_tab %>%
  separate(PBS_ID, c("arxh_pbs", "telos_pbs"), sep = "\\.", remove=FALSE)

whole_tab$arxh_pbs=gsub(".*_", "", whole_tab$arxh_pbs)
whole_tab$telos_pbs=gsub("_[^_]+$", "", whole_tab$telos_pbs)
whole_tab$telos_pbs=gsub("_D", "", whole_tab$telos_pbs)

whole_tab=whole_tab %>%
  separate(PPT_ID, c("arxh_ppt", "telos_ppt"), sep="\\.", remove=FALSE)

whole_tab$arxh_ppt=gsub(".*_", "", whole_tab$arxh_ppt)
whole_tab$telos_ppt=gsub("_[^_]+$", "", whole_tab$telos_ppt)
whole_tab$telos_ppt=gsub("_D", "", whole_tab$telos_ppt)


### change for P

# prosthetoume 8 sthn arxh tou PBS kai afairoume 5 apo to telos tou PPT gia to D strand

whole_tab$arxh_pbs=as.numeric(as.character(whole_tab$arxh_pbs)) + 5
whole_tab$telos_ppt=as.numeric(as.character(whole_tab$telos_ppt)) - 5

summary_tab=whole_tab %>%
  select(Ext_Inter_ID, Chr, Full_Start, Full_End, Start_Mer, End_Mer, arxh_pbs, telos_ppt, Distance_PBS_PPT, Strand)

### change for P

LTR5_coords=summary_tab %>%
  select(Chr, Full_Start, arxh_pbs, Start_Mer, Ext_Inter_ID, Strand)

LTR3_coords=summary_tab %>%
  select(Chr, telos_ppt, Full_End, End_Mer, Ext_Inter_ID, Strand)

INT_coords=summary_tab %>%
  select(Chr, arxh_pbs, telos_ppt, Ext_Inter_ID, Strand)

# LTR5_coords_smt=LTR5_coords %>%
#   mutate("start_plus"=as.numeric(as.character(Full_Start)) + as.numeric(input3))
# 
# LTR5_coords_smt=LTR5_coords_smt %>%
#   select(Chr, Full_Start, start_plus, Start_5mer, Ext_Inter_ID, Strand)

# LTR5_coords_smt2=LTR5_coords %>%
#   mutate("start_minus"=as.numeric(as.character(Full_Start)) - as.numeric(input3), 
#          "incl_st"=as.numeric(as.character(Full_Start)) + 5)

# LTR5_coords_smt2=LTR5_coords_smt2 %>%
#   select(Chr, start_minus, incl_st, Start_5mer, Ext_Inter_ID, Strand)

# LTR3_coords_smt=LTR3_coords %>%
#   mutate("end_minus"=as.numeric(as.character(Full_End)) - as.numeric(input3))
# 
# LTR3_coords_smt=LTR3_coords_smt %>%
#   select(Chr, end_minus, Full_End, End_5mer, Ext_Inter_ID, Strand)

# LTR3_coords_smt2=LTR3_coords %>%
#   mutate("end_plus"=as.numeric(as.character(Full_End)) + as.numeric(input3),
#          "incl_end"=as.numeric(as.character(Full_End)) - 5)

# LTR3_coords_smt2=LTR3_coords_smt2 %>%
#   select(Chr, incl_end, end_plus, End_5mer, Ext_Inter_ID, Strand)


write.table(LTR5_coords, input3, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(LTR3_coords, input4, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(INT_coords, input5, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR5_coords_smt, input7, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR5_coords_smt2, input8, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR3_coords_smt, input8, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR3_coords_smt2, input10, sep = "\t", row.names = F, quote= F, col.names = F)






##### P STRAND

ext_tab=read.table(input6, sep="\t")

# ext_tab=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_P_EXTENDED_INTERNAL_CLEAN.bed", sep="\t")

ext_tab=ext_tab %>%
  select(V1, V2, V3, V4, V9, V14)


ext_tab=unite(ext_tab, half_id, c(V2, V3), remove=TRUE, sep="-")
ext_tab=unite(ext_tab, full_id, c(V1, half_id), remove=TRUE, sep=":")
ext_tab$full_id=paste0(ext_tab$full_id, "_P_element")
ext_tab$full_id=paste("rc", ext_tab$full_id, sep="_")


full_tab=read.table(input7, sep="\t")

# full_tab=read.table("t2t-col.20210610.fasta.F2B.P_FULLLENGTH.bed", sep="\t")

# full_tab$V4=paste0(full_tab$V4, "_P_element")
# full_tab$V4=paste("rc", full_tab$V4, sep="_")


whole_tab=merge(full_tab, ext_tab, by.x = c("V4"), by.y = c("full_id"), all.x = TRUE)


### change for P
colnames(whole_tab)=c("Ext_Inter_ID", "Chr", "Full_Start", "Full_End", "Start_Mer",
                      "End_Mer", "Strand", "PPT_ID", "PBS_ID", "Distance_PPT_PBS")

whole_tab=whole_tab %>%
  separate(PBS_ID, c("arxh_pbs", "telos_pbs"), sep = "\\.", remove=FALSE)

whole_tab$arxh_pbs=gsub(".*_", "", whole_tab$arxh_pbs)
whole_tab$telos_pbs=gsub("_[^_]+$", "", whole_tab$telos_pbs)
whole_tab$telos_pbs=gsub("_P", "", whole_tab$telos_pbs)

whole_tab=whole_tab %>%
  separate(PPT_ID, c("arxh_ppt", "telos_ppt"), sep="\\.", remove=FALSE)

whole_tab$arxh_ppt=gsub(".*_", "", whole_tab$arxh_ppt)
whole_tab$telos_ppt=gsub("_[^_]+$", "", whole_tab$telos_ppt)
whole_tab$telos_ppt=gsub("_P", "", whole_tab$telos_ppt)


### change for P

# afairoume 8 apo to telos tou PBS kai prosthetoume 5 sto telos tou PPT gia to P strand

whole_tab$telos_pbs=as.numeric(as.character(whole_tab$telos_pbs)) - 5
whole_tab$arxh_ppt=as.numeric(as.character(whole_tab$arxh_ppt)) + 5

summary_tab=whole_tab %>%
  select(Ext_Inter_ID, Chr, Full_Start, Full_End, Start_Mer, End_Mer, arxh_ppt, telos_pbs, Distance_PPT_PBS, Strand)

### change for P

LTR5_coords=summary_tab %>%
  select(Chr, telos_pbs, Full_End, Start_Mer, Ext_Inter_ID, Strand)

LTR3_coords=summary_tab %>%
  select(Chr, Full_Start, arxh_ppt, End_Mer, Ext_Inter_ID, Strand)

INT_coords=summary_tab %>%
  select(Chr, arxh_ppt, telos_pbs, Ext_Inter_ID, Strand)

# LTR5_coords_smt=LTR5_coords %>%
#   mutate("start_minus"=as.numeric(as.character(Full_End)) - as.numeric(input3))
# 
# LTR5_coords_smt=LTR5_coords_smt %>%
#   select(Chr, start_minus, Full_End, Start_5mer, Ext_Inter_ID, Strand)

# LTR5_coords_smt2=LTR5_coords %>%
#   mutate("start_plus"=as.numeric(as.character(Full_End)) + as.numeric(input3),
#          "incl_st"=as.numeric(as.character(Full_End)) - 5)

# LTR5_coords_smt2=LTR5_coords_smt2 %>%
#   select(Chr, incl_st, start_plus, Start_5mer, Ext_Inter_ID, Strand)

# LTR3_coords_smt=LTR3_coords %>%
#   mutate("end_plus"=as.numeric(as.character(Full_Start)) + as.numeric(input3))
# 
# LTR3_coords_smt=LTR3_coords_smt %>%
#   select(Chr, Full_Start, end_plus, End_5mer, Ext_Inter_ID, Strand)

# LTR3_coords_smt2=LTR3_coords %>%
#   mutate("end_minus"=as.numeric(as.character(Full_Start)) - as.numeric(input3),
#          "incl_end"=as.numeric(as.character(Full_Start)) + 5)

# LTR3_coords_smt2=LTR3_coords_smt2 %>%
#   select(Chr, end_minus, incl_end, End_5mer, Ext_Inter_ID, Strand)




write.table(LTR5_coords, input8, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(LTR3_coords, input9, sep = "\t", row.names = F, quote= F, col.names = F)
write.table(INT_coords, input10, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR5_coords_smt, input14, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR5_coords_smt2, input17, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR3_coords_smt, input15, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(LTR3_coords_smt, input19, sep = "\t", row.names = F, quote= F, col.names = F)




