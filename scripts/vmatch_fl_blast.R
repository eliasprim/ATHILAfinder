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
input6 = args[6]
input7 = args[7]
input8 = args[8]
input9 = args[9]
input10 = args[10]
input11 = args[11]
input12 = args[12]
input13 = args[13]
input14 = args[14]
input15 = args[15]
input16 = args[16]


############ LENGTH TABLE FOR BLAST ELEMENTS

length_table=read.table(input1, sep="\t")

# length_table=read.table("Evs-0.ragtag_scaffolds.fa.DP_FULLLENGTH_BLAST_NONOVERLAPPING_LENGTH.txt", sep="\t")

############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PBS EXCEPT 4 and 8

vmatch_tab_a=read.table(input2, fill = T)

# vmatch_tab_a=read.table("Evs-0.ragtag_scaffolds.fa.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PBS_EXCEPT4n8.vmatch", fill=T)

flt1_a=vmatch_tab_a[which(vmatch_tab_a$V1!="Query:" & vmatch_tab_a$V1!="!" & vmatch_tab_a$V1!="!!"),]

flt2_start_a=flt1_a[which(flt1_a$V1=="21" | flt1_a$V1=="22"),]

colnames(flt2_start_a)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_end_a=flt1_a[which(flt1_a$V1=="Sbjct:"),]

colnames(flt3_end_a)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt4_a=cbind(flt2_start_a, flt3_end_a)

final_flt_a=flt4_a %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

final_flt_a$Mismatches=gsub("-", "", final_flt_a$Mismatches)

colnames(final_flt_a)[8]="Length"



############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PBS 4

vmatch_tab_b=read.table(input3, fill = T)

# vmatch_tab_b=read.table("Evs-0.ragtag_scaffolds.fa.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PBS_ATHILA4.vmatch", fill=T)

flt1_b=vmatch_tab_b[which(vmatch_tab_b$V1!="Query:" & vmatch_tab_b$V1!="!" & vmatch_tab_b$V1!="!!"),]

flt2_start_b=flt1_b[which(flt1_b$V1=="21" | flt1_b$V1=="22"),]

colnames(flt2_start_b)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_end_b=flt1_b[which(flt1_b$V1=="Sbjct:"),]

colnames(flt3_end_b)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt4_b=cbind(flt2_start_b, flt3_end_b)

final_flt_b=flt4_b %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

final_flt_b$Mismatches=gsub("-", "", final_flt_b$Mismatches)

colnames(final_flt_b)[8]="Length"



############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PBS 8

vmatch_tab_c=read.table(input4, fill = T)

# vmatch_tab_c=read.table("Evs-0.ragtag_scaffolds.fa.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PBS_ATHILA8.vmatch", fill=T)

flt1_c=vmatch_tab_c[which(vmatch_tab_c$V1!="Query:" & vmatch_tab_c$V1!="!" & vmatch_tab_c$V1!="!!"),]

flt2_start_c=flt1_c[which(flt1_c$V1=="21" | flt1_c$V1=="22"),]

colnames(flt2_start_c)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_end_c=flt1_c[which(flt1_c$V1=="Sbjct:"),]

colnames(flt3_end_c)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt4_c=cbind(flt2_start_c, flt3_end_c)

final_flt_c=flt4_c %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

final_flt_c$Mismatches=gsub("-", "", final_flt_c$Mismatches)

colnames(final_flt_c)[8]="Length"


pbs_table=rbind(final_flt_a, final_flt_b, final_flt_c)

pbs_table$Start=as.numeric(pbs_table$Start)
pbs_table$End=as.numeric(pbs_table$End)
pbs_table$Mismatches=as.numeric(pbs_table$Mismatches)

# ltr_len1=500
# ltr_len2=2500

ltr_len1=as.numeric(input5)
ltr_len2=as.numeric(input6)

pbs_length=pbs_table[which(pbs_table$Start>=ltr_len1 & pbs_table$Start<=ltr_len2),]

pbs_length_new=pbs_length[order(pbs_length$Mismatches),]

pbs_length_new=pbs_length_new %>%
  distinct(Chr, Start, .keep_all = TRUE)

pbs_length_new=as.data.table(pbs_length_new)

pbs_final=pbs_length_new[pbs_length_new[, .I[which.min(Mismatches)], by=Chr]$V1]

pbs_final=as.data.frame(pbs_final)


############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PPT

vmatch_tab_d=read.table(input7, fill = T)

# vmatch_tab_d=read.table("Evs-0.ragtag_scaffolds.fa.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PPT.vmatch", fill=T)

flt1_d=vmatch_tab_d[which(vmatch_tab_d$V1!="Query:" & vmatch_tab_d$V1!="!" & vmatch_tab_d$V1!="!!"),]

flt2_start_d=flt1_d[which(flt1_d$V1=="21" | flt1_d$V1=="22"),]

colnames(flt2_start_d)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_end_d=flt1_d[which(flt1_d$V1=="Sbjct:"),]

colnames(flt3_end_d)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt4_d=cbind(flt2_start_d, flt3_end_d)

final_flt_d=flt4_d %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

final_flt_d$Mismatches=gsub("-", "", final_flt_d$Mismatches)

colnames(final_flt_d)[8]="Length"


ppt_and_length=merge(final_flt_d, length_table, by.x = c("Chr"), by.y =  c("V1"), all.x = TRUE)

ppt_and_length=ppt_and_length %>%
  mutate(ltr3_1=as.numeric(as.character(V2)) - ltr_len1,
         ltr3_2=as.numeric(as.character(V2)) - ltr_len2)

ppt_and_length$Start=as.numeric(ppt_and_length$Start)
ppt_and_length$End=as.numeric(ppt_and_length$End)
ppt_and_length$Mismatches=as.numeric(ppt_and_length$Mismatches)

ppt_and_length_filter=ppt_and_length[which(ppt_and_length$Start>=ppt_and_length$ltr3_2 & ppt_and_length$Start<=ppt_and_length$ltr3_1),]

ppt_and_length_filter_new=ppt_and_length_filter[order(ppt_and_length_filter$Mismatches),]

ppt_and_length_filter_new=ppt_and_length_filter_new %>%
  distinct(Chr, Start, .keep_all = TRUE)

ppt_and_length_filter_new=as.data.table(ppt_and_length_filter_new)

ppt_final=ppt_and_length_filter_new[ppt_and_length_filter_new[, .I[which.min(Mismatches)], by=Chr]$V1]

ppt_final=as.data.frame(ppt_final)


pbs_ppt_final=merge(pbs_final, ppt_final, by.x = c("Chr"), by.y =  c("Chr"), all.x = TRUE)


# pbs_ppt_final=pbs_ppt_final[which(pbs_ppt_final$Mismatches.x<=5 & pbs_ppt_final$Mismatches.y<=5),]

pbs_ppt_final=pbs_ppt_final[which(pbs_ppt_final$Mismatches.x<=as.numeric(input8) & pbs_ppt_final$Mismatches.y<=as.numeric(input9)),]

print(pbs_ppt_final)

pbs_ppt_final=pbs_ppt_final %>%
  select(-Strand.x, -Strand.y, -Mismatches.x, -Mismatches.y, -Length.x, -Length.y, -ltr3_1, -ltr3_2)

colnames(pbs_ppt_final)=c("ID", "PBS_Seq", "PBS_Type", "PBS_Start", "PBS_End", "PPT_Seq", "PPT_Type", "PPT_Start", "PPT_End", "Length")

pbs_ppt_final=pbs_ppt_final %>%
  mutate(ID_2_destroy=ID)


# remove_ID="Atha_Evs-0"
# remove_ID=paste0(remove_ID, "_")

remove_ID=as.character(input10)
remove_ID=paste0(remove_ID, "_")


pbs_ppt_final$ID_2_destroy=gsub(remove_ID, "", pbs_ppt_final$ID_2_destroy)

# pbs_ppt_final=pbs_ppt_final %>% 
#   mutate_at("ID_2_destroy", str_replace, remove_ID, "")


pbs_ppt_final=pbs_ppt_final %>%
  separate(ID_2_destroy, c("prwto", "deutero"), sep = "-", remove=FALSE)

pbs_ppt_final=pbs_ppt_final %>%
  separate(prwto, c("chr", "coord1"), sep="\\.", remove=FALSE)


pbs_ppt_final$deutero=gsub("_Athila", "", pbs_ppt_final$deutero)

pbs_ppt_final=pbs_ppt_final %>%
  separate(deutero, c("coord2", "strand"), sep="_", remove=FALSE)


pbs_ppt_final$strand=gsub("D", "+", pbs_ppt_final$strand)
pbs_ppt_final$strand=gsub("P", "-", pbs_ppt_final$strand)

pbs_ppt_final=pbs_ppt_final %>%
  select(-prwto, -deutero)

# pbs_ppt_final$prwto=gsub(".*\\.", "", pbs_ppt_final$prwto)
# pbs_ppt_final$deutero=gsub("_.*", "", pbs_ppt_final$deutero)

D_strand=pbs_ppt_final[which(pbs_ppt_final$strand=="+"),]
P_strand=pbs_ppt_final[which(pbs_ppt_final$strand=="-"),]


if (dim(D_strand)[1] == 0) {
  
  write.table(D_strand_5ltr, input11, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(D_strand_3ltr, input12, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(D_strand_internal, input13, sep = "\t", row.names = F, quote= F, col.names = F)
  
} else {
  
  D_strand$PBS_Start=as.numeric(as.character(D_strand$PBS_Start)) + 5
  D_strand$PPT_End=as.numeric(as.character(D_strand$PPT_End)) - 5
  
  
  D_strand=D_strand %>%
    mutate(ltr5_end=as.numeric(as.character(coord1)) + as.numeric(as.character(PBS_Start)),
           ltr3_start=as.numeric(as.character(coord1)) + as.numeric(as.character(PPT_End)))
  
  
  D_strand_5ltr=D_strand %>%
    select(chr, coord1, ltr5_end, ID, PBS_Seq, PPT_Seq, strand)
  
  D_strand_5ltr$ID=paste(D_strand_5ltr$ID, "5prime", sep="_")
  
  D_strand_3ltr=D_strand %>%
    select(chr, ltr3_start, coord2, ID, PBS_Seq, PPT_Seq, strand)
  
  D_strand_3ltr$ID=paste(D_strand_3ltr$ID, "3prime", sep="_")
  
  D_strand_internal=D_strand %>%
    select(chr, ltr5_end, ltr3_start, ID, PBS_Seq, PPT_Seq, strand)
  
  D_strand_internal$ID=paste(D_strand_internal$ID, "internal", sep="_")
  
  
  write.table(D_strand_5ltr, input11, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(D_strand_3ltr, input12, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(D_strand_internal, input13, sep = "\t", row.names = F, quote= F, col.names = F)
  
}


# D_strand$PBS_Start=as.numeric(as.character(D_strand$PBS_Start)) + 5
# D_strand$PPT_End=as.numeric(as.character(D_strand$PPT_End)) - 5
# 
# 
# D_strand=D_strand %>%
#   mutate(ltr5_end=as.numeric(as.character(coord1)) + as.numeric(as.character(PBS_Start)),
#          ltr3_start=as.numeric(as.character(coord1)) + as.numeric(as.character(PPT_End)))
# 
# 
# D_strand_5ltr=D_strand %>%
#   select(chr, coord1, ltr5_end, ID, PBS_Seq, PPT_Seq, strand)
# 
# D_strand_5ltr$ID=paste(D_strand_5ltr$ID, "5prime", sep="_")
# 
# D_strand_3ltr=D_strand %>%
#   select(chr, ltr3_start, coord2, ID, PBS_Seq, PPT_Seq, strand)
# 
# D_strand_3ltr$ID=paste(D_strand_3ltr$ID, "3prime", sep="_")
# 
# D_strand_internal=D_strand %>%
#   select(chr, ltr5_end, ltr3_start, ID, PBS_Seq, PPT_Seq, strand)
# 
# D_strand_internal$ID=paste(D_strand_internal$ID, "internal", sep="_")
# 
# 
# write.table(D_strand_5ltr, input11, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(D_strand_3ltr, input12, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(D_strand_internal, input13, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(D_strand_5ltr, "t2t-col.20210610.fasta.D_5LTR_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(D_strand_3ltr, "t2t-col.20210610.fasta.D_3LTR_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(D_strand_internal, "t2t-col.20210610.fasta.D_INTERNAL_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)




if (dim(P_strand)[1] == 0) {
  
  write.table(P_strand, input14, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(P_strand, input15, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(P_strand, input16, sep = "\t", row.names = F, quote= F, col.names = F)
  
} else {
  P_strand$PBS_Start=as.numeric(as.character(P_strand$PBS_Start)) + 5
  P_strand$PPT_End=as.numeric(as.character(P_strand$PPT_End)) - 5
  
  
  P_strand=P_strand %>%
    mutate(ltr5_end=as.numeric(as.character(coord2)) - as.numeric(as.character(PBS_Start)),
           ltr3_start=as.numeric(as.character(coord2)) - as.numeric(as.character(PPT_End)))
  
  
  P_strand_5ltr=P_strand %>%
    select(chr, ltr5_end, coord2, ID, PBS_Seq, PPT_Seq, strand)
  
  P_strand_5ltr$ID=paste(P_strand_5ltr$ID, "5prime", sep="_")
  
  P_strand_3ltr=P_strand %>%
    select(chr, coord1, ltr3_start, ID, PBS_Seq, PPT_Seq, strand)
  
  P_strand_3ltr$ID=paste(P_strand_3ltr$ID, "3prime", sep="_")
  
  P_strand_internal=P_strand %>%
    select(chr, ltr3_start, ltr5_end, ID, PBS_Seq, PPT_Seq, strand)
  
  P_strand_internal$ID=paste(P_strand_internal$ID, "internal", sep="_")
  
  write.table(P_strand_5ltr, input14, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(P_strand_3ltr, input15, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(P_strand_internal, input16, sep = "\t", row.names = F, quote= F, col.names = F)
}


# P_strand$PBS_Start=as.numeric(as.character(P_strand$PBS_Start)) + 5
# P_strand$PPT_End=as.numeric(as.character(P_strand$PPT_End)) - 5
# 
# 
# P_strand=P_strand %>%
#   mutate(ltr5_end=as.numeric(as.character(coord2)) - as.numeric(as.character(PBS_Start)),
#          ltr3_start=as.numeric(as.character(coord2)) - as.numeric(as.character(PPT_End)))
# 
# 
# P_strand_5ltr=P_strand %>%
#   select(chr, ltr5_end, coord2, ID, PBS_Seq, PPT_Seq, strand)
# 
# P_strand_5ltr$ID=paste(P_strand_5ltr$ID, "5prime", sep="_")
# 
# P_strand_3ltr=P_strand %>%
#   select(chr, coord1, ltr3_start, ID, PBS_Seq, PPT_Seq, strand)
# 
# P_strand_3ltr$ID=paste(P_strand_3ltr$ID, "3prime", sep="_")
# 
# P_strand_internal=P_strand %>%
#   select(chr, ltr3_start, ltr5_end, ID, PBS_Seq, PPT_Seq, strand)
# 
# P_strand_internal$ID=paste(P_strand_internal$ID, "internal", sep="_")
#
# 
# write.table(P_strand_5ltr, input14, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(P_strand_3ltr, input15, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(P_strand_internal, input16, sep = "\t", row.names = F, quote= F, col.names = F)


# write.table(P_strand_5ltr, "t2t-col.20210610.fasta.P_5LTR_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(P_strand_3ltr, "t2t-col.20210610.fasta.P_3LTR_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(P_strand_internal, "t2t-col.20210610.fasta.P_INTERNAL_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)

