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
input17 = args[17]


############ LENGTH TABLE FOR BLAST ELEMENTS

length_table=read.table(input1, sep="\t")

# length_table=read.table("Col-CEN_v1.2.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING_LENGTH.txt", sep="\t")

############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PBS EXCEPT 4 and 8

vmatch_tab_a=read.table(input2, fill = T)

# vmatch_tab_a=read.table("Col-CEN_v1.2.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PBS.vmatch", fill=T)

flt1_a=vmatch_tab_a[which(vmatch_tab_a$V1!="Query:" & vmatch_tab_a$V1!="!" & vmatch_tab_a$V1!="!!"),]

# flt2_start_a=flt1_a[which(flt1_a$V1=="20" | flt1_a$V1=="21"),]

suppressWarnings({flt2_start_a=flt1_a[!is.na(as.numeric(flt1_a$V1)),]})

colnames(flt2_start_a)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_end_a=flt1_a[which(flt1_a$V1=="Sbjct:"),]

colnames(flt3_end_a)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt4_a=cbind(flt2_start_a, flt3_end_a)

final_flt_a=flt4_a %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

final_flt_a$Mismatches=gsub("-", "", final_flt_a$Mismatches)

colnames(final_flt_a)[8]="Length"


pbs_table=final_flt_a

pbs_table$Start=as.numeric(pbs_table$Start)
pbs_table$End=as.numeric(pbs_table$End)
pbs_table$Mismatches=as.numeric(pbs_table$Mismatches)

ltr_len1=as.numeric(input3)
ltr_len2=as.numeric(input4)

# ltr_len1=500
# ltr_len2=2500

pbs_length=pbs_table[which(pbs_table$Start>=ltr_len1 & pbs_table$Start<=ltr_len2),]

pbs_length_new=pbs_length[order(pbs_length$Mismatches),]

pbs_length_new=pbs_length_new %>%
  distinct(Chr, Start, .keep_all = TRUE)

pbs_length_new=as.data.table(pbs_length_new)

pbs_final=pbs_length_new[pbs_length_new[, .I[which.min(Mismatches)], by=Chr]$V1]

pbs_final=as.data.frame(pbs_final)


############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS --- PPT

vmatch_tab_d=read.table(input5, fill = T)

# vmatch_tab_d=read.table("Col-CEN_v1.2.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING_PPT.vmatch", fill=T)

flt1_d=vmatch_tab_d[which(vmatch_tab_d$V1!="Query:" & vmatch_tab_d$V1!="!" & vmatch_tab_d$V1!="!!"),]

# flt2_start_d=flt1_d[which(flt1_d$V1=="20" | flt1_d$V1=="21"),]

suppressWarnings({flt2_start_d=flt1_d[!is.na(as.numeric(flt1_d$V1)),]})

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

pbs_ppt_final=pbs_ppt_final[which(pbs_ppt_final$Mismatches.x<=as.numeric(input6) & pbs_ppt_final$Mismatches.y<=as.numeric(input7)),]

# pbs_ppt_final=pbs_ppt_final[which(pbs_ppt_final$Mismatches.x<=5 & pbs_ppt_final$Mismatches.y<=5),]

# print(pbs_ppt_final)

pbs_ppt_final=pbs_ppt_final %>%
  select(-Strand.x, -Strand.y, -Length.x, -Length.y, -ltr3_1, -ltr3_2)

colnames(pbs_ppt_final)=c("ID", "PBS_Seq", "PBS_Type", "PBS_Mismatches", "PBS_Start", "PBS_End", "PPT_Seq", "PPT_Type", "PPT_Mismatches", 
                          "PPT_Start", "PPT_End", "Length")

pbs_ppt_final=pbs_ppt_final %>%
  mutate(ID_2_destroy=ID)


remove_ID=as.character(input8)
remove_ID=paste0(remove_ID, "_")

# remove_ID="Atha_Col-0"
# remove_ID=paste0(remove_ID, "_")

pbs_ppt_final$ID_2_destroy=gsub(remove_ID, "", pbs_ppt_final$ID_2_destroy)

# pbs_ppt_final=pbs_ppt_final %>% 
#   mutate_at("ID_2_destroy", str_replace, remove_ID, "")


pbs_ppt_final=pbs_ppt_final %>%
  separate(ID_2_destroy, c("prwto", "deutero"), sep = "-", remove=FALSE)

pbs_ppt_final=pbs_ppt_final %>%
  separate(prwto, c("chr", "coord1"), sep="\\.", remove=FALSE)

pbs_ppt_final$deutero=gsub("(_[PD]).*$", "\\1", pbs_ppt_final$deutero)

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

oligomer=as.numeric(input9) - 5

# oligomer=20 - 5


D_strand=D_strand %>%
  mutate(pos_st_st=as.numeric(as.character(coord1)) + as.numeric(as.character(PBS_Start)) - oligomer,
         pos_st_en=as.numeric(as.character(coord1)) + as.numeric(as.character(PBS_Start)) + 5,
         pos_en_st=as.numeric(as.character(coord1)) + as.numeric(as.character(PPT_End)) - 5,
         pos_en_en=as.numeric(as.character(coord1)) + as.numeric(as.character(PPT_End)) + oligomer,
         ext_st=as.numeric(as.character(coord1)) + as.numeric(as.character(PBS_Start)) + 5 - 2000,
         ext_en=as.numeric(as.character(coord1)) + as.numeric(as.character(PPT_End)) - 5 + 2000)

D_strand=unite(D_strand, start_or_end_id, c(ext_st, ext_en), remove=FALSE, sep="-")
D_strand=unite(D_strand, start_or_end_id, c(start_or_end_id, ID), remove=FALSE, sep="$")



D_strand_pos_st=D_strand %>%
  select(chr, pos_en_st, pos_en_en, start_or_end_id, PBS_Type, PBS_Seq, PBS_Mismatches, strand)

write.table(D_strand_pos_st, input10, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(D_strand_pos_st, "Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_POTENTIAL_STARTS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

D_strand_pos_en=D_strand %>%
  select(chr, pos_st_st, pos_st_en, start_or_end_id, PPT_Type, PPT_Seq, PPT_Mismatches, strand)

write.table(D_strand_pos_en, input11, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(D_strand_pos_en, "Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_POTENTIAL_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

D_strand_ext=D_strand %>%
  select(chr, ext_st, ext_en, ID, PBS_Type, PBS_Seq, PBS_Mismatches, PPT_Type, PPT_Seq, PPT_Mismatches, strand) # I removed Length column

write.table(D_strand_ext, input12, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(D_strand_ext, "Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_EXTENDED.bed", sep = "\t", row.names = F, quote= F, col.names = F)

D_strand_int=D_strand %>%
  select(chr, pos_st_en, pos_en_st, ID, PBS_Type, PBS_Seq, PBS_Mismatches, PPT_Type, PPT_Seq, PPT_Mismatches, strand) # I removed Length column

write.table(D_strand_int, input13, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(D_strand_int, "Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_INTERNAL.bed", sep = "\t", row.names = F, quote= F, col.names = F)



P_strand=P_strand %>%
  mutate(pos_st_st=as.numeric(as.character(coord2)) - as.numeric(as.character(PBS_Start)) - 5,
         pos_st_en=as.numeric(as.character(coord2)) - as.numeric(as.character(PBS_Start)) + oligomer,
         pos_en_st=as.numeric(as.character(coord2)) - as.numeric(as.character(PPT_End)) - oligomer,
         pos_en_en=as.numeric(as.character(coord2)) - as.numeric(as.character(PPT_End)) + 5,
         ext_st=as.numeric(as.character(coord2)) - as.numeric(as.character(PBS_Start)) - 5 + 2000,
         ext_en=as.numeric(as.character(coord2)) - as.numeric(as.character(PPT_End)) + 5 - 2000)

P_strand=unite(P_strand, start_or_end_id, c(ext_en, ext_st), remove=FALSE, sep="-")
P_strand=unite(P_strand, start_or_end_id, c(start_or_end_id, ID), remove=FALSE, sep="$")


P_strand_pos_st=P_strand %>%
  select(chr, pos_en_st, pos_en_en, start_or_end_id, PBS_Type, PBS_Seq, PBS_Mismatches, strand)

write.table(P_strand_pos_st, input14, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(P_strand_pos_st, "Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_POTENTIAL_STARTS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

P_strand_pos_en=P_strand %>%
  select(chr, pos_st_st, pos_st_en, start_or_end_id, PPT_Type, PPT_Seq, PPT_Mismatches, strand)

write.table(P_strand_pos_en, input15, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(P_strand_pos_en, "Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_POTENTIAL_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

P_strand_ext=P_strand %>%
  select(chr, ext_en, ext_st, ID, PPT_Type, PPT_Seq, PPT_Mismatches, PBS_Type, PBS_Seq, PBS_Mismatches, strand)

write.table(P_strand_ext, input16, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(P_strand_ext, "Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_EXTENDED.bed", sep = "\t", row.names = F, quote= F, col.names = F)

P_strand_int=P_strand %>%
  select(chr, pos_en_en, pos_st_st, ID, PPT_Type, PPT_Seq, PPT_Mismatches, PBS_Type, PBS_Seq, PBS_Mismatches, strand)

write.table(P_strand_int, input17, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(P_strand_int, "Col-CEN_v1.2.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_INTERNAL.bed", sep = "\t", row.names = F, quote= F, col.names = F)

