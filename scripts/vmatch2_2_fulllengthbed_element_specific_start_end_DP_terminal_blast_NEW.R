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


################################# F O R  D  S T R A N D #################################


######### START D STRAND

dokimh_st_d=read.table(input1, fill = T)

# dokimh_st_d=read.table("Col-CEN_v1.2_clean.fasta_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_START_D_mis15.vmatch", fill = T)

flt1_st_d=dokimh_st_d[which(dokimh_st_d$V1!="Query:" & dokimh_st_d$V1!="!" & dokimh_st_d$V1!="!!" & dokimh_st_d$V1!="!!!" & 
                            dokimh_st_d$V1!="!!!!"),]

flt2_st_d=flt1_st_d[which(flt1_st_d$V1=="15" | flt1_st_d$V1=="20" | flt1_st_d$V1=="25"),]

colnames(flt2_st_d)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_st_d=flt1_st_d[which(flt1_st_d$V1=="Sbjct:"),]

colnames(flt3_st_d)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_st_d$Start=as.numeric(as.character(flt2_st_d$Start))
flt3_st_d$End=as.numeric(as.character(flt3_st_d$End))

flt4_st_d=cbind(flt2_st_d, flt3_st_d)

final_flt_st_d=flt4_st_d %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_st_d)[8]="Length"



telos=final_flt_st_d

ltr_min=2000-as.numeric(input2)
ltr_max=2000-as.numeric(input3)

# ltr_max=2000-2500
# ltr_min=2000-500

telos_new=telos[which(telos$Start>=ltr_max & telos$Start<=ltr_min),] 

telos_new1=telos_new %>%
  mutate(ID_2_destroy=Chr)

telos_new1$ID_2_destroy=gsub(".*:", "", telos_new1$ID_2_destroy)
telos_new1$ID_2_destroy=gsub("_P_element", "", telos_new1$ID_2_destroy)

telos_new2=telos_new1 %>%
  separate(Type, c("EXT_ID", "BLAST_ID"), "\\$")

telos_new2$BLAST_ID=gsub("_P_element", "", telos_new2$BLAST_ID)
telos_new2$EXT_ID=gsub("rc_", "", telos_new2$EXT_ID)

final_telos=telos_new2[which(telos_new2$EXT_ID==telos_new2$ID_2_destroy),]

final_telos$Mismatches=as.numeric(final_telos$Mismatches)
final_telos=final_telos[order(-final_telos$Mismatches),]
final_telos$Mismatches=abs(final_telos$Mismatches)

tr=final_telos %>%
  arrange(Chr, Mismatches, Start) %>%
  group_by(Chr) %>%
  slice_head(n=1)

tr=as.data.frame(tr)
# tr
# print("start")


######### END D STRAND

dokimh_en_d=read.table(input4, fill = T)

# dokimh_en_d=read.table("Col-CEN_v1.2_clean.fasta_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_END_D_mis15.vmatch", fill = T)

flt1_en_d=dokimh_en_d[which(dokimh_en_d$V1!="Query:" & dokimh_en_d$V1!="!" & dokimh_en_d$V1!="!!" & dokimh_en_d$V1!="!!!" & 
                            dokimh_en_d$V1!="!!!!"),]

flt2_en_d=flt1_en_d[which(flt1_en_d$V1=="15" | flt1_en_d$V1=="20" | flt1_en_d$V1=="25"),]

colnames(flt2_en_d)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_en_d=flt1_en_d[which(flt1_en_d$V1=="Sbjct:"),]

colnames(flt3_en_d)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_en_d$Start=as.numeric(as.character(flt2_en_d$Start))
flt3_en_d$End=as.numeric(as.character(flt3_en_d$End))

flt4_en_d=cbind(flt2_en_d, flt3_en_d)

final_flt_en_d=flt4_en_d %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_en_d)[8]="Length"


len_tab=read.table(input5, sep="\t")

# len_tab=read.table("Col-CEN_v1.2_clean.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING_EXTENDED_LENGTH.txt", sep="\t")

temp1=merge(final_flt_en_d, len_tab, by.x = c("Chr"), by.y=c("V1"), all.x=TRUE)

temp2=temp1 %>%
  mutate(three_LTR=as.numeric(as.character(V2)) - 2000)

final_flt_en_d_new=temp2[which(temp2$Start>=temp2$three_LTR),]



end=final_flt_en_d_new

end_new=end %>%
  mutate(LTR3_a=as.numeric(as.character(V2)) - ltr_min,
         LTR3_b=as.numeric(as.character(V2)) - ltr_max)

end_new=end_new[which(end_new$Start>=end_new$LTR3_a & end_new$Start<=end_new$LTR3_b),]

end_new2=end_new %>%
  mutate(ID_2_destroy=Chr)

end_new2$ID_2_destroy=gsub(".*:", "", end_new2$ID_2_destroy)
end_new2$ID_2_destroy=gsub("_P_element", "", end_new2$ID_2_destroy)

end_new3=end_new2 %>%
  separate(Type, c("EXT_ID", "BLAST_ID"), "\\$")

end_new3$BLAST_ID=gsub("_P_element", "", end_new3$BLAST_ID)
end_new3$EXT_ID=gsub("rc_", "", end_new3$EXT_ID)


final_end=end_new3[which(end_new3$EXT_ID==end_new3$ID_2_destroy),]

final_end$Mismatches=as.numeric(final_end$Mismatches)
final_end1=final_end[order(-final_end$Mismatches),]
final_end1$Mismatches=abs(final_end1$Mismatches)

skata=final_end1 %>%
  arrange(Chr, Mismatches, Start) %>%
  group_by(Chr) %>%
  slice_head(n=1)

# skata <- final_end1 %>%
#   arrange(Chr, Mismatches, Start) %>%
#   group_by(Chr) %>%
#   top_n(1, wt = Start)

skata2=as.data.frame(skata)
# skata2
# rt
# print("end")


######################## SUPER MERGING

mega_tab=merge(tr, skata2, by.x = c("Chr"), by.y=c("Chr"), all.x=TRUE)
# mega_tab=na.omit(mega_tab)
# mega_tab
mega_tab1=mega_tab %>%
  select(Chr, Seq.x, Start.x, Seq.y, End.y, BLAST_ID.x)

colnames(mega_tab1)=c("Lineage_ID", "Seq_Start", "Lineage_Start", "Seq_End", "Lineage_End", "BLAST_ID")

mega_tab2=mega_tab1 %>%
  select(Lineage_ID, Seq_Start, Seq_End, Lineage_Start, Lineage_End, BLAST_ID)

new_mega_tab=mega_tab2 %>%
  separate(Lineage_ID, c("chr", "coords"), sep = ":", remove=FALSE)

new_mega_tab1=new_mega_tab %>%
  separate(coords, c("Start_Elem", "End_Elem"), sep="-", remove=FALSE)

new_mega_tab1$chr=gsub("rc_", "", new_mega_tab1$chr)
new_mega_tab1$End_Elem=gsub("_P_element", "", new_mega_tab1$End_Elem)

P_STRAND=new_mega_tab1[new_mega_tab1$Lineage_ID %like% "_P_element", ]

suppressMessages({D_STRAND=anti_join(new_mega_tab1, P_STRAND)})

NEW_D_STRAND=D_STRAND %>%
  mutate(New_Start_Elem=as.numeric(as.character(Start_Elem)) + as.numeric(as.character(Lineage_Start)),
         New_End_Elem=as.numeric(as.character(Start_Elem)) + as.numeric(as.character(Lineage_End)))

NEW_P_STRAND=P_STRAND %>%
  mutate(New_End_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(Lineage_Start)),
         New_Start_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(Lineage_End)))


add_ID=as.character(input6)

# add_ID="Atha_Col-0"

NEW_D_STRAND1=NEW_D_STRAND %>%
  mutate(Strand="+", New_ID_1=add_ID, New_ID_2=paste0("D_", as.character(input7)))

NEW_P_STRAND1=NEW_P_STRAND %>%
  mutate(Strand="-", New_ID_1=add_ID, New_ID_2=paste0("P_", as.character(input7)))

new2_mega_tab_D=NEW_D_STRAND1 %>%
  select(chr, New_Start_Elem, New_End_Elem, Lineage_ID, Seq_Start, Seq_End, Strand, New_ID_1, New_ID_2, BLAST_ID)

new2_mega_tab_D=unite(new2_mega_tab_D, ID_NEW, c(New_Start_Elem, New_End_Elem), remove=FALSE, sep="-")
new2_mega_tab_D=unite(new2_mega_tab_D, ID_NEW, c(chr, ID_NEW), remove=FALSE, sep=".")
new2_mega_tab_D=unite(new2_mega_tab_D, ID_NEW, c(New_ID_1, ID_NEW), remove=FALSE, sep="_")
new2_mega_tab_D=unite(new2_mega_tab_D, ID_NEW, c(ID_NEW, New_ID_2), remove=FALSE, sep="_")

new2_mega_tab_D2=new2_mega_tab_D %>%
  select(chr, New_Start_Elem, New_End_Elem, ID_NEW, Seq_Start, Seq_End, BLAST_ID, Strand)

# new2_mega_tab_D2

new2_mega_tab_P=NEW_P_STRAND1 %>%
  select(chr, New_Start_Elem, New_End_Elem, Lineage_ID, Seq_Start, Seq_End, Strand, New_ID_1, New_ID_2, BLAST_ID)

new2_mega_tab_P=unite(new2_mega_tab_P, ID_NEW, c(New_Start_Elem, New_End_Elem), remove=FALSE, sep="-")
new2_mega_tab_P=unite(new2_mega_tab_P, ID_NEW, c(chr, ID_NEW), remove=FALSE, sep=".")
new2_mega_tab_P=unite(new2_mega_tab_P, ID_NEW, c(New_ID_1, ID_NEW), remove=FALSE, sep="_")
new2_mega_tab_P=unite(new2_mega_tab_P, ID_NEW, c(ID_NEW, New_ID_2), remove=FALSE, sep="_")

new2_mega_tab_P2=new2_mega_tab_P %>%
  select(chr, New_Start_Elem, New_End_Elem, ID_NEW, Seq_Start, Seq_End, BLAST_ID, Strand)

# new2_mega_tab_P2

# write.table(new2_mega_tab_D2, "Col-CEN_v1.2.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_SPECIFIC_ENDS.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new2_mega_tab_D2, input8, sep = "\t", row.names = F, quote= F, col.names = F)

