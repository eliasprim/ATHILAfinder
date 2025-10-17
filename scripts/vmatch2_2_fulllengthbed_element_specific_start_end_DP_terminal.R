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
#
# input13 = args[13]
# input14 = args[14]


# # STEP ARGUMENT FOR BASH SCRIPT STEP 5

################################# F O R  D  S T R A N D #################################


######### START D STRAND

dokimh_st_d=read.table(input1, fill = T)

# dokimh_st_d=read.table("t2t-col.20210610.fasta.F2B_ATHILA_SPECIFIC_START_D_mis4.vmatch", fill = T)

flt1_st_d=dokimh_st_d[which(dokimh_st_d$V1!="Query:" & dokimh_st_d$V1!="!" & dokimh_st_d$V1!="!!" & dokimh_st_d$V1!="!!!" & 
                            dokimh_st_d$V1!="!!!!"),]

# flt2_st_d=flt1_st_d[which(flt1_st_d$V1=="15" | flt1_st_d$V1=="20" | flt1_st_d$V1=="25"),]

suppressWarnings({flt2_st_d=flt1_st_d[!is.na(as.numeric(flt1_st_d$V1)),]})

colnames(flt2_st_d)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_st_d=flt1_st_d[which(flt1_st_d$V1=="Sbjct:"),]

colnames(flt3_st_d)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_st_d$Start=as.numeric(as.character(flt2_st_d$Start))
flt3_st_d$End=as.numeric(as.character(flt3_st_d$End))

flt4_st_d=cbind(flt2_st_d, flt3_st_d)

final_flt_st_d=flt4_st_d %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_st_d)[8]="Length"


mis01=final_flt_st_d[which(final_flt_st_d$Start<as.numeric(input4)),] ### 92/109 20mer/4mis, 105/109 15mer/3mis, 90/109 25mer/5mis
mis01_new=mis01[order(-mis01$Mismatches),]

new_mis01=mis01_new %>%
  distinct(Chr, Start, .keep_all = TRUE)




######### START P STRAND

dokimh_st_p=read.table(input2, fill = T)

# dokimh_st_p=read.table("t2t-col.20210610.fasta.F2B_ATHILA_SPECIFIC_START_P_mis4.vmatch", fill = T)

flt1_st_p=dokimh_st_p[which(dokimh_st_p$V1!="Query:" & dokimh_st_p$V1!="!" & dokimh_st_p$V1!="!!" & dokimh_st_p$V1!="!!!" & 
                              dokimh_st_p$V1!="!!!!"),]

# flt2_st_p=flt1_st_p[which(flt1_st_p$V1=="15" | flt1_st_p$V1=="20" | flt1_st_p$V1=="25"),]

suppressWarnings({flt2_st_p=flt1_st_p[!is.na(as.numeric(flt1_st_p$V1)),]})

colnames(flt2_st_p)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_st_p=flt1_st_p[which(flt1_st_p$V1=="Sbjct:"),]

colnames(flt3_st_p)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_st_p$Start=as.numeric(as.character(flt2_st_p$Start))
flt3_st_p$End=as.numeric(as.character(flt3_st_p$End))

flt4_st_p=cbind(flt2_st_p, flt3_st_p)

final_flt_st_p=flt4_st_p %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_st_p)[8]="Length"


mis02=final_flt_st_p[which(final_flt_st_p$Start<as.numeric(input4)),] ### 77/97 20mer/4mis, 90/97 15mer/3mis, 75/97 25mer/5mer
mis02_new=mis02[order(-mis02$Mismatches),]

new_mis02=mis02_new %>%
  distinct(Chr, Start, .keep_all = TRUE)



aporipsh_arxhs=rbind(final_flt_st_d, final_flt_st_p)

telos=rbind(new_mis01, new_mis02)



# d=as.data.frame(table(telos$Mismatches))
# d=d[order(d$Freq, decreasing = TRUE),]
# d
# 
# telos$Mismatches=gsub("0", "0_118", telos$Mismatches)
# telos$Mismatches=gsub("-3", "3_4", telos$Mismatches)
# telos$Mismatches=gsub("-4", "4_3", telos$Mismatches)
# telos$Mismatches=gsub("-2", "2_14", telos$Mismatches)
# telos$Mismatches=gsub("-1", "1_30", telos$Mismatches)
# telos$Mismatches=gsub("-5", "5_18", telos$Mismatches)
# 
# ggplot(telos, aes(Mismatches, Start, fill=Mismatches)) +
#   geom_violin(width=1) +
#   geom_boxplot(width=0.2, color="black", alpha=0.8) +
#   theme(legend.position="none", plot.title = element_text(size=11)) +
#   ggtitle("") +
#   ylab("Start of k-mer") +
#   xlab("Mismatches_Times") +
#   ggtitle("Element-specific 25mer starts with up to 5 mismatches 165/206 elements")




### telos_new=telos[which(telos$Start>=as.numeric(input3) & telos$Start<=as.numeric(input4)),]

ltr_max=2000-as.numeric(input4)
ltr_min=2000-as.numeric(input3)

# ltr_max=2000-1900
# ltr_min=2000-1100

telos_new=telos[which(telos$Start>=ltr_max & telos$Start<=ltr_min),] ### 159/169 20mer/4mis, 172/195 elements 15mer/3mis, 156/165 25mer/5mis


duplicate_telos_new=telos_new %>%
  group_by(Chr) %>%
    filter(n()>1)

duplicate_telos_new=as.data.frame(duplicate_telos_new)



suppressMessages({unique_telos_new=anti_join(telos_new, duplicate_telos_new)})  ### merge to the final start table



skt=read.table(input5, sep='\t')

# skt=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_POTENTIAL_STARTS.txt", sep='\t')

super_test=merge(duplicate_telos_new, skt, by.x = c("Chr"), by.y = c("V1"), all.x = TRUE)

super_test$Seq=as.character(super_test$Seq)
super_test$V2=as.character(super_test$V2)

ert=super_test[which(super_test$Seq==super_test$V2),]
ert=as.data.table(ert)

new_ert=ert[ert[, .I[which.min(Mismatches)], by=Chr]$V1] 
new_ert=as.data.frame(new_ert)

new2_ert=new_ert %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length) ### merge to the final start table

suppressMessages({skta=anti_join(super_test, ert)})

jk=intersect(ert$Chr, skta$Chr)
jk=as.data.frame(jk)

clear_skta=skta[!skta$Chr %in% jk$jk, , drop = FALSE]
clear_skta=as.data.table(clear_skta)

sure_start_st=clear_skta[clear_skta[, .I[which.min(Mismatches)], by=Chr]$V1]

# sure_start_st=NULL;
# if (as.character(input5)=="maximum"){
#   sure_start_st=telos_new[passed_start[, .I[which.min(Start)], by=Chr]$V1]
# } else {
#   sure_start_st=telos_new[passed_start[, .I[which.max(Start)], by=Chr]$V1]
# }

sure_start_st=as.data.frame(sure_start_st)

sure_start_st=sure_start_st %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length) ### merge to the final start table

final_start_table=rbind(unique_telos_new, new2_ert, sure_start_st)


# d1=as.data.frame(table(sure_start_st$Mismatches))
# d1=d1[order(d1$Freq, decreasing = TRUE),]
# d1
# 
# sure_start_st$Mismatches=gsub("0", "0_128", sure_start_st$Mismatches)
# sure_start_st$Mismatches=gsub("-4", "4_16", sure_start_st$Mismatches)
# sure_start_st$Mismatches=gsub("-2", "2_7", sure_start_st$Mismatches)
# sure_start_st$Mismatches=gsub("-3", "3_3", sure_start_st$Mismatches)
# sure_start_st$Mismatches=gsub("-1", "1_32", sure_start_st$Mismatches)
# 
# ggplot(sure_start_st, aes(Mismatches, Start, fill=Mismatches)) +
#   geom_violin(width=1) +
#   geom_boxplot(width=0.2, color="black", alpha=0.8) +
#   theme(legend.position="none", plot.title = element_text(size=11)) +
#   ggtitle("") +
#   ylab("Start of k-mer") +
#   xlab("Mismatches_Times") +
#   ggtitle("Element-specific 20mer starts with up to 4 mismatches 159/169 elements, 1100-1900 LTR length")






######### END D STRAND

dokimh_en_d=read.table(input6, fill = T)

# dokimh_en_d=read.table("t2t-col.20210610.fasta.F2B_ATHILA_SPECIFIC_END_D_mis4.vmatch", fill = T)

flt1_en_d=dokimh_en_d[which(dokimh_en_d$V1!="Query:" & dokimh_en_d$V1!="!" & dokimh_en_d$V1!="!!" & dokimh_en_d$V1!="!!!" & 
                            dokimh_en_d$V1!="!!!!"),]

# flt2_en_d=flt1_en_d[which(flt1_en_d$V1=="15" | flt1_en_d$V1=="20" | flt1_en_d$V1=="25"),]

suppressWarnings({flt2_en_d=flt1_en_d[!is.na(as.numeric(flt1_en_d$V1)),]})

colnames(flt2_en_d)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_en_d=flt1_en_d[which(flt1_en_d$V1=="Sbjct:"),]

colnames(flt3_en_d)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_en_d$Start=as.numeric(as.character(flt2_en_d$Start))
flt3_en_d$End=as.numeric(as.character(flt3_en_d$End))

flt4_en_d=cbind(flt2_en_d, flt3_en_d)

final_flt_en_d=flt4_en_d %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_en_d)[8]="Length"


final_flt_en_d=final_flt_en_d %>%
  separate(Chr, c("chr", "coords"), sep = ":", remove=FALSE)

final_flt_en_d=final_flt_en_d %>%
  separate(coords, c("Start_Elem", "End_Elem"), sep="-", remove=TRUE)

final_flt_en_d=final_flt_en_d %>%
  mutate(Length_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(Start_Elem)))

final_flt_en_d=final_flt_en_d %>%
  mutate(three_LTR=as.numeric(as.character(Length_Elem)) - as.numeric(input4)) 

mis03=final_flt_en_d[which(final_flt_en_d$Start>=final_flt_en_d$three_LTR),] ### 92/109 20mer/4mis, 105/109 15mer/3mis, 92/109 25mer/mis

mis03_new=mis03[order(-mis03$Mismatches),]

new_mis03=mis03_new %>%
  distinct(Chr, Start, .keep_all = TRUE)



######### END P STRAND

dokimh_en_p=read.table(input7, fill = T)

# dokimh_en_p=read.table("t2t-col.20210610.fasta.F2B_ATHILA_SPECIFIC_END_P_mis4.vmatch", fill = T)

flt1_en_p=dokimh_en_p[which(dokimh_en_p$V1!="Query:" & dokimh_en_p$V1!="!" & dokimh_en_p$V1!="!!" & dokimh_en_p$V1!="!!!" & 
                            dokimh_en_p$V1!="!!!!"),]

# flt2_en_p=flt1_en_p[which(flt1_en_p$V1=="15" | flt1_en_p$V1=="20" | flt1_en_p$V1=="25"),]

suppressWarnings({flt2_en_p=flt1_en_p[!is.na(as.numeric(flt1_en_p$V1)),]})

colnames(flt2_en_p)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")

flt3_en_p=flt1_en_p[which(flt1_en_p$V1=="Sbjct:"),]

colnames(flt3_en_p)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")

flt2_en_p$Start=as.numeric(as.character(flt2_en_p$Start))
flt3_en_p$End=as.numeric(as.character(flt3_en_p$End))

flt4_en_p=cbind(flt2_en_p, flt3_en_p)

final_flt_en_p=flt4_en_p %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)

colnames(final_flt_en_p)[8]="Length"


final_flt_en_p=final_flt_en_p %>%
  separate(Chr, c("chr", "coords"), sep = ":", remove=FALSE)

final_flt_en_p=final_flt_en_p %>%
  separate(coords, c("Start_Elem", "End_Elem"), sep="-", remove=TRUE)

final_flt_en_p$End_Elem=gsub("_P_element", "", final_flt_en_p$End_Elem)

final_flt_en_p=final_flt_en_p %>%
  mutate(Length_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(Start_Elem)))

final_flt_en_p=final_flt_en_p %>%
  mutate(three_LTR=as.numeric(as.character(Length_Elem)) - as.numeric(input4)) 

mis04=final_flt_en_p[which(final_flt_en_p$Start>=final_flt_en_p$three_LTR),] ### 88/97 20mer/4mis, 95/97 15mer/3mis, 80/97 25mer/5mis

mis04_new=mis04[order(-mis04$Mismatches),]

new_mis04=mis04_new %>%
  distinct(Chr, Start, .keep_all = TRUE)



aporipsh_telous=rbind(final_flt_en_d, final_flt_en_p)

end=rbind(new_mis03, new_mis04)

# ltr_max=2000-as.numeric(input4) 1900
# ltr_min=2000-as.numeric(input3) 1100

### end=end %>%
###   mutate(LTR3_a=as.numeric(as.character(Length_Elem)) - as.numeric(input4),
###          LTR3_b=as.numeric(as.character(Length_Elem)) - as.numeric(input3)) ### arguments here

end=end %>%
  mutate(LTR3_a=as.numeric(as.character(Length_Elem)) - ltr_min,
         LTR3_b=as.numeric(as.character(Length_Elem)) - ltr_max) ### arguments here


end_new=end[which(end$Start>=end$LTR3_a & end$Start<=end$LTR3_b),] ### 167/206 20mer/4mis, 187/206 15mer/3mis

duplicate_end_new=end_new %>%
  group_by(Chr) %>%
  filter(n()>1)

duplicate_end_new=as.data.frame(duplicate_end_new)


suppressMessages({unique_end_new=anti_join(end_new, duplicate_end_new)})  ### merge to the final start table


unique_end_new=unique_end_new %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length)



skt_end=read.table(input8, sep='\t')

# skt_end=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_POTENTIAL_ENDS.txt", sep='\t')

super_test_end=merge(duplicate_end_new, skt_end, by.x = c("Chr"), by.y =  c("V1"), all.x = TRUE)

super_test_end$Seq=as.character(super_test_end$Seq)
super_test_end$V2=as.character(super_test_end$V2)

ert_end=super_test_end[which(super_test_end$Seq==super_test_end$V2),] ### merge to the final start table
ert_end=as.data.table(ert_end)

new_ert_end=ert_end[ert_end[, .I[which.min(Mismatches)], by=Chr]$V1]

new2_ert_end=new_ert_end %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length)


suppressMessages({skta_end=anti_join(super_test_end, ert_end)})

jk_end=intersect(ert_end$Chr, skta_end$Chr)
jk_end=as.data.frame(jk_end)

clear_skta_end=skta_end[!skta_end$Chr %in% jk_end$jk_end, , drop = FALSE]
clear_skta_end=as.data.table(clear_skta_end)

sure_end_en=clear_skta_end[clear_skta_end[, .I[which.max(Mismatches)], by=Chr]$V1]

# sure_start_st=NULL;
# if (as.character(input5)=="maximum"){
#   sure_start_st=telos_new[passed_start[, .I[which.max(Start)], by=Chr]$V1]
# } else {
#   sure_start_st=telos_new[passed_start[, .I[which.min(Start)], by=Chr]$V1]
# }

sure_end_en=as.data.frame(sure_end_en)### merge to the final start table

sure_end_en=sure_end_en %>%
  select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length)

final_end_table=rbind(unique_end_new, new2_ert_end, sure_end_en)



dfg=setdiff(final_start_table$Chr, final_end_table$Chr)
dfg=as.data.frame(dfg)
colnames(dfg)=c("tautothta")

# dokimh=setdiff(sure_start_st$Chr, sure_end_en$Chr) # 27 elements (20mer/4mis), 16 elements (15mer/3mis), 30 elements (25mer/5mis) that have start and not end
# dokimh
# 
# dokimh1=setdiff(sure_end_en$Chr, sure_start_st$Chr) # 35 elements (20mer/4mis), 31 elements (15mer/3mis), 34 elements (25mer/5mis) that have end and not start
# dokimh1


######################## SUPER MERGING

mega_tab=merge(final_start_table, final_end_table, by.x = c("Chr"), by.y=c("Chr"), all.x=TRUE)
mega_tab=na.omit(mega_tab) 

mega_tab=mega_tab %>%
  select(Chr, Seq.x, Start.x, Seq.y, End.y)

colnames(mega_tab)=c("ID", "Seq_Start", "Start", "Seq_End", "End")

mega_tab=mega_tab %>%
  select(ID, Seq_Start, Seq_End, Start, End)

new_mega_tab=mega_tab %>%
  separate(ID, c("chr", "coords"), sep = ":", remove=FALSE)

new_mega_tab=new_mega_tab %>%
  separate(coords, c("Start_Elem", "End_Elem"), sep="-", remove=FALSE)

new_mega_tab$chr=gsub("rc_", "", new_mega_tab$chr)
new_mega_tab$End_Elem=gsub("_P_element", "", new_mega_tab$End_Elem)


P_STRAND=new_mega_tab[new_mega_tab$ID %like% "_P_element", ]

suppressMessages({D_STRAND=anti_join(new_mega_tab, P_STRAND)})

NEW_D_STRAND=D_STRAND %>%
  mutate(New_Start_Elem=as.numeric(as.character(Start_Elem)) + as.numeric(as.character(Start)),
         New_End_Elem=as.numeric(as.character(Start_Elem)) + as.numeric(as.character(End)))

NEW_P_STRAND=P_STRAND %>%
  mutate(New_End_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(Start)),
         New_Start_Elem=as.numeric(as.character(End_Elem)) - as.numeric(as.character(End)))

NEW_D_STRAND=NEW_D_STRAND %>%
  mutate(Strand="+")

NEW_P_STRAND=NEW_P_STRAND %>%
  mutate(Strand="-")

new2_mega_tab_D=NEW_D_STRAND %>%
  select(chr, New_Start_Elem, New_End_Elem, ID, Seq_Start, Seq_End, Strand)

new2_mega_tab_P=NEW_P_STRAND %>%
  select(chr, New_Start_Elem, New_End_Elem, ID, Seq_Start, Seq_End, Strand)

# dim(new2_mega_tab_D)
# dim(new2_mega_tab_P)

# full_new2_mega_tab=rbind(new2_mega_tab_D, new2_mega_tab_P)

# 20mer element-specific 4 mis --> 132 elements
# 15mer element-specific 3 mis --> 156 elements
# 25mer element-specific 5 mis --> 126 elements

# write.table(new2_mega_tab_D, "t2t-col.20210610.fasta.F2B.D_FULLLENGTH.bed", sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(new2_mega_tab_P, "t2t-col.20210610.fasta.F2B.P_FULLLENGTH.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new2_mega_tab_D, input9, sep = "\t", row.names = F, quote= F, col.names = F)
# 
write.table(new2_mega_tab_P, input10, sep = "\t", row.names = F, quote= F, col.names = F)



###### rejected elements

idi_tab=merge(aporipsh_arxhs, final_start_table, by.x = c("Chr"), by.y=c("Chr"), all.x=TRUE)
id_tab_full=merge(idi_tab, final_end_table, by.x = c("Chr"), by.y=c("Chr"), all.x=TRUE)


final_id_tab_nas=id_tab_full[is.na(id_tab_full$Start.x) | is.na(id_tab_full$Start.y),]
final_id_tab_nas=unique(final_id_tab_nas$Chr)
final_id_tab_nas=as.data.frame(final_id_tab_nas)
colnames(final_id_tab_nas)=c("tautothta")

missing_elements=rbind(final_id_tab_nas, dfg)

rjct_elem=read.table(input11, sep = "\t")

# rjct_elem=read.table("t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_EXTENDED_INTERNAL_CLEAN.bed", sep = "\t")


rjct_elem=unite(rjct_elem, id_half, c(V2, V3), remove=FALSE, sep="-")
rjct_elem=unite(rjct_elem, id_full, c(V1, id_half), remove=FALSE, sep=":")

rjct_elem_DDD=rjct_elem[which(rjct_elem$V17=="+"),]


rjct_elem_PPP=rjct_elem[which(rjct_elem$V17=="-"),]

rjct_elem_PPP$id_full=paste0(rjct_elem_PPP$id_full, "_P_element")
rjct_elem_PPP$id_full=paste("rc", rjct_elem_PPP$id_full, sep="_")


MERGED_rjct_elem=rbind(rjct_elem_DDD, rjct_elem_PPP)

 
new_rjct_elem=merge(missing_elements, MERGED_rjct_elem, by.x=c("tautothta"), by.y=c("id_full"), all.x=TRUE)

new2_rjct_elem=new_rjct_elem %>%
   select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17)

# write.table(new2_rjct_elem, "t2t-col.20210610.fasta.F2B.ALL.PBSNPPT_DP_EXTENDED_INTERNAL_CLEAN_REJECTED.bed", sep = "\t", row.names = F, quote= F, col.names = F)
 
write.table(new2_rjct_elem, input12, sep = "\t", row.names = F, quote= F, col.names = F)



