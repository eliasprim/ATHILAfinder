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


############ OUR FULLLENGTH ELEMENTS

d_v2_left=read.table(input1, sep='\t')

# d_v2_left=read.table("Evs-0.ragtag_scaffolds.fa.F2B.D_FULLLENGTH_NEW_PLUS_TSD_LEFTV2.txt", sep='\t')

d_v3_right=read.table(input2, sep='\t')

# d_v3_right=read.table("Evs-0.ragtag_scaffolds.fa.F2B.D_FULLLENGTH_NEW_PLUS_TSD_RIGHTV3.txt", sep='\t')

p_v2_left=read.table(input3, sep='\t')

# p_v2_left=read.table("Evs-0.ragtag_scaffolds.fa.F2B.P_FULLLENGTH_NEW_PLUS_TSD_LEFTV2.fasta_sl.txt", sep='\t')

p_v3_right=read.table(input4, sep='\t')

# p_v3_right=read.table("Evs-0.ragtag_scaffolds.fa.F2B.P_FULLLENGTH_NEW_PLUS_TSD_RIGHTV3.fasta_sl.txt", sep='\t')


d_elem=merge(d_v2_left, d_v3_right, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
colnames(d_elem)=c("fl_id", "tsd_left", "tsd_right")


p_elem=merge(p_v2_left, p_v3_right, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
colnames(p_elem)=c("fl_id", "tsd_left", "tsd_right")
p_elem$fl_id=gsub("_P_element", "", p_elem$fl_id)
p_elem$fl_id=gsub("rc_", "", p_elem$fl_id)


dp_elem=rbind(d_elem, p_elem)


mismatches=adist(dp_elem$tsd_left, dp_elem$tsd_right)
mismatches=diag(mismatches)
mismatches=as.data.frame(mismatches)

mis_dp_elem=cbind(dp_elem, mismatches)


same_tsd=mis_dp_elem[which(mis_dp_elem$mismatches==0 | mis_dp_elem$mismatches==1),]

same_tsd=same_tsd %>%
  mutate(same_tsd="Y", tsd_seq=tsd_left)


not_same_tsd=mis_dp_elem[which(mis_dp_elem$mismatches!=0 & mis_dp_elem$mismatches!=1),]

not_same_tsd=not_same_tsd %>%
  mutate(same_tsd="N", tsd_seq=tsd_right)


ult_tsd_tab=rbind(same_tsd, not_same_tsd)


write.table(ult_tsd_tab, input5, sep = "\t", row.names = F, quote= F, col.names = F)





############ BLAST FULLLENGTH ELEMENTS

# d_v2_left_bl=read.table(input6, sep='\t')

# d_v2_left_bl=read.table("Evs-0.ragtag_scaffolds.fa.D_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.txt", sep='\t')


mtry6=try(read.table(input6, sep = "\t"), 
          silent = TRUE)

if (class(mtry6) != "try-error") {
  d_v2_left_bl=read.table(input6, sep = "\t")
} else {
  message("File is empty, please check")
  d_v2_left_bl=data.frame()
}


# d_v3_right_bl=read.table(input7, sep='\t')

# d_v3_right_bl=read.table("Evs-0.ragtag_scaffolds.fa.D_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.txt", sep='\t')


mtry7=try(read.table(input7, sep = "\t"), 
          silent = TRUE)

if (class(mtry7) != "try-error") {
  d_v3_right_bl=read.table(input7, sep = "\t")
} else {
  message("File is empty, please check")
  d_v3_right_bl=data.frame()
}

# p_v2_left_bl=read.table(input8, sep='\t')

# p_v2_left_bl=read.table("Evs-0.ragtag_scaffolds.fa.P_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.fasta_sl.txt", sep='\t')


mtry8=try(read.table(input8, sep = "\t"), 
            silent = TRUE)

if (class(mtry8) != "try-error") {
  p_v2_left_bl=read.table(input8, sep = "\t")
} else {
  message("File is empty, please check")
  p_v2_left_bl=data.frame()
}


# p_v3_right_bl=read.table(input9, sep='\t')

# p_v3_right_bl=read.table("Evs-0.ragtag_scaffolds.fa.P_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.fasta_sl.txt", sep='\t')


mtry9=try(read.table(input9, sep = "\t"), 
          silent = TRUE)

if (class(mtry9) != "try-error") {
  p_v3_right_bl=read.table(input9, sep = "\t")
} else {
  message("File is empty, please check")
  p_v3_right_bl=data.frame()
}




if (dim(p_v2_left_bl)[1] == 0 && dim(p_v3_right_bl)[1] == 0) {
  
  d_elem_bl=merge(d_v2_left_bl, d_v3_right_bl, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
  colnames(d_elem_bl)=c("fl_id", "tsd_left", "tsd_right")
  
  # dp_elem_bl=rbind(d_elem_bl, p_elem_bl)
  
  dp_elem_bl=d_elem_bl
  
  
  mismatches_bl=adist(dp_elem_bl$tsd_left, dp_elem_bl$tsd_right)
  mismatches_bl=diag(mismatches_bl)
  mismatches_bl=as.data.frame(mismatches_bl)
  
  mis_dp_elem_bl=cbind(dp_elem_bl, mismatches_bl)
  
  
  same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches==0 | mis_dp_elem_bl$mismatches==1),]
  
  same_tsd_bl=same_tsd_bl %>%
    mutate(same_tsd="Y", tsd_seq=tsd_left)
  
  
  not_same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches!=0 & mis_dp_elem_bl$mismatches!=1),]
  
  not_same_tsd_bl=not_same_tsd_bl %>%
    mutate(same_tsd="N", tsd_seq=tsd_right)
  
  
  ult_tsd_tab_bl=rbind(same_tsd_bl, not_same_tsd_bl)
  
} else {
  
  d_elem_bl=merge(d_v2_left_bl, d_v3_right_bl, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
  colnames(d_elem_bl)=c("fl_id", "tsd_left", "tsd_right")
  
  
  p_elem_bl=merge(p_v2_left_bl, p_v3_right_bl, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
  colnames(p_elem_bl)=c("fl_id", "tsd_left", "tsd_right")
  p_elem_bl$fl_id=gsub("_P_element", "", p_elem_bl$fl_id)
  p_elem_bl$fl_id=gsub("rc_", "", p_elem_bl$fl_id)
  
  
  dp_elem_bl=rbind(d_elem_bl, p_elem_bl)
  
  
  mismatches_bl=adist(dp_elem_bl$tsd_left, dp_elem_bl$tsd_right)
  mismatches_bl=diag(mismatches_bl)
  mismatches_bl=as.data.frame(mismatches_bl)
  
  mis_dp_elem_bl=cbind(dp_elem_bl, mismatches_bl)
  
  
  same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches==0 | mis_dp_elem_bl$mismatches==1),]
  
  same_tsd_bl=same_tsd_bl %>%
    mutate(same_tsd="Y", tsd_seq=tsd_left)
  
  
  not_same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches!=0 & mis_dp_elem_bl$mismatches!=1),]
  
  not_same_tsd_bl=not_same_tsd_bl %>%
    mutate(same_tsd="N", tsd_seq=tsd_right)
  
  
  ult_tsd_tab_bl=rbind(same_tsd_bl, not_same_tsd_bl)
  
}


# d_elem_bl=merge(d_v2_left_bl, d_v3_right_bl, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
# colnames(d_elem_bl)=c("fl_id", "tsd_left", "tsd_right")
# 
# 
# p_elem_bl=merge(p_v2_left_bl, p_v3_right_bl, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
# colnames(p_elem_bl)=c("fl_id", "tsd_left", "tsd_right")
# p_elem_bl$fl_id=gsub("_P_element", "", p_elem_bl$fl_id)
# p_elem_bl$fl_id=gsub("rc_", "", p_elem_bl$fl_id)
# 
# 
# dp_elem_bl=rbind(d_elem_bl, p_elem_bl)
# 
# 
# mismatches_bl=adist(dp_elem_bl$tsd_left, dp_elem_bl$tsd_right)
# mismatches_bl=diag(mismatches_bl)
# mismatches_bl=as.data.frame(mismatches_bl)
# 
# mis_dp_elem_bl=cbind(dp_elem_bl, mismatches_bl)
# 
# 
# same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches==0 | mis_dp_elem_bl$mismatches==1),]
# 
# same_tsd_bl=same_tsd_bl %>%
#   mutate(same_tsd="Y", tsd_seq=tsd_left)
# 
# 
# not_same_tsd_bl=mis_dp_elem_bl[which(mis_dp_elem_bl$mismatches!=0 & mis_dp_elem_bl$mismatches!=1),]
# 
# not_same_tsd_bl=not_same_tsd_bl %>%
#   mutate(same_tsd="N", tsd_seq=tsd_right)
# 
# 
# ult_tsd_tab_bl=rbind(same_tsd_bl, not_same_tsd_bl)


write.table(ult_tsd_tab_bl, input10, sep = "\t", row.names = F, quote= F, col.names = F)





############ BLAST SOLO LTRs 

d_v2_left_solo=read.table(input11, sep='\t')

# d_v2_left_solo=read.table("Evs-0.ragtag_scaffolds.fa.D_SOLO_LTR_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.txt", sep='\t')

d_v3_right_solo=read.table(input12, sep='\t')

# d_v3_right_solo=read.table("Evs-0.ragtag_scaffolds.fa.D_SOLO_LTR_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.txt", sep='\t')

p_v2_left_solo=read.table(input13, sep='\t')

# p_v2_left_solo=read.table("Evs-0.ragtag_scaffolds.fa.P_SOLO_LTR_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.fasta_sl.txt", sep='\t')

p_v3_right_solo=read.table(input14, sep='\t')

# p_v3_right_solo=read.table("Evs-0.ragtag_scaffolds.fa.P_SOLO_LTR_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.fasta_sl.txt", sep='\t')


d_elem_solo=merge(d_v2_left_solo, d_v3_right_solo, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
colnames(d_elem_solo)=c("fl_id", "tsd_left", "tsd_right")


p_elem_solo=merge(p_v2_left_solo, p_v3_right_solo, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)
colnames(p_elem_solo)=c("fl_id", "tsd_left", "tsd_right")
p_elem_solo$fl_id=gsub("_P_element", "", p_elem_solo$fl_id)
p_elem_solo$fl_id=gsub("rc_", "", p_elem_solo$fl_id)


dp_elem_solo=rbind(d_elem_solo, p_elem_solo)


mismatches_solo=adist(dp_elem_solo$tsd_left, dp_elem_solo$tsd_right)
mismatches_solo=diag(mismatches_solo)
mismatches_solo=as.data.frame(mismatches_solo)

mis_dp_elem_solo=cbind(dp_elem_solo, mismatches_solo)


same_tsd_solo=mis_dp_elem_solo[which(mis_dp_elem_solo$mismatches==0 | mis_dp_elem_solo$mismatches==1),]

same_tsd_solo=same_tsd_solo %>%
  mutate(same_tsd="Y", tsd_seq=tsd_left)


not_same_tsd_solo=mis_dp_elem_solo[which(mis_dp_elem_solo$mismatches!=0 & mis_dp_elem_solo$mismatches!=1),]

not_same_tsd_solo=not_same_tsd_solo %>%
  mutate(same_tsd="N", tsd_seq=tsd_right)


ult_tsd_tab_solo=rbind(same_tsd_solo, not_same_tsd_solo)


write.table(ult_tsd_tab_solo, input15, sep = "\t", row.names = F, quote= F, col.names = F)


clear_tsd_tab_solo=ult_tsd_tab_solo[which(ult_tsd_tab_solo$same_tsd=="Y"),]


write.table(clear_tsd_tab_solo, input16, sep = "\t", row.names = F, quote= F, col.names = F)


no_tsd_tab_solo=ult_tsd_tab_solo[which(ult_tsd_tab_solo$same_tsd=="N"),]


write.table(no_tsd_tab_solo, input17, sep = "\t", row.names = F, quote= F, col.names = F)