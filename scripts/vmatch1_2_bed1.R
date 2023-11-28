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

# STEP ARGUMENT FOR BASH SCRIPT STEP 1

############ EDIT VMATCH TABLES AND PREPARE BED INPUT FOR BEDTOOLS

# vmatch_tab=read.table(input1, fill = T)

# vmatch_tab=read.table("t2t-col.20210610.fasta.F2B.PPT.vmatch", fill=T)

# flt1=vmatch_tab[which(vmatch_tab$V1!="Query:" & vmatch_tab$V1!="!" & vmatch_tab$V1!="!!"),]
# 
# flt2_start=flt1[which(flt1$V1=="21" | flt1$V1=="22"),]
# 
# colnames(flt2_start)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")
# 
# flt3_end=flt1[which(flt1$V1=="Sbjct:"),]
# 
# colnames(flt3_end)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")
# 
# flt4=cbind(flt2_start, flt3_end)
# 
# final_flt=flt4 %>%
#   select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)
# 
# final_flt$Mismatches=gsub("-", "", final_flt$Mismatches)
# 
# colnames(final_flt)[8]="Length"
# 
# final_flt_d=final_flt[which(final_flt$Strand=="D"),]
# 
# final_flt_d=final_flt_d %>%
#   mutate(RevComp_Seq=Seq)
# 
# final_flt_p=final_flt[which(final_flt$Strand=="P"),]


# ####################### REVERSE COMPLEMENT FUNCTIONS

seq_rev <- function(char) {
  alphabets <- strsplit(char, split = "")[[1]]
  return(rev(alphabets))
}

seq_compl <- function(seq) {
  # Check if there's "T" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  cmplvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A", "N" = "N")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A", "N" = "N")
    }
  })
  return(paste(cmplvec, collapse = ""))
}

revcom <- function(input) {
  # Make sure the input is character and in upper case
  input <- as.character(input)
  input <- toupper(input)
  # Use regular expression to check if there's only legal bases
  # present in the sequence
  legal_char <- Reduce(`&`, grepl("^[A,T,C,G,U,N]*$", input))
  if (!legal_char) {
    stop("revcom() only applies to DNA/RNA sequences, and only A/T/C/G/U/N is allowed")
  }
  rev <- seq_rev(input)
  return(seq_compl(rev))
}



rtry=try(read.table(input1, fill = T), 
          silent = TRUE)

if (class(rtry) != "try-error") {
  
  vmatch_tab=read.table(input1, fill = T)
  
  flt1=vmatch_tab[which(vmatch_tab$V1!="Query:" & vmatch_tab$V1!="!" & vmatch_tab$V1!="!!"),]
  
  flt2_start=flt1[which(flt1$V1=="20" | flt1$V1=="21"),]
  
  colnames(flt2_start)=c("Length_1", "Chr", "Start", "Strand", "Length_2", "Type", "Zeros", "Mismatches")
  
  flt3_end=flt1[which(flt1$V1=="Sbjct:"),]
  
  colnames(flt3_end)=c("Sbjct", "Seq", "End", "Chaos_1", "NA_1", "Chaos_2", "NA_2", "NA_3")
  
  flt4=cbind(flt2_start, flt3_end)
  
  final_flt=flt4 %>%
    select(Chr, Seq, Type, Mismatches, Start, End, Strand, Length_1)
  
  final_flt$Mismatches=gsub("-", "", final_flt$Mismatches)
  
  colnames(final_flt)[8]="Length"
  
  final_flt_d=final_flt[which(final_flt$Strand=="D"),]
  
  final_flt_d=final_flt_d %>%
    mutate(RevComp_Seq=Seq)
  
  final_flt_p=final_flt[which(final_flt$Strand=="P"),]
  
  y=NULL;
  for (i in final_flt_p$Seq)
  {
    tmp=revcom(i)
    y=rbind(y, tmp)
  }
  
  new_final_flt_p=cbind(final_flt_p, y)
  
  colnames(new_final_flt_p)[9]="RevComp_Seq"
  
  ult_flt=rbind(final_flt_d, new_final_flt_p)
  
  new_ult_flt=ult_flt %>%
    select(Chr, Start, End, Seq, RevComp_Seq, Type, Mismatches, Length, Strand)
  
  new_ult_flt$Type=gsub("junction", "", new_ult_flt$Type)
  # new_ult_flt$Chr=gsub("C", "c", new_ult_flt$Chr)
  
  new_ult_flt=new_ult_flt %>%
    mutate(rn=row_number())
  
  new_ult_flt=unite(new_ult_flt, id_0, c(Type, Chr, Mismatches, Start), remove=FALSE, sep="_")
  new_ult_flt=unite(new_ult_flt, id_1, c(End, Strand, rn), remove=FALSE, sep="_")
  new_ult_flt=unite(new_ult_flt, ID, c(id_0, id_1), remove=TRUE, sep=".")
  
  new2_ult_flt=new_ult_flt %>%
    select(Chr, Start, End, ID, Seq, RevComp_Seq, Type, Mismatches, Length, Strand)
  
  new2_ult_flt$Strand=gsub("D", "+", new2_ult_flt$Strand)
  new2_ult_flt$Strand=gsub("P", "-", new2_ult_flt$Strand)
  
  new2_ult_flt_d=new2_ult_flt[which(new2_ult_flt$Strand=="+"),]
  new2_ult_flt_p=new2_ult_flt[which(new2_ult_flt$Strand=="-"),]
  
  write.table(new2_ult_flt_d, input2, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(new2_ult_flt_p, input3, sep = "\t", row.names = F, quote= F, col.names = F)
  
} else {
  message("File is empty, please check")
  
  vmatch_tab_d=data.frame()
  vmatch_tab_p=data.frame()
  
  write.table(vmatch_tab_d, input2, sep = "\t", row.names = F, quote= F, col.names = F)
  write.table(vmatch_tab_p, input3, sep = "\t", row.names = F, quote= F, col.names = F)
  
}



# new_final_flt_p=cbind(final_flt_p, y)
# 
# colnames(new_final_flt_p)[9]="RevComp_Seq"
# 
# ult_flt=rbind(final_flt_d, new_final_flt_p)
# 
# new_ult_flt=ult_flt %>%
#   select(Chr, Start, End, Seq, RevComp_Seq, Type, Mismatches, Length, Strand)
# 
# new_ult_flt$Type=gsub("junction", "", new_ult_flt$Type)
# # new_ult_flt$Chr=gsub("C", "c", new_ult_flt$Chr)
# 
# new_ult_flt=new_ult_flt %>%
#   mutate(rn=row_number())
# 
# new_ult_flt=unite(new_ult_flt, id_0, c(Type, Chr, Mismatches, Start), remove=FALSE, sep="_")
# new_ult_flt=unite(new_ult_flt, id_1, c(End, Strand, rn), remove=FALSE, sep="_")
# new_ult_flt=unite(new_ult_flt, ID, c(id_0, id_1), remove=TRUE, sep=".")
# 
# new2_ult_flt=new_ult_flt %>%
#   select(Chr, Start, End, ID, Seq, RevComp_Seq, Type, Mismatches, Length, Strand)
# 
# new2_ult_flt$Strand=gsub("D", "+", new2_ult_flt$Strand)
# new2_ult_flt$Strand=gsub("P", "-", new2_ult_flt$Strand)
# 
# new2_ult_flt_d=new2_ult_flt[which(new2_ult_flt$Strand=="+"),]
# new2_ult_flt_p=new2_ult_flt[which(new2_ult_flt$Strand=="-"),]
# 
# write.table(new2_ult_flt_d, input2, sep = "\t", row.names = F, quote= F, col.names = F)
# write.table(new2_ult_flt_p, input3, sep = "\t", row.names = F, quote= F, col.names = F)

