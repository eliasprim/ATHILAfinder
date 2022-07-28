rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))


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


ourfl=read.table(input1, sep="\t")

# ourfl=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_RENAMED.bed", sep="\t")

ourfl=ourfl %>%
  select(V1, V2, V3, V4, V7)

colnames(ourfl)=c("chr", "genome_left_coord_FL", "genome_right_coord_FL", "TE_ID", "direction")

ourfl=ourfl %>%
  mutate(quality="intact", length_FL=abs(as.numeric(as.character(genome_right_coord_FL)) - as.numeric(as.character(genome_left_coord_FL))))




our5ltr=read.table(input2, sep="\t")

# our5ltr=read.table("t2t-col.20210610.fasta.F2B.DP_5LTR_RENAMED.bed", sep="\t")

colnames(our5ltr)=c("chr", "genome_5LTR_start_coord", "genome_5LTR_stop_coord", "PRIME_ID", "start_mer", "direction")

our5ltr=our5ltr %>%
  mutate(LTR5_length=abs(as.numeric(as.character(genome_5LTR_stop_coord)) - as.numeric(as.character(genome_5LTR_start_coord))))

our5ltr$PRIME_ID=gsub("_5prime", "", our5ltr$PRIME_ID)

our5ltr=our5ltr %>%
  select(genome_5LTR_start_coord, genome_5LTR_stop_coord, PRIME_ID, LTR5_length)



our3ltr=read.table(input3, sep="\t")

# our3ltr=read.table("t2t-col.20210610.fasta.F2B.DP_3LTR_RENAMED.bed", sep="\t")

colnames(our3ltr)=c("chr", "genome_3LTR_start_coord", "genome_3LTR_stop_coord", "PRIME_ID", "end_mer", "direction")

our3ltr=our3ltr %>%
  mutate(LTR3_length=abs(as.numeric(as.character(genome_3LTR_stop_coord)) - as.numeric(as.character(genome_3LTR_start_coord))))

our3ltr$PRIME_ID=gsub("_3prime", "", our3ltr$PRIME_ID)

our3ltr=our3ltr %>%
  select(genome_3LTR_start_coord, genome_3LTR_stop_coord, PRIME_ID, LTR3_length)





fl5ltr=merge(ourfl, our5ltr, by.x = c("TE_ID"), by.y = c("PRIME_ID"), all.x = TRUE)

fl5ltr3ltr=merge(fl5ltr, our3ltr, by.x = c("TE_ID"), by.y = c("PRIME_ID"), all.x = TRUE)

fl5ltr3ltr=fl5ltr3ltr %>%
  select(chr, genome_left_coord_FL, genome_right_coord_FL, TE_ID, direction, quality, length_FL, genome_5LTR_start_coord,
         genome_5LTR_stop_coord, genome_3LTR_start_coord, genome_3LTR_stop_coord, LTR5_length, LTR3_length)

fl5ltr3ltr=fl5ltr3ltr %>%
  mutate(origin="initial_pipeline")




blastfl=read.table(input4, sep="\t")

# blastfl=read.table("t2t-col.20210610.fasta.DP_FULLLENGTH_BLAST_NONOVERLAPPING.bed", sep="\t")

blastfl=blastfl %>%
  select(V1, V2, V3, V4, V7)

colnames(blastfl)=c("chr", "genome_left_coord_FL", "genome_right_coord_FL", "TE_ID", "direction")

blastfl=blastfl %>%
  mutate(quality="intact", length_FL=abs(as.numeric(as.character(genome_right_coord_FL)) - as.numeric(as.character(genome_left_coord_FL))))




blast5ltr=read.table(input5, sep="\t")

# blast5ltr=read.table("t2t-col.20210610.fasta.DP_5LTR_BLAST_NONOVERLAPPING.bed", sep="\t")

blast5ltr=blast5ltr %>%
  select(V1, V2, V3, V4, V7)

colnames(blast5ltr)=c("chr", "genome_5LTR_start_coord", "genome_5LTR_stop_coord", "PRIME_ID", "direction")

blast5ltr=blast5ltr %>%
  mutate(LTR5_length=abs(as.numeric(as.character(genome_5LTR_stop_coord)) - as.numeric(as.character(genome_5LTR_start_coord))))

blast5ltr$PRIME_ID=gsub("_5prime", "", blast5ltr$PRIME_ID)

blast5ltr=blast5ltr %>%
  select(genome_5LTR_start_coord, genome_5LTR_stop_coord, PRIME_ID, LTR5_length)


blast3ltr=read.table(input6, sep="\t")

# blast3ltr=read.table("t2t-col.20210610.fasta.DP_3LTR_BLAST_NONOVERLAPPING.bed", sep="\t")

blast3ltr=blast3ltr %>%
  select(V1, V2, V3, V4, V7)

colnames(blast3ltr)=c("chr", "genome_3LTR_start_coord", "genome_3LTR_stop_coord", "PRIME_ID", "direction")

blast3ltr=blast3ltr %>%
  mutate(LTR3_length=abs(as.numeric(as.character(genome_3LTR_stop_coord)) - as.numeric(as.character(genome_3LTR_start_coord))))

blast3ltr$PRIME_ID=gsub("_3prime", "", blast3ltr$PRIME_ID)

blast3ltr=blast3ltr %>%
  select(genome_3LTR_start_coord, genome_3LTR_stop_coord, PRIME_ID, LTR3_length)



blastfl5ltr=merge(blastfl, blast5ltr, by.x = c("TE_ID"), by.y = c("PRIME_ID"), all.x = TRUE)

blastfl5ltr3ltr=merge(blastfl5ltr, blast3ltr, by.x = c("TE_ID"), by.y = c("PRIME_ID"), all.x = TRUE)

blastfl5ltr3ltr=blastfl5ltr3ltr %>%
  select(chr, genome_left_coord_FL, genome_right_coord_FL, TE_ID, direction, quality, length_FL, genome_5LTR_start_coord,
         genome_5LTR_stop_coord, genome_3LTR_start_coord, genome_3LTR_stop_coord, LTR5_length, LTR3_length)

blastfl5ltr3ltr=blastfl5ltr3ltr %>%
  mutate(origin="BLAST")



allfl=rbind(fl5ltr3ltr, blastfl5ltr3ltr)




pidtab1=read.table(input7, sep="\t", header = T)

# pidtab1=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID1_TABLE.txt", sep="\t", header = T)

allfl_pid1=merge(allfl, pidtab1, by.x = c("TE_ID"), by.y = c("FL_ID"), all.x = TRUE)


pidtab2=read.table(input8, sep="\t", header = T)

# pidtab2=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID2_TABLE.txt", sep="\t", header = T)

allfl_pid2=merge(allfl_pid1, pidtab2, by.x = c("TE_ID"), by.y = c("FL_ID"), all.x = TRUE)


pidtab3=read.table(input9, sep="\t", header = T)

# pidtab3=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID3_TABLE.txt", sep="\t", header = T)

allfl_pid3=merge(allfl_pid2, pidtab3, by.x = c("TE_ID"), by.y = c("FL_ID"), all.x = TRUE)


pidtab4=read.table(input10, sep="\t", header = T)

# pidtab4=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_BLAST_TOGETHER_PID4_TABLE.txt", sep="\t", header = T)

allfl_pid4=merge(allfl_pid3, pidtab4, by.x = c("TE_ID"), by.y = c("FL_ID"), all.x = TRUE)












solo=read.table(input11, sep="\t")

# solo=read.table("t2t-col.20210610.fasta.DP_SOLO_LTR_BLAST_NONOVERLAPPING_CLEAR_PBSPPT.bed", sep="\t")

solo=solo %>%
  select(V1, V2, V3, V4, V7)

colnames(solo)=c("chr", "genome_left_coord_FL", "genome_right_coord_FL", "TE_ID", "direction")

solo=solo %>%
  mutate(quality="solo", length_FL=abs(as.numeric(as.character(genome_right_coord_FL)) - as.numeric(as.character(genome_left_coord_FL))),
         genome_5LTR_start_coord="NA", genome_5LTR_stop_coord="NA", genome_3LTR_start_coord="NA", genome_3LTR_stop_coord="NA",
         LTR5_length="NA", LTR3_length="NA", origin="BLAST", perc_ident1="NA", perc_ident2="NA", perc_ident3="NA", perc_ident4="NA")



allflsolo=rbind(allfl_pid4, solo)




tsd=read.table(input12, sep="\t")

# tsd=read.table("t2t-col.20210610.fasta.F2B.DP_ULTIMATE_TSD_TABLE.txt", sep="\t")

tsd=tsd %>%
  select(V1, V5, V6)

colnames(tsd)=c("TE_ID", "TSD", "TSD_seq")

allflsolotsd=merge(allflsolo, tsd, by.x = c("TE_ID"), by.y = c("TE_ID"), all.x = TRUE)



allflsolotsd=allflsolotsd %>%
  select(chr, genome_left_coord_FL, genome_right_coord_FL, TE_ID, direction, quality, length_FL, TSD, TSD_seq, perc_ident1, perc_ident2,
         perc_ident3, perc_ident4, origin, genome_5LTR_start_coord, genome_5LTR_stop_coord, genome_3LTR_start_coord, 
         genome_3LTR_stop_coord, LTR5_length, LTR3_length)


allflsolotsd=allflsolotsd[order(allflsolotsd$quality),]

allflsolotsd=allflsolotsd %>%
  mutate(species=as.character(input13))


write.table(allflsolotsd, input14, sep = "\t", row.names = F, quote= F)


# write.table(allflsolotsd, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_AND_SOLO_SUMMARY_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)




hmm=read.table(input15, sep="\t", header=T)

# hmm=read.table("t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_SL.fasta.hmmscanned.E001incE001.domtbl_HMMTABLE.txt", sep="\t", header=T)

hmm=hmm %>%
  select(-contains("FROM"), -contains("DIST"), -contains("RATIO"))

names(hmm)=gsub(x=names(hmm), pattern = "_envTO", replacement = "")

hmm[is.na(hmm)]=0

hmm=hmm %>%
  mutate_at(vars(-1), ~ str_replace(., "^0$", "N"))
            
hmm=hmm %>%
  mutate_at(vars(-1), ~ str_replace(., "[^N]+", "Y"))

hmm$query_name=gsub("_internal", "", hmm$query_name)

intact_fl=allflsolotsd[which(allflsolotsd$quality=="intact"),]

intact_fl=intact_fl %>%
  mutate(species=as.character(input13))

allflhmm=merge(intact_fl, hmm, by.x = c("TE_ID"), by.y = c("query_name"), all.x = TRUE)



write.table(allflhmm, input16, sep = "\t", row.names = F, quote= F)


# write.table(allflhmm, "t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_AND_HMM_SUMMARY_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)


