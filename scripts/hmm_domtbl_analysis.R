rm(list=ls())

options(max.print=10000000)
options(scipen = 999)


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(tibble))

# setwd('~/Desktop/')


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


domtbl=read.table(input1, sep="", fill = T, header=F, row.names=NULL)

# domtbl=read.table("Aalpi_genome.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_SL.fasta.hmmscanned.E001incE001.domtbl", sep="", fill=T, row.names = NULL, header = F)

domtbl=domtbl[complete.cases(domtbl),]

domtbl=domtbl %>%
  select(V1:V22)

colnames(domtbl)=c("target_name", "accession", "tlen", "query_name", "accession_minus", "qlen", "FS_evalue",
                   "FS_score", "FS_bias", "TD_number", "TD_of", "TD_c_evalue", "TD_i_evalue", "TD_score", "TD_bias",
                   "hmm_coords_from", "hmm_coords_to", "alignment_coords_from", "alignment_coords_to", "env_coords_from",
                   "env_coords_to", "acc")

# domtbl=domtbl[which(domtbl$target_name!="Sirevirus_ENVmonocot-mod.stock" & domtbl$target_name!="Sirevirus_ENVdicot-mod.stock"),]


keep_AA=domtbl %>%
  select(query_name, env_coords_from, env_coords_to, target_name, accession)

keep_AA=keep_AA %>%
  mutate(neo_id=query_name)

keep_AA$neo_id=gsub("_[^_]+$", "", keep_AA$neo_id)
keep_AA$neo_id=gsub("_internal", "", keep_AA$neo_id)

keep_AA_new=unite(keep_AA, hmm_id, c(neo_id, target_name), remove=FALSE, sep="~~~~~~")

keep_AA_new=keep_AA_new %>%
  select(query_name, env_coords_from, env_coords_to, hmm_id, target_name, accession)

keep_AA_new=keep_AA_new %>%
  mutate(strand="+")

write.table(keep_AA_new, input2, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(keep_AA_new, "t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_sl.fasta.hmmscanned.E001incE001.domtbl_HMM_AA.bed", sep = "\t", row.names = F, quote= F, col.names = F)


aa_2_nn=domtbl %>%
  mutate(# hmm_coords_from_NN=(as.numeric(as.character(hmm_coords_from)) * 3) - 2,
         # hmm_coords_to_NN=(as.numeric(as.character(hmm_coords_to)) * 3) - 2,
         alignment_coords_from_NN=(as.numeric(as.character(alignment_coords_from)) * 3) - 1,
         alignment_coords_to_NN=(as.numeric(as.character(alignment_coords_to)) * 3) - 1,
         env_coords_from_NN=(as.numeric(as.character(env_coords_from)) * 3) - 1,
         env_coords_to_NN=(as.numeric(as.character(env_coords_to)) * 3) - 1,
         qlen_NN=(as.numeric(as.character(qlen)) * 3) - 1)


orf300_ids=read.table(input3, fill = T)

# orf300_ids=read.table("Aalpi_genome.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_IDs.txt", sep = "")

orf300_ids=orf300_ids %>%
  select(V1, V2, V4)

orf300_ids$V2=gsub("[[]", "", orf300_ids$V2)
orf300_ids$V4=gsub("[]]", "", orf300_ids$V4)

colnames(orf300_ids)=c("Internal_ID", "ORF_Start", "ORF_End")

new_domtbl_orf300=merge(aa_2_nn, orf300_ids, by.x=c("query_name"), by.y=c("Internal_ID"), all.x=TRUE)

right_coords_domtbl_orf300=new_domtbl_orf300 %>%
  mutate(# corr_coord_hmm_coords_from_NN=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(hmm_coords_from_NN)),
         # corr_coord_hmm_coords_to_NN=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(hmm_coords_to_NN)),
         alignFROM=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(alignment_coords_from_NN)),
         alignTO=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(alignment_coords_to_NN)),
         envAAAAA=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(env_coords_from_NN)),
         envBBBBB=as.numeric(as.character(ORF_Start)) + as.numeric(as.character(env_coords_to_NN)))

updated_coords_domtbl_orf300=right_coords_domtbl_orf300 %>%
  select(target_name, accession, tlen, query_name, qlen, qlen_NN, FS_evalue, FS_score, FS_bias, TD_c_evalue,
         TD_i_evalue, TD_score, TD_bias, hmm_coords_from, hmm_coords_to, ORF_Start, ORF_End, 
         alignFROM, alignTO, envAAAAA, envBBBBB)


query_name=unique(updated_coords_domtbl_orf300$query_name)
query_name=as.data.frame(query_name)

for (i in query_name)
{
  updated_coords_sorted=updated_coords_domtbl_orf300[order(updated_coords_domtbl_orf300$envAAAAA),]
}

dokimh_final=updated_coords_sorted %>%
  select(target_name, query_name, qlen_NN, ORF_Start, ORF_End,
         alignFROM, alignTO, envAAAAA, envBBBBB)


###### run the following command if you want the results for the whole internal and not just the ORFs
dokimh_final$query_name=sub("_[^_]+$", "", dokimh_final$query_name)

dokimh_final$envAAAAA=as.numeric(as.character(dokimh_final$envAAAAA))
dokimh_final$envBBBBB=as.numeric(as.character(dokimh_final$envBBBBB))

dokimh_final=dokimh_final %>%
  mutate(envCCCCC=as.numeric(as.character(envBBBBB)) - as.numeric(as.character(envAAAAA)))

dokimh_final_4_table=dokimh_final



############################ SCRIPT PART TO GET THE NT AND AA SEQUENCES BY USING THE BED FILES BELOW

hmm_aa_length=read.table(input4, sep="\t")

# hmm_aa_length=read.table("orfis.gypsy.updated.hmmdb_NAME_AND_LENGTH.txt", sep="\t")

hmm_aa_length=hmm_aa_length %>%
  mutate(nucl_length=as.numeric(as.character(V2)) * 3)

hmm_aa_length=hmm_aa_length %>%
  select(-V2)

dokimh_final=merge(dokimh_final, hmm_aa_length, by.x = c("target_name"), by.y =  c("V1"), all.x = TRUE)

dokimh_final=dokimh_final %>%
  mutate(envDDDDD=as.numeric(as.character(envCCCCC)) / as.numeric(as.character(nucl_length)))

dokimh_final$envDDDDD=round(dokimh_final$envDDDDD, digits = 3)


count_targets=as.data.frame(table(dokimh_final$target_name))
count_targets=count_targets[order(count_targets$Freq, decreasing = TRUE),]

write.table(count_targets, input5, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(count_targets, "t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_sl.fasta.hmmscanned.E001incE001.domtbl_COUNT_TABLE.txt", sep = "\t", row.names = F, quote= F, col.names = F)



mhkos_5ltr=read.table(input6, sep="\t")

# mhkos_5ltr=read.table("NT1_030222.fasta.F2B.DP_5LTR_BLAST_TOGETHER_LENGTH.txt", sep="\t")

mhkos_5ltr$V1=gsub("_5prime", "", mhkos_5ltr$V1)

dokimh_final=dokimh_final %>%
  mutate(ID_2_destroy=query_name)

dokimh_final$ID_2_destroy=gsub(as.character(input7), "", dokimh_final$ID_2_destroy)

# dokimh_final$ID_2_destroy=gsub(".*Chr", "", dokimh_final$ID_2_destroy)

# dokimh_final$ID_2_destroy=paste0("Chr", dokimh_final$ID_2_destroy)

dokimh_final_new=dokimh_final %>%
  separate(ID_2_destroy, c("part1", "part2"), sep = "-", remove=FALSE)

dokimh_final_new$query_name=gsub("_internal", "", dokimh_final_new$query_name)

dokimh_final_new=merge(dokimh_final_new, mhkos_5ltr, by.x = c("query_name"), by.y =  c("V1"), all.x = TRUE)

### dokimh_final_new$part1=gsub(".*_", "", dokimh_final_new$part1)
dokimh_final_new$part2=gsub("_Athila_internal", "", dokimh_final_new$part2)

dokimh_final_new=dokimh_final_new %>%
  separate(part1, c("chr", "arxh"), sep = "\\.", remove=TRUE)

dokimh_final_new=dokimh_final_new %>%
  separate(part2, c("telos", "strand"), sep = "_", remove=TRUE)

dokimh_final_new$strand=gsub("D", "+", dokimh_final_new$strand)
dokimh_final_new$strand=gsub("P", "-", dokimh_final_new$strand)

road_2_bed=unite(dokimh_final_new, query_and_target, c(query_name, target_name), remove=TRUE, sep="~~~~~~")

road_2_bed=road_2_bed %>%
  select(-qlen_NN, -ORF_Start, -ORF_End, -alignFROM, -alignTO)

road_2_bed_D=road_2_bed[which(road_2_bed$strand=="+"),]
road_2_bed_P=road_2_bed[which(road_2_bed$strand=="-"),]

road_2_bed_D=road_2_bed_D %>%
  mutate(NEA_ARXH=as.numeric(as.character(arxh)) + as.numeric(as.character(V2)) + as.numeric(as.character(envAAAAA)),
         NEO_TELOS=as.numeric(as.character(arxh)) + as.numeric(as.character(V2)) + as.numeric(as.character(envBBBBB)))

road_2_bed_D=road_2_bed_D %>%
  select(chr, NEA_ARXH, NEO_TELOS, query_and_target, envCCCCC, strand)


write.table(road_2_bed_D, input8, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(road_2_bed_D, "t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_sl.fasta.hmmscanned.E001incE001.domtbl_D_HMM_BED.bed", sep = "\t", row.names = F, quote= F, col.names = F)


road_2_bed_P=road_2_bed_P %>%
  mutate(NEO_TELOS=as.numeric(as.character(telos)) - as.numeric(as.character(V2)) - as.numeric(as.character(envAAAAA)),
         NEA_ARXH=as.numeric(as.character(telos)) - as.numeric(as.character(V2)) - as.numeric(as.character(envBBBBB)))

road_2_bed_P=road_2_bed_P %>%
  select(chr, NEA_ARXH, NEO_TELOS, query_and_target, envCCCCC, strand)


write.table(road_2_bed_P, input9, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(road_2_bed_P, "t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_sl.fasta.hmmscanned.E001incE001.domtbl_P_HMM_BED.bed", sep = "\t", row.names = F, quote= F, col.names = F)





############################ SCRIPT PART TO MAKE THE FULL HMM TABLE

multiple_minimum=dokimh_final_4_table %>%
  distinct(target_name, query_name, .keep_all = TRUE)

multiple_minimum=unite(multiple_minimum, target_query, c(target_name, query_name), remove=FALSE, sep="~~~~~~")

multiple_minimum=multiple_minimum %>%
  select(-target_name, -query_name, -qlen_NN, -ORF_Start, -ORF_End, -alignFROM, -alignTO, -envBBBBB)

multiple_maximum=dokimh_final[order(-dokimh_final$envBBBBB),]

multiple_maximum=multiple_maximum %>%
  distinct(target_name, query_name, .keep_all = TRUE)

multiple_maximum=multiple_maximum[order(multiple_maximum$envAAAAA),]

multiple_maximum=unite(multiple_maximum, target_query, c(target_name, query_name), remove=FALSE, sep="~~~~~~")

multiple_maximum=multiple_maximum %>%
  select(-envAAAAA)

multiple_all=merge(multiple_minimum, multiple_maximum, by.x = c("target_query"), by.y =  c("target_query"))

multiple_all=multiple_all %>%
  select(target_name, query_name, qlen_NN, ORF_Start, ORF_End, alignFROM, alignTO, envAAAAA, envBBBBB)

multiple_all=multiple_all[order(multiple_all$envAAAAA),]

multiple_all=multiple_all %>%
  mutate(envCCCCC=as.numeric(as.character(envBBBBB)) - as.numeric(as.character(envAAAAA)))

multiple_all=merge(multiple_all, hmm_aa_length, by.x = c("target_name"), by.y =  c("V1"), all.x = TRUE)

multiple_all=multiple_all %>%
  mutate(envDDDDD=as.numeric(as.character(envCCCCC)) / as.numeric(as.character(nucl_length)))

multiple_all$envDDDDD=round(multiple_all$envDDDDD, digits = 3)





############################ MAKE THE FULL HMM TABLE -- FINAL PART

df1=multiple_all %>%
  mutate(target_type = paste0(target_name, "target_type")) %>%
  pivot_wider(id_cols = query_name, names_from = target_name, names_glue = "{target_name}_{.value}", names_sort = TRUE, 
              values_from=c(envAAAAA, envBBBBB, envCCCCC, envDDDDD)) %>%
  left_join(multiple_all, by = c("query_name" = "query_name"))

df1=as.data.frame(df1)


clean_df1=df1 %>%
  distinct(query_name, .keep_all = TRUE) %>%
  select(-target_name, -alignFROM, -alignTO, -envAAAAA, -envBBBBB, -qlen_NN, -ORF_Start, -ORF_End, -envCCCCC, -nucl_length, -envDDDDD)


orf300_ids$Internal_ID=gsub("_[^_]+$", "", orf300_ids$Internal_ID)
unique_orf300_ids=unique(orf300_ids$Internal_ID)
unique_orf300_ids=as.data.frame(unique_orf300_ids)


rest_elem=setdiff(unique_orf300_ids$unique_orf300_ids, clean_df1$query_name)
rest_elem=as.data.frame(rest_elem)
colnames(rest_elem)=c("query_name")


diff_hmms=colnames(clean_df1)
diff_hmms=diff_hmms[-1]


if (dim(rest_elem)[1] == 0) {
  new_rest_elem=clean_df1 
} else {
  new_rest_elem=cbind(rest_elem, setNames(lapply(diff_hmms, function(x) x=0), diff_hmms))
}


# new_rest_elem=cbind(rest_elem, setNames(lapply(diff_hmms, function(x) x=0), diff_hmms))


new_clean_df1=rbind(clean_df1, new_rest_elem)

colnames(new_clean_df1)[1]="111111111111"

sorted_df1=new_clean_df1[ ,order(names(new_clean_df1))]

# names(sorted_df1)=gsub(x = names(sorted_df1), pattern = "AAAAA", replacement = "FROM")
# names(sorted_df1)=gsub(x = names(sorted_df1), pattern = "BBBBB", replacement = "TO")
# names(sorted_df1)=gsub(x = names(sorted_df1), pattern = "CCCCC", replacement = "DIST")
# names(sorted_df1)=gsub(x = names(sorted_df1), pattern = "DDDDD", replacement = "RATIO")


rest_hmm=setdiff(hmm_aa_length$V1, multiple_all$target_name)



#### new code

rest_hmmFROM=paste0(rest_hmm, "_envAAAAA")
rest_hmmTO=paste0(rest_hmm, "_envBBBBB")
rest_hmmDIST=paste0(rest_hmm, "_envCCCCC")
rest_hmmRATIO=paste0(rest_hmm, "_envDDDDD")


final_hmm_table1=cbind(sorted_df1, setNames(lapply(rest_hmmFROM, function(x) x=0), rest_hmmFROM))
final_hmm_table2=cbind(final_hmm_table1, setNames(lapply(rest_hmmTO, function(x) x=0), rest_hmmTO))
final_hmm_table3=cbind(final_hmm_table2, setNames(lapply(rest_hmmDIST, function(x) x=0), rest_hmmDIST))
final_hmm_table4=cbind(final_hmm_table3, setNames(lapply(rest_hmmRATIO, function(x) x=0), rest_hmmRATIO))


final_hmm_table5=final_hmm_table4[,order(colnames(final_hmm_table4))]


names(final_hmm_table5)=gsub(x = names(final_hmm_table5), pattern = "AAAAA", replacement = "FROM")
names(final_hmm_table5)=gsub(x = names(final_hmm_table5), pattern = "BBBBB", replacement = "TO")
names(final_hmm_table5)=gsub(x = names(final_hmm_table5), pattern = "CCCCC", replacement = "DIST")
names(final_hmm_table5)=gsub(x = names(final_hmm_table5), pattern = "DDDDD", replacement = "RATIO")

### end of new code

colnames(final_hmm_table5)[1]="query_name"



write.table(final_hmm_table5, input10, sep = "\t", row.names = F, quote= F)


# write.table(final_hmm_table, "t2t-col.20210610.fasta.F2B.DP_INTERNAL_BLAST_TOGETHER.fasta.ORF300_sl.fasta.hmmscanned.E001incE001.domtbl_HMMTABLE.txt", sep = "\t", row.names = F, quote= F, col.names = T)

