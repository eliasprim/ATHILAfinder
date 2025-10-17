rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("plyranges")

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(plyranges))


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
input11 = args[11]
input12 = args[12]


# our server 59565
# our server remove -F 59561
# nikos 1 45386
# nikos 2 old db 45407


blast_tab=read.table(input1, sep='\t')

# blast_tab=read.table("t2t-col.20210610.fasta.FL_AGAINST_MG_BLAST_OUTPUT.txt", sep='\t')

fl_length=read.table(input2, sep='\t')

# fl_length=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH_RENAMED_LENGTH.txt", sep='\t')

new_blast_tab=merge(blast_tab, fl_length, by.x = c("V1"), by.y =  c("V1"), all.x=TRUE)

colnames(new_blast_tab)=c("query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "evalue", "bitscore", "fulllength_length")

new_blast_tab=new_blast_tab %>%
  mutate(query_length=as.numeric(as.character(qend)) - as.numeric(as.character(qstart)),
         coverage=as.numeric(as.character(query_length)) / as.numeric(as.character(fulllength_length)))

# high_quality_elem=new_blast_tab[which(new_blast_tab$coverage>=0.99),]

high_quality_elem=new_blast_tab[which(new_blast_tab$coverage>=as.numeric(input3)),]

test=high_quality_elem %>% dplyr::select(sseqid, sstart, send, coverage, query) 

test=test %>% 
  mutate(i_start = case_when(sstart < send ~ sstart, send < sstart ~ send))

test=test %>% mutate(i_end = case_when(sstart < send ~ send, send < sstart ~ sstart))

test=test %>% mutate(strand = case_when(sstart < send ~ "+", send < sstart ~ "-"))

test=test %>% select(sseqid, i_start, i_end, coverage, query, strand)

colnames(test)=c("seqnames", "start", "end", "coverage", "query", "strand")

test_irange=test %>% 
  as_granges()

test_disjoin=reduce(test_irange,with.revmap=TRUE)

list_revmap=as.data.frame(mcols(test_disjoin))

filtered_data=c()
for(i in 1:nrow(list_revmap)){
  filtered_data=c(filtered_data, (slice(list_revmap, i) %>% 
                                    unlist(use.names=FALSE))[which.max(slice(test,  slice(list_revmap, i) 
                                                                             %>% unlist(use.names=FALSE))$coverage)])
}

best_hits=slice(test, filtered_data)

# -1 sthn arxh twn d alla kai twn p
best_hits$start=as.numeric(as.character(best_hits$start)) - 1

d_strand=best_hits[which(best_hits$strand=="+"),]
p_strand=best_hits[which(best_hits$strand=="-"),]

d_strand=d_strand %>%
  # mutate(lineage="D_Athila", species_code=as.character(input5)) %>%
  mutate(lineage = paste0("D_", as.character(input4)), species_code=as.character(input5))


  #### Edw input 4 lineage

# d_strand=d_strand %>%
#   mutate(atha="D_Athila", species_code="Atha")

d_strand=unite(d_strand, coords, c(start, end), remove=FALSE, sep="-")
d_strand=unite(d_strand, chr_coords, c(seqnames, coords), remove=FALSE, sep=".")
d_strand=unite(d_strand, full_id, c(chr_coords, lineage), remove=FALSE, sep="_")
d_strand=unite(d_strand, sp_full_id, c(species_code, full_id), remove=FALSE, sep="_")

d_strand=d_strand %>%
  select(seqnames, start, end, sp_full_id, coverage, query, strand)

dblu=setdiff(d_strand$sp_full_id, fl_length$V1)
dblu=as.data.frame(dblu)

d_strand_new=subset(d_strand, sp_full_id %in% dblu$dblu)

# write.table(d_strand_new, "t2t-col.20210610.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(d_strand_new, input6, sep = "\t", row.names = F, quote= F, col.names = F)



p_strand=p_strand %>%
  # mutate(lineage="P_Athila", species_code=as.character(input5)) %>%
  mutate(lineage = paste0("P_", as.character(input4)), species_code=as.character(input5))

# p_strand=p_strand %>%
#    mutate(atha="P_Athila", species_code="Atha")

p_strand=unite(p_strand, coords, c(start, end), remove=FALSE, sep="-")
p_strand=unite(p_strand, chr_coords, c(seqnames, coords), remove=FALSE, sep=".")
p_strand=unite(p_strand, full_id, c(chr_coords, lineage), remove=FALSE, sep="_")
p_strand=unite(p_strand, sp_full_id, c(species_code, full_id), remove=FALSE, sep="_")

p_strand=p_strand %>%
  select(seqnames, start, end, sp_full_id, coverage, query, strand)

pblu=setdiff(p_strand$sp_full_id, fl_length$V1)
pblu=as.data.frame(pblu)

p_strand_new=subset(p_strand, sp_full_id %in% pblu$pblu)

# write.table(p_strand_new, "t2t-col.20210610.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(p_strand_new, input7, sep = "\t", row.names = F, quote= F, col.names = F)



#########  BED FOR TSDs

tsd_bed_non=rbind(d_strand, p_strand)

tsd_bed_non=tsd_bed_non %>%
  mutate(TSD_start=as.numeric(as.character(start)) - 5, TSD_end=as.numeric(as.character(end)) + 5)

tsd_bed_V2=tsd_bed_non %>%
  select(seqnames, TSD_start, start, sp_full_id, coverage, strand)

tsd_bed_V2_D=tsd_bed_V2[which(tsd_bed_V2$strand=="+"),]
tsd_bed_V2_P=tsd_bed_V2[which(tsd_bed_V2$strand=="-"),]


# write.table(tsd_bed_V2_D, "t2t-col.20210610.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V2_D, input8, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(tsd_bed_V2_P, "t2t-col.20210610.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_LEFTV2.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V2_P, input9, sep = "\t", row.names = F, quote= F, col.names = F)



tsd_bed_V3=tsd_bed_non %>%
  select(seqnames, end, TSD_end, sp_full_id, coverage, strand)


tsd_bed_V3_D=tsd_bed_V3[which(tsd_bed_V3$strand=="+"),]
tsd_bed_V3_P=tsd_bed_V3[which(tsd_bed_V3$strand=="-"),]


# write.table(tsd_bed_V3_D, "t2t-col.20210610.fasta.D_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V3_D, input10, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(tsd_bed_V3_P, "t2t-col.20210610.fasta.P_FULLLENGTH_BLAST_NONOVERLAPPING_PLUS_TSD_RIGHTV3.bed", sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V3_P, input11, sep = "\t", row.names = F, quote= F, col.names = F)



######### OVERLAPPING ELEMENTS


overlapping=high_quality_elem %>%
  select(sseqid, sstart, send, coverage)

overlapping=overlapping %>%
  mutate(distance_ss=as.numeric(as.character(send)) - as.numeric(as.character(sstart)))

overlapping_D=overlapping[which(overlapping$distance_ss>0),]
overlapping_P=overlapping[which(overlapping$distance_ss<0),]

overlapping_D=overlapping_D %>%
  # mutate(strand="+", lineage="D_Athila") %>%
  mutate(strand="+", lineage = paste0("D_", as.character(input4)))

overlapping_P=overlapping_P %>%
  # mutate(strand="-", lineage="P_Athila") %>%
  mutate(strand="-", lineage = paste0("P_", as.character(input4)))

overlapping_D=unite(overlapping_D, coords, c(sstart, send), remove=FALSE, sep="-")
overlapping_D=unite(overlapping_D, chr_coords, c(sseqid, coords), remove=FALSE, sep=".")
overlapping_D=unite(overlapping_D, full_id, c(chr_coords, lineage), remove=FALSE, sep="_")

overlapping_D=overlapping_D %>%
  select(sseqid, sstart, send, full_id, coverage, strand)

overlapping_P=unite(overlapping_P, coords, c(send, sstart), remove=FALSE, sep="-")
overlapping_P=unite(overlapping_P, chr_coords, c(sseqid, coords), remove=FALSE, sep=".")
overlapping_P=unite(overlapping_P, full_id, c(chr_coords, lineage), remove=FALSE, sep="_")

overlapping_P=overlapping_P %>%
  select(sseqid, send, sstart, full_id, coverage, strand)

colnames(overlapping_P)=c("sseqid", "sstart", "send", "full_id", "coverage", "strand")


full_overlapping=rbind(overlapping_D, overlapping_P)


write.table(full_overlapping, input12, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(full_overlapping, "t2t-col.20210610.fasta.DP_FULLLENGTH_BLAST_OVERLAPPING.bed", sep = "\t", row.names = F, quote= F, col.names = F)



####### check overlaps based on intersect

# try=read.table("NONOVER_AND_OVER.bed", sep="\t")
# head(try)
# 
# try=try %>%
#   mutate(dist_st=as.numeric(as.character(V2)) - as.numeric(as.character(V8)),
#          dist_en=as.numeric(as.character(V3)) - as.numeric(as.character(V9)),
#          skata=as.numeric(as.character(V3)) - as.numeric(as.character(V8)))
# 
# try[,13:15]
# 
# int_try=try[which(try$dist_st!=0 & try$dist_en!=0),]
# int_try


### 8 overlaps --> 2/66 and 3/199 have a partial overlap (size of a LTR) --> shared LTR

