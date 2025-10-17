rm(list=ls())

options(max.print=10000000)
options(scipen = 999)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))

# setwd('~/Desktop/athilafinder_ubuntu_args/')

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
input18 = args[18]

###### UNIVERSAL NAMES FOR FULLLENGTH, 5' AND 3' LTRs

full_bed=read.table(input1, sep="\t")

# full_bed=read.table("t2t-col.20210610.fasta.F2B.DP_FULLLENGTH.bed", sep="\t")

# full_bed=full_bed %>%
#   mutate(species_code="Atha")

full_bed=full_bed %>%
  mutate(species_code=as.character(input2))

full_bed=unite(full_bed, half_new_id, c(V2, V3), remove=FALSE, sep="-")
full_bed=unite(full_bed, new_id, c(V1, half_new_id), remove=FALSE, sep=".")
full_bed=unite(full_bed, fasta_id, c(V1, half_new_id), remove=FALSE, sep=":")

lineage=as.character(input3)

full_bed_D=full_bed[which(full_bed$V7=="+"),]
full_bed_D=full_bed_D %>%
  # mutate(strand="D_Athila") %>%
  mutate(strand=paste0("D_", lineage))
full_bed_D=unite(full_bed_D, new_full_id, c(new_id, strand), remove=TRUE, sep="_")
full_bed_D=unite(full_bed_D, ult_new_full_id, c(species_code, new_full_id), remove=TRUE, sep="_")


full_bed_P=full_bed[which(full_bed$V7=="-"),]
full_bed_P=full_bed_P %>%
  # mutate(strand="P_Athila") %>%
  mutate(strand=paste0("P_", lineage))
full_bed_P=unite(full_bed_P, new_full_id, c(new_id, strand), remove=TRUE, sep="_")
full_bed_P=unite(full_bed_P, ult_new_full_id, c(species_code, new_full_id), remove=TRUE, sep="_")


# full_bed_P$fasta_id=paste0(full_bed_P$fasta_id, " P_element")
# full_bed_P$fasta_id=paste("rc", full_bed_P$fasta_id, sep="_")

ultimate_bed=rbind(full_bed_D, full_bed_P)
ultimate_bed=ultimate_bed %>%
  select(V1, V2, V3, V4, fasta_id, ult_new_full_id, V5, V6, V7)

to_rename_ltrs=ultimate_bed %>%
  select(V4, ult_new_full_id)

to_rename=ultimate_bed %>%
  select(fasta_id, ult_new_full_id)

new_ultimate_bed=ultimate_bed %>%
  select(V1, V2, V3, ult_new_full_id, V5, V6, V7)

write.table(to_rename, input4, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new_ultimate_bed, input5, sep = "\t", row.names = F, quote= F, col.names = F)


#########  BED FOR TSDs

tsd_bed=new_ultimate_bed

tsd_bed=tsd_bed %>%
  mutate(TSD_V2=as.numeric(as.character(V2)) - 5, TSD_V3=as.numeric(as.character(V3)) + 5)

tsd_bed_V2=tsd_bed %>%
  select(V1, TSD_V2, V2, ult_new_full_id, V5, V6, V7)

tsd_bed_V2_D=tsd_bed_V2[which(tsd_bed_V2$V7=="+"),]
tsd_bed_V2_P=tsd_bed_V2[which(tsd_bed_V2$V7=="-"),]

write.table(tsd_bed_V2_D, input6, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V2_P, input7, sep = "\t", row.names = F, quote= F, col.names = F)



tsd_bed_V3=tsd_bed %>%
  select(V1, V3, TSD_V3, ult_new_full_id, V5, V6, V7)


tsd_bed_V3_D=tsd_bed_V3[which(tsd_bed_V3$V7=="+"),]
tsd_bed_V3_P=tsd_bed_V3[which(tsd_bed_V3$V7=="-"),]



write.table(tsd_bed_V3_D, input8, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(tsd_bed_V3_P, input9, sep = "\t", row.names = F, quote= F, col.names = F)




ltr5_bed=read.table(input10, sep="\t")

rename_5ltr=merge(to_rename_ltrs, ltr5_bed, by.x = c("V4"), by.y =  c("V5"), all.x = TRUE)
rename_5ltr=unite(rename_5ltr, half_5ltr_id, c(V2, V3), remove=FALSE, sep="-")
rename_5ltr=unite(rename_5ltr, fasta_5ltr_id, c(V1, half_5ltr_id), remove=FALSE, sep=":")
rename_5ltr=rename_5ltr %>%
  mutate(prime="5prime")

rename_5ltr=unite(rename_5ltr, new_full_id_prime, c(ult_new_full_id, prime), remove=TRUE, sep="_")

to_5ltr_rename=rename_5ltr %>%
  select(fasta_5ltr_id, new_full_id_prime)

new_ltr5_bed=rename_5ltr %>%
  select(V1, V2, V3, new_full_id_prime, V4.y, V6)

write.table(to_5ltr_rename, input11, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new_ltr5_bed, input12, sep = "\t", row.names = F, quote= F, col.names = F)

ltr3_bed=read.table(input13, sep="\t")

rename_3ltr=merge(to_rename_ltrs, ltr3_bed, by.x = c("V4"), by.y =  c("V5"), all.x = TRUE)
rename_3ltr=unite(rename_3ltr, half_3ltr_id, c(V2, V3), remove=FALSE, sep="-")
rename_3ltr=unite(rename_3ltr, fasta_3ltr_id, c(V1, half_3ltr_id), remove=FALSE, sep=":")
rename_3ltr=rename_3ltr %>%
  mutate(prime="3prime")

rename_3ltr=unite(rename_3ltr, new_full_id_prime, c(ult_new_full_id, prime), remove=TRUE, sep="_")

to_3ltr_rename=rename_3ltr %>%
  select(fasta_3ltr_id, new_full_id_prime)

new_ltr3_bed=rename_3ltr %>%
  select(V1, V2, V3, new_full_id_prime, V4.y, V6)

write.table(to_3ltr_rename, input14, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new_ltr3_bed, input15, sep = "\t", row.names = F, quote= F, col.names = F)

int_bed=read.table(input16, sep="\t")

rename_int=merge(to_rename_ltrs, int_bed, by.x = c("V4"), by.y =  c("V4"), all.x = TRUE)
rename_int=unite(rename_int, half_int_id, c(V2, V3), remove=FALSE, sep="-")
rename_int=unite(rename_int, fasta_int_id, c(V1, half_int_id), remove=FALSE, sep=":")
rename_int=rename_int %>%
  mutate(prime="internal")

rename_int=unite(rename_int, new_full_id_prime, c(ult_new_full_id, prime), remove=TRUE, sep="_")

to_int_rename=rename_int %>%
  select(fasta_int_id, new_full_id_prime)

new_int_bed=rename_int %>%
  select(V1, V2, V3, new_full_id_prime, V5)

write.table(to_int_rename, input17, sep = "\t", row.names = F, quote= F, col.names = F)

write.table(new_int_bed, input18, sep = "\t", row.names = F, quote= F, col.names = F)
