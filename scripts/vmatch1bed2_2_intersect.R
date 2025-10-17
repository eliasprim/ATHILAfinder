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
input4 = args[4]
input5 = args[5]
input6 = args[6]
input7 = args[7]
input8 = args[8]
input9 = args[9]

# STEP ARGUMENT FOR BASH SCRIPT STEP 2

############ FILTER BEDTOOLS OUTPUT AND PREPARE BED FOR BEDTOOLS INTERSECT ANALYSIS

after_window1=read.table(input1, sep="\t")

# after_window1=read.table("t2t-col.20210610.fasta.F2B.PBS_2_PPT_FULL_D_r10k.txt", sep="\t")

after_window1=after_window1 %>%
  mutate(distance=as.numeric(as.character(V12)) - as.numeric(as.character(V3)))

after_window1=after_window1 %>%
  add_count(V4, name="first_signature_count") %>%
  add_count(V14, name="second_signature_count")

after_window2=read.table(input2, sep="\t")

# after_window2=read.table("t2t-col.20210610.fasta.F2B.PBS_2_PPT_FULL_D_r10k_v.txt", sep="\t")

after_window2=after_window2 %>%
  mutate(V11="NA", V12="NA", V13="NA", V14="NA", V15="NA", V16="NA", V17="NA", V18="NA", V19="NA", V20="NA", distance="NA", first_signature_count=1, second_signature_count=0)

after_window3=read.table(input3, sep="\t")

# after_window3=read.table("t2t-col.20210610.fasta.F2B.PPT_2_PBS_FULL_D_l10k_v.txt", sep="\t")

colnames(after_window3)=c("V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20")

after_window3=after_window3 %>%
  mutate(V1="NA", V2="NA", V3="NA", V4="NA", V5="NA", V6="NA", V7="NA", V8="NA", V9="NA", V10="NA", distance="NA", first_signature_count=0, second_signature_count=1)

after_window3=after_window3 %>%
  select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, distance, first_signature_count, second_signature_count)

ultimate_table=rbind(after_window1, after_window2, after_window3)

ultimate_table$first_signature_count[ultimate_table$first_signature_count>=2]="two_or_more_times"
ultimate_table$first_signature_count[ultimate_table$first_signature_count==0]="no"
ultimate_table$first_signature_count[ultimate_table$first_signature_count==1]="one_time"

ultimate_table$second_signature_count[ultimate_table$second_signature_count>=2]="two_or_more_times"
ultimate_table$second_signature_count[ultimate_table$second_signature_count==0]="no"
ultimate_table$second_signature_count[ultimate_table$second_signature_count==1]="one_time"


# write everything summary table

write.table(ultimate_table, input4, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(ultimate_table, "t2t-col.20210610.fasta.F2B.PBS_2_PPT_D_EVERYTHING.txt", sep = "\t", row.names = F, quote= F, col.names = F)



rejected_everything=ultimate_table[which(ultimate_table$first_signature_count!="one_time" | ultimate_table$second_signature_count!="one_time"),]

rejected_everything=rejected_everything %>%
  select(V1, V2, V3, V4, V5, V6, V8, V7, V14, V15, V16, V18, V17, distance, first_signature_count, second_signature_count, V10)

colnames(rejected_everything)[3]="V13"


single_ultimate_table=ultimate_table[which(ultimate_table$first_signature_count=="one_time" & ultimate_table$second_signature_count=="one_time"),]

single_ultimate_table=single_ultimate_table %>%
  select(V1, V2, V13, V4, V5, V6, V8, V7, V14, V15, V16, V18, V17, distance, first_signature_count, second_signature_count, V10)

single_ultimate_table$distance=as.numeric(as.character(single_ultimate_table$distance))


distance=single_ultimate_table %>%
  select(V4, V7, V14, V17, distance, V10)

write.table(distance, input5, sep = "\t", row.names = F, quote= F, col.names = F)


single_ultimate_table1=single_ultimate_table[which(single_ultimate_table$distance>=as.numeric(input6) & single_ultimate_table$distance<=as.numeric(input7)),]

write.table(single_ultimate_table1, input8, sep = "\t", row.names = F, quote= F, col.names = F)


rejected_elements=single_ultimate_table[which(single_ultimate_table$distance<=as.numeric(input6) | single_ultimate_table$distance>=as.numeric(input7)),]

all_rejected=rbind(rejected_elements, rejected_everything)

write.table(all_rejected, input9, sep = "\t", row.names = F, quote= F, col.names = F)

# write.table(single_ultimate_table, "t2t-col.20210610.fasta.F2B.PBS_2_PPT_D_INTERSECT.bed", sep = "\t", row.names = F, quote= F, col.names = F)



### try_kati

# pinakas=read.table("t2t-col.20210610.fasta.F2B.ALL_PBSPPT_DP_WINDOW2_SORTED.bed", sep='\t')
# 
# head(pinakas)
# 
# try_kati=pinakas %>%
#   select(V1, V2, V3, V4, V9, V14, V15)
# 
# write.table(try_kati, "t2t-col.20210610.fasta.TRY_KATI.bed", sep = "\t", row.names = F, quote= F, col.names = F)





# pin=read.table("t2t-col.20210610.fasta.F2B.PBS_PPT_DP_DISTANCE.txt", sep = "\t")
# 
# plus_pin=pin[which(pin$V6=="+"),]
# 
# plus_pin=plus_pin %>%
#   select(V2, V5, V6)
# 
# minus_pin=pin[which(pin$V6=="-"),]
# 
# minus_pin=minus_pin %>%
#   select(V4, V5, V6)
# 
# colnames(minus_pin)=c("V2", "V5", "V6")
# 
# asd=rbind(plus_pin, minus_pin)
# 
# asd$V2=as.factor(asd$V2)
# xlabs=paste(levels(asd$V2), "_",table(asd$V2), sep="")
# 
# ggplot(asd, aes(V2, V5, fill=V2)) +
#   geom_violin(width=1) +
#   geom_boxplot(width=0.2, color="orange", alpha=0.6) +
#   scale_fill_viridis(discrete = TRUE) +
#   theme(legend.position="none") +
#   ggtitle("") +
#   ylab("Internal Length") +
#   xlab("PBS Type_Counts") +
#   theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold"), plot.title = element_text(size=1)) +
#   scale_x_discrete(labels=xlabs)


