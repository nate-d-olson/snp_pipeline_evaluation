##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Demonstration figures for variant calling algorithm validation review
##
## written by: Nate Olson
## NIST
## June 2014
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%

library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)

#######
## Import of vcf files from GATK variant caller
## input file names need to be the following format or change the input loop below
## ORG_LAB_PLATFORM_pileup.expand.tsv
#######

snp_vcf_files <- list.files()[grep("^S.*vcf$", c(list.files()), value = F)]

#%#%#% generating a dataframe combining the pileup parse files

snp_caller_gatk = data.frame()
for(i in 1:length(snp_vcf_files)){
  fileTable <- read.table(snp_vcf_files[i], header = F, sep = "\t")
  colnames(fileTable) <- c("Control", "Location", "ID", "reference","alternative","qual", 
                           "filter","INFO","FORMAT","other")
  fileTable$dataset <- snp_vcf_files[i]#gsub(".sort.realign.vcf", "",snp_vcf_files[i])
  snp_caller_gatk <- rbind(snp_caller_gatk, fileTable)
}
remove(fileTable, i, snp_vcf_files)

snp_caller_gatk$dataset <- gsub("_L001_","",snp_caller_gatk$dataset)
snp_caller_gatk$dataset <- gsub("R1_001","",snp_caller_gatk$dataset)


# defining ture variants as variants with median quality scores of greater than 2000
snp_summary <- ddply(snp_caller_gatk, .(Control, Location), summarize, count = length(dataset), medqual = median(qual))
snp_summary$ref_call <- 1
snp_summary$ref_call[snp_summary$medqual < 2500] <- 0
#resulting in 48 positive and 41 negatives

snp_caller_gatk$dataset <- str_replace(snp_caller_gatk$dataset,".vcf",replacement="")
snp_caller_gatk$dataset <- str_replace(snp_caller_gatk$dataset,".sort",replacement="")
snp_caller_gatk$dataset <- str_replace(snp_caller_gatk$dataset,".realign",replacement="")
snp_calls <- subset(snp_caller_gatk, 
                    select = c(Control, Location, qual, dataset))
snp_calls$class <- "UG"
snp_calls$class[grep("hc", snp_calls$dataset)] <- "HC"

#setting max value of qual to 5000
snp_calls$qual[snp_calls$qual > 5000] <- 5000
snp_calls <- dcast(snp_calls, Control*Location~dataset, value.var="qual", fill=0)
snp_calls <- join(snp_calls, snp_summary)
snp_calls <- subset(snp_calls, select = -c(count,medqual, Control, Location))
# Need to figure out how to incorporate the random and perfect datasets into the ROCR analysis
write.csv(snp_calls, "~/Desktop/snp_demo_data.csv", row.names=F)