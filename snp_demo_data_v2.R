##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%
##
## Demostration figures for variant calling algorithm validation review
##
## written by: Nate Olson
## NIST
## June 2014
## modified to try and incorporate uncertainties
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#########################%%%%%%%

library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)

setwd("demo_vcfs/")



#######
## Import of vcf files from GATK variant caller
## input file names need to be the following format or change the input loop below
## ORG_LAB_PLATFORM_pileup.expand.tsv
#######

snp_vcf_files <- list.files()[grep("^S.*vcf$", c(list.files()), value = F)]

#%#%#% generating a dataframe combining the pileup parse files
parse_vcf <- function(file){
    fileTable <- read.table(file, header = F, sep = "\t")
    colnames(fileTable) <- c("Control", "Location", "ID", "reference","alternative","qual", 
                             "filter","INFO","FORMAT","other")
    fileTable$dataset <- file
    return(fileTable)
}
snp_calls <- ldply(snp_vcf_files,.fun = parse_vcf)


snp_calls$dataset <- gsub("_L001_","",snp_calls$dataset)
snp_calls$dataset <- gsub("R1_001","",snp_calls$dataset)
head(snp_calls)
library(plyr)
library(stringr)
# defining ture variants as variants with median quality scores of greater than 2000
snp_summary <- ddply(snp_calls, .(Control, Location), summarize, count = length(dataset), medqual = median(qual))
snp_summary$ref_call <- 1
snp_summary$ref_call[snp_summary$medqual < 2500] <- 0


snp_calls$dataset <- str_replace(snp_calls$dataset,".vcf",replacement="")
snp_calls$dataset <- str_replace(snp_calls$dataset,".sort",replacement="")
snp_calls$dataset <- str_replace(snp_calls$dataset,".realign",replacement="")

snp_calls <- subset(snp_calls, select = c(Control, Location, qual, dataset))

#setting max value of qual to 5000
snp_calls$qual[snp_calls$qual > 5000] <- 5000
snp_calls <- dcast(snp_calls, Control*Location~dataset, value.var="qual", fill=0)
snp_calls <- join(snp_calls, snp_summary)
snp_calls <- subset(snp_calls, select = -c(count,medqual, Control, Location))

write.csv(snp_calls, "../snp_demo_data.csv", row.names=F)
