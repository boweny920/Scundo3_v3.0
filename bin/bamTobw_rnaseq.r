#!/usr/bin/env Rscript

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(csaw)

args <- commandArgs(trailingOnly=TRUE)
bam.files <- args[1]
strand <- args[2]


#bam.files     <- list.files(getwd(), pattern="primary_hits.bam$", full.names=T)
condition_name <- bam.files
condition_name <- gsub(".Aligned.sortedByCoord.out.bam","",condition_name)


#bam_stat <- read.delim("bam_stats.txt")
#total <- bam_stat$total_alignments

cur<-readGAlignments(bam.files)
size <- width(cur)[1]
norm.factor<-10^6
total <- length(cur)
#####RNAseq tracks#####
if(strand == "reverse"){
	#seperate strand for stranded RNA-seq tracks
	cur.pos<-cur[strand(cur) == "-"]
	cur.neg<-cur[strand(cur) == "+"]

    cvg.pos<-coverage(cur.pos)
    cvg.neg<-coverage(cur.neg)
    #total<- sum(as.numeric(sum(cvg)))
    cvglist.pos<-cvg.pos/total*norm.factor
    cvglist.neg<-cvg.neg/total*norm.factor*(-1)

    export.bw(cvglist.pos,paste0(condition_name,".pos.rpm.bw"))
	export.bw(cvglist.neg,paste0(condition_name,".neg.rpm.bw"))

}else if(strand == "no"){

	cvg<-coverage(cur)
	cvglist<-cvg/total*norm.factor
	export.bw(cvglist,paste0(condition_name,".rpm.bw"))

}
