#!/usr/bin/env Rscript

library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(csaw)

args <- commandArgs(trailingOnly=TRUE)
bam_files <- args[1]
type <- args[2]
sample_type <- args[3]
genome <- args[4]
index_genome <- args[5]

print(type)
print(sample_type)
print(genome)

condition_name <- sub(".sorted.bam","",bam_files)


if(type=="single"){
    param <- readParam(minq=20)
}else if(type=="pair"){
    param <- readParam(max.frag=400,pe="both",minq=20)
}

max.delay <- 500
dedup.on <- reform(param,dedup=TRUE)

##cross-correlation Plot
x <- correlateReads(bam_files,max.delay,param=dedup.on)

data <- as.data.frame(x)
colnames(data) <- condition_name
write.table(data,paste0(condition_name,"_cross_corr.data"),quote=F,row.names=F,col.names=T)

jpeg(paste0(condition_name,"_cross_corr.jpg"))
    plot(0:max.delay,x,type="l",ylab="CCF",xlab=paste0(condition_name," Delay(bp)"))
dev.off()

if(sample_type=="Input" | type=="pair"){

	cur<-as(readGAlignments(bam_files), "GRanges")
	size <- width(cur)[1]
	norm.factor<-10^6
    cvg<-coverage(cur)
    cvglist<-cvg/length(cur)*norm.factor
    # export.bw(cvglist,paste0("../bigwig_tracks/",condition_name,".rpm.bw"))
    export.bw(cvglist,paste0("./",condition_name,".rpm.bw"))

}else{
	size <- maximizeCcf(x)
	norm.factor<-10^6

    cur<-as(readGAlignments(bam_files), "GRanges")
    cur.size<-resize(cur, size)
    cvg<-coverage(cur.size)
    cvglist<-cvg/length(cur)*norm.factor
    export.bw(cvglist,paste0("./",condition_name,".ex",size,".rpm.bw"))

}

#chance plot
all_files <- list.files(paste0("/n/analysis/genomes/",index_genome))
fa_fai_files <- grep("\\.fa\\.fai$", all_files, value = TRUE)
fa_fai_file_path <- paste0("/n/analysis/genomes/",index_genome, "/", fa_fai_files)
info <- read.delim(fa_fai_file_path,header=F)
x <- Seqinfo(seqnames=as.character(info[,1]), seqlengths=info[,2], genome=genome)
tiles <- unlist(tileGenome(x, ntile=10000))
binned <- countOverlaps(tiles, cur)
chance <- cumsum(as.numeric(sort(binned)))
chance <-chance/max(chance)
chance <- as.data.frame(chance)
colnames(chance) <- condition_name
write.table(chance,paste0(condition_name,".chancedata"),quote=F,row.names=F,col.names=T)




#chrname <- names(cvglist)

#fetchExtendedChromInfoFromUCSC(genome)
#calculate_bin <- function(chrname){
#	all_binned <- NULL
#	for(i in 1:length(chrname)){
#		name <- eval(parse(text=paste0("cvglist$",chrname[i])))
#		bin <- floor(length(as.vector(name))/1000)

#		binned <- NULL

#		for(j in 1:1000){

#			temp <- sum(as.vector(name)[j:bin*j])

#			binned <- c(binned,temp)
#		}

#		all_binned <- c(all_binned,binned)

#	}
#	return(all_binned)
#}

#data <- mclapply(chrname,calculate_bin,mc.cores=length(chrname))

#meta <- as.data.frame(meta)
#colnames(meta) <- condition_name
#write.table(meta,paste0(condition_name,".metadata"),quote=F,row.names=F,col.names=T)


####MetaGene Plot
#bam <- windowCounts(bam_files,width=50,spacing=50,filter=20,param=dedup.on)
#rwsum <- rowSums(assay(bam))
#maxed <- findMaxima(rowRanges(bam), range=1000, metric=rwsum)
#meta <- profileSites(bam_files, rowRanges(bam)[maxed],param=dedup.on, weight=1/rwsum[maxed])
#meta <- as.data.frame(meta)
#colnames(meta) <- condition_name
#write.table(meta,paste0(condition_name,".metadata"),quote=F,row.names=F,col.names=T)


#Chance QC Plot

#bam <- windowCounts(bam_files,width=50000,param=dedup.on)

#x=1:nrow(assay(bam))/nrow(assay(bam))
#chance=cumsum(as.numeric(sort(assay(bam))))
#chance=y/max(y)

#chance <- cbind(as.data.frame(x),as.data.frame(y))
#colnames(chance) <- c(paste0(condition_name,"_x"),paste0(condition_name,"_y"))
#write.table(chance,paste0(condition_name,"_chance.data"),quote=F,row.names=F,col.names=T,sep="\t")
#save(bam,file=paste0(condition_name,"_bam.rda"))

#find . -iname "run_secundo.sh" | xargs grep "Rscript" >command.txt
