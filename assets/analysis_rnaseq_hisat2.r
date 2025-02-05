library(rjson)
library(parallel)
library(matrixStats)
library(grid)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(plyr)
library(dplyr)
library(pander)
library(dygraphs)
library(knitr)
library(rmarkdown)

theme_x <- theme_bw() + theme(text = element_text(size=25))
setwd("/n/core/Bioinformatics/analysis/CompBio/boweny/nf-Pipeline/Scundo3_v4/assets/tmp_data")


# Alignment stats
bam_stats_name  <- function(x) {
    tail(unlist(strsplit(dirname(x), "/")),1)
}

options(stringsAsFactors=FALSE, width=80, knitr.figure_dir = "secundo", mc.cores=4)

figure_path <- function(filename="") {
  paste0(getOption("knitr.figure_dir"), "/", filename)
}

load_align_summary <- function(align.summary.file, path.prefix="") {
    # parse the input line from HISAT2's align summary to determine the total
    # number of reads
    align.summary.raw     <- read.delim(paste0(path.prefix, align.summary.file), sep = " ", header = FALSE)

    # align.summary         <- unlist(lapply(align.summary.raw, function(a) { strsplit(a," ") }))
    data.frame(sample_name=bam_stats_name(align.summary.file),
            total=as.numeric(align.summary.raw[1,1]),
            unique=as.numeric(align.summary.raw[4,5]),
            uniquePercent=as.numeric(gsub("[()%]", "", align.summary.raw[4,6]))/100,
            multi=as.numeric(align.summary.raw[5,5]),
            multiPercent=as.numeric(gsub("[()%]", "", align.summary.raw[5,6]))/100
            )
}

align.stats.files <- list.files(getwd(), pattern="*.sam.log", recursive=TRUE)

align.stats.list  <- lapply(align.stats.files, load_align_summary)

align.stats       <- do.call(rbind, align.stats.list)
align.stats$unaligned <- align.stats$total - align.stats$unique - align.stats$multi
align.stats$unalignedPercent <- align.stats$unaligned/align.stats$total
align.stats$aligned <- align.stats$unique + align.stats$multi
align.stats$alignedPercent <- align.stats$aligned/align.stats$total
align.stats$totalPercent <- 1

write.table(align.stats, file=figure_path("align_stats.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

bam.stats.d <- align.stats[,c(1,2,9,3,5,7)]
bam.stats.p <- align.stats[,c(1,11,10,4,6,8)]

colnames(bam.stats.d) <- c("sample_name","total","all_aligned","unique_mapped_reads","multi_mapped_reads","unaligned")
colnames(bam.stats.p) <- c("sample_name","total","all_aligned","unique_mapped_reads","multi_mapped_reads","unaligned")

bam.stats.ar          <- melt(bam.stats.d, id.vars="sample_name")
bam.stats.ar$variable <- cap_first(bam.stats.ar$variable)
bam.stats.ar$variable <- factor(bam.stats.ar$variable, levels=alignment.metric.orders.d)

g <- ggplot(bam.stats.ar, aes(x=sample_name, y=value, fill=variable)) +
     geom_bar(stat="identity",position="dodge") +
     scale_fill_brewer("Metric", palette="Set2") +
     scale_y_continuous(name="Reads", labels = comma) +
     theme_x +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           legend.position="top") +
     guides(fill=guide_legend(title.position="top")) +
     xlab("Sample Name") 
    #  theme(plot.margin = margin(1, 1, 1, 5, "cm"))

g

folder_path <- "./"

file_path <- file.path(folder_path, "alignment_reads.pdf")

# Save the plot with specified size
ggsave(filename = file_path, plot = g, device = "pdf", width = 15, height = 10, units = "in", limitsize = FALSE)




picard.ext <- "\\.rnaseq\\.stats"

read_picard_stats <- function(picard.stats.file) {
    read.delim(picard.stats.file, skip=6, nrows=1)
}

format_name <- function(x, picard_ext=picard.ext) {
    gsub(picard_ext, "", basename(x))
}


picard.stats.files     <- list.files("/n/core/Bioinformatics/analysis/CompBio/boweny/secundo/HISAT2_Bacterial_Test", pattern=picard.ext, full.names=TRUE, recursive=TRUE)
picard.stats.dl        <- lapply(picard.stats.files, read_picard_stats)
names(picard.stats.dl) <- sapply(picard.stats.files, format_name, USE.NAMES=FALSE)
picard.stats.dx        <- do.call(rbind, picard.stats.dl)
names(picard.stats.dx) <- tolower(names(picard.stats.dx))
exclude.cols           <- c("sample","group")
include.cols           <- names(picard.stats.dx)[!(names(picard.stats.dx) %in% exclude.cols)]
picard.stats.dx        <- picard.stats.dx[, include.cols]
picard.stats.dx$sample <- as.character(rownames(picard.stats.dx))
picard.stats.dx        <- picard.stats.dx[order(picard.stats.dx$sample), ]

# plot percentages
perc.cols.i            <- c("pct_coding_bases","pct_utr_bases","pct_intronic_bases","pct_intergenic_bases","pct_ribosomal_bases")
p.stats.perc.x         <- picard.stats.dx[,c("sample",perc.cols.i)]

# add whole genome data
genome.files <- tail(system(paste("find", "/n/core/Bioinformatics/analysis/CompBio/boweny/secundo/HISAT2_Bacterial_Test", "-iname '*.GenomeBpTypes.txt' | sort"), intern=TRUE))

reformat_genome <- function(genome) {
	condition_name <- sub("/n/core/Bioinformatics/analysis/CompBio/boweny/secundo/HISAT2_Bacterial_Test","",genome)
	condition_name <- sub(".GenomeBpTypes.txt","",condition_name)
	condition_name <- sub("/","",condition_name)
	genome <- read.delim(genome)
    genome$pct_cds <- genome$PC.CDS/genome$TOTAL
    #genome$pct_utr <- (genome$PC.UTR+genome$NC.EXON)/genome$TOTAL
    genome$pct_utr <- (genome$PC.UTR)/genome$TOTAL
    genome$pct_intron <- genome$INTRON/genome$TOTAL
    #genome$pct_inter <- genome$INTERGENE/genome$TOTAL 
    genome$pct_inter <- (genome$INTERGENE+genome$NC.EXON)/genome$TOTAL
    genome$pct_rib <- genome$RIBOSOME/genome$TOTAL
    genome[1,1] <- condition_name
    rownames(genome) <- genome[,1]
    genome <- genome[,c(1,9:13)]
    colnames(genome) <- colnames(p.stats.perc.x)

    genome[1,]
}

genome.data <- mclapply(genome.files,reformat_genome)
genome.data <- do.call("rbind",genome.data)



p.stats.perc.x <- rbind(p.stats.perc.x, genome.data)

p.stats.perc.x$group  <- length(p.stats.perc.x$sample)
p.stats.perc.xm        <- melt(p.stats.perc.x, id.vars=c("sample","group"))

# plot stranded-ness
strand.cols.i          <- c("correct_strand_reads","incorrect_strand_reads")
p.stats.strand.x       <- picard.stats.dx[,c("sample",strand.cols.i)]
p.stats.strand.x$reads <- rowSums(p.stats.strand.x[,-1])

p.stats.strand.x$perc_correct_strand_reads   <- NA
p.stats.strand.x$perc_incorrect_strand_reads <- NA

p.stats.strand.x <- transform(p.stats.strand.x,
                              perc_correct_strand_reads = (correct_strand_reads / reads),
                              perc_incorrect_strand_reads = (incorrect_strand_reads / reads))

x.strand.cols <- c("sample", names(p.stats.strand.x)[grep("^perc", names(p.stats.strand.x))])

p.stats.strand.xm <- melt(p.stats.strand.x[x.strand.cols], id.vars=c("sample"))

perc.xm.lbls <- c("pct_coding_bases"="CDS",
				  "pct_utr_bases"="UTR",
                  "pct_intronic_bases"="Intronic",
                  "pct_intergenic_bases"="Intergenic",
                  "pct_ribosomal_bases"="Ribosome")

p.stats.perc.xm$variable <- factor(p.stats.perc.xm$variable, levels=names(perc.xm.lbls))

p.stats.perc.xm          <- p.stats.perc.xm[with(p.stats.perc.xm, order(sample, variable)), ]

# manually split by group and plot with flipped coordinates.
# facet_wrap doesn't work with coord_flip

plot_feat_dist <- function(group, p.stats.all) {
    p.stats.group        <- p.stats.all[p.stats.all$group == group,]
    p.stats.group$sample <- factor(p.stats.group$sample, levels=rev(unique(p.stats.group$sample)))
    ps <- ggplot(p.stats.group, aes(x=sample, y=value, fill=variable)) +
        geom_bar(stat="identity") +
        scale_y_continuous("Percentage", labels=percent) +
        labs(x="Samples") +
        scale_fill_brewer(name="Base distribution", labels=perc.xm.lbls, palette="Set2") +
        theme_x +
        coord_flip() +
        theme(strip.background = element_rect(colour='NA', fill='NA'),
                strip.text.y = element_text(colour='white'),
                legend.position="top") +
        guides(fill=guide_legend(title.position="top"))
    ps

    file_path <- file.path("./", paste0("picard_feature_distribution_", group, ".png"))
    ggsave(filename = file_path, plot = ps, width = 10, height = 10, units = "in")
}

unique_groups <- unique(p.stats.perc.xm$group)
psx           <- lapply(unique_groups, plot_feat_dist, p.stats.perc.xm)


#Strand Specificity

strand.xm.lbls <- c("perc_correct_strand_reads"="Sense Strand Reads",
                    "perc_incorrect_strand_reads"="Antisense Strand Reads")

p.st <- ggplot(p.stats.strand.xm, aes(x=sample, y=value, fill=variable)) +
        geom_bar(stat="identity") +
        scale_y_continuous(labels=percent) +
        scale_fill_brewer(name="Strand correctness", labels=strand.xm.lbls, palette="Set2") +
        labs(x="Sample", y="Percentage") +
        theme_x +
        theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=30),legend.position="top")

p.st

file_path <- file.path("./", "Strand_Specificity.png")

ggsave(filename = file_path, plot = p.st, width = 15, height = 15, units = "in")

# Picard coverage plot -- 5' to 3' tx coverage

require(dygraphs)
# read the bins...
picard.bins.dl <- lapply(picard.stats.files, function(x) {
  sample.name  <- format_name(x)
  dx           <- read.delim(x, skip=11, header=FALSE, col.names=c("bin","cov"))
  colnames(dx)    <- c("bin",sample.name)
  dx[sample.name]
})

picard.bins.dx <- do.call(cbind, picard.bins.dl)

picard.bins.dx <- data.frame(bins=seq(0,100,1),picard.bins.dx)

dygraph(picard.bins.dx, main="RNA-Seq coverage vs. Transcript Position",
		xlab="Normalized Distance Along Transcript", ylab="Normalized Coverage") %>%
		dyHighlight(highlightSeriesOpts=list(strokeWidth=3)) %>%
		dyLegend(show="always",width=950) %>%
		dyRangeSelector() %>%
		dyCSS("dygraph.css")

g <- dygraph(picard.bins.dx, main="RNA-Seq coverage vs. Transcript Position",
		xlab="Normalized Distance Along Transcript", ylab="Normalized Coverage") %>%
		dyHighlight(highlightSeriesOpts=list(strokeWidth=3)) %>%
		dyLegend(show="always",width=950) %>%
		dyRangeSelector() %>%
		dyCSS("dygraph.css")

file_path <- file.path("./", "picard_coverage_plot.html")
htmlwidgets::saveWidget(g, file_path)


#TPM FPKM cal

get_featCount_files <- function() {
  system(paste("find", getwd(), "-iname '*featCount' | sort"), intern=TRUE)
}

reformat_featCount_TPM <- function(df, label) {
    names(df)[ncol(df)] <- "counts"
    df$gene_length_kb <- df$Length / 1000
    df$rpk <- df$counts / df$gene_length_kb
    sum_rpk <- sum(df$rpk)
    df$TPM <- (df$rpk / sum_rpk) * 1e6

    df <- df[, c("Geneid", "TPM")]
    df$sample <- label
    df
}

reformat_featCount_FPKM <- function(df, label) {
    names(df)[ncol(df)] <- "counts"
    total_aligned_reads <- sum(df$counts)

    df <- df %>%
        mutate(FPKM = (counts/(Length/1000))/((total_aligned_reads / 1e6)))

    df <- df[, c("Geneid", "FPKM")]
    df$sample <- label
    df
}

load_featCount_TPM <- function(rsem.dir=getwd()) {
  rsem.files <- get_featCount_files()
  rsem.dfx.list   <- mclapply(rsem.files, function(rsem.file, rsem.dir) {
    rsem.data  <- read.delim(rsem.file, stringsAsFactors=F, header=T, comment.char = "#")
    rsem.label <- sub("^/", "",dirname(gsub(rsem.dir, "", rsem.file)))
    rsem.label <- tail(unlist(strsplit(rsem.label, "/")),1)
    reformat_featCount_TPM(rsem.data, rsem.label)
  }, rsem.dir)
  rsem.dfx.list
}

load_featCount_FPKM <- function(rsem.dir=getwd()) {
  rsem.files <- get_featCount_files()
  rsem.dfx.list   <- mclapply(rsem.files, function(rsem.file, rsem.dir) {
    rsem.data  <- read.delim(rsem.file, stringsAsFactors=F, header=T, comment.char = "#")
    rsem.label <- sub("^/", "",dirname(gsub(rsem.dir, "", rsem.file)))
    rsem.label <- tail(unlist(strsplit(rsem.label, "/")),1)
    reformat_featCount_FPKM(rsem.data, rsem.label)
  }, rsem.dir)
  rsem.dfx.list
}

#RSEM FPKM values
rsem.fpkms.list            <- load_featCount_FPKM()
rsem.fpkms.tidy            <- do.call(rbind, rsem.fpkms.list)
rsem.fpkms.tidy.c          <- rsem.fpkms.tidy %>% group_by(Geneid, sample) %>% summarise(FPKM=sum(FPKM))
rsem.fpkms.wide            <- dcast(rsem.fpkms.tidy.c, Geneid ~ sample, value.var="FPKM")

#RSEM TPM values
rsem.tpms.list            <- load_featCount_TPM()
rsem.tpms.tidy            <- do.call(rbind, rsem.tpms.list)
rsem.tpms.tidy.c          <- rsem.tpms.tidy %>% group_by(Geneid, sample) %>% summarise(TPM=sum(TPM))
rsem.tpms.wide            <- dcast(rsem.tpms.tidy.c, Geneid ~ sample, value.var="TPM")

rsem.tpms.wide.orig <- rsem.tpms.wide

# convert to log2+1 for correlation
rsem.tpms.wide[,-1]        <- log2(rsem.tpms.wide[,-1]+ 0.1)

rsem.tpms.wide[,-1][is.infinite(as.matrix(rsem.tpms.wide[,-1]))] <- 0
rsem.tpms.wide[,-1][is.nan(as.matrix(rsem.tpms.wide[,-1]))]      <- 0
rsem.tpms.wide[,-1][is.na(as.matrix(rsem.tpms.wide[,-1]))]       <- 0

rsem.tpms.log <- rsem.tpms.wide[,-1][which(rowSums(rsem.tpms.wide[,-1])>0),]

calculate_density <- function(x){
	data = density(x,n=200,from=-2,to=10.5)
	return(data$y)
}

density <- sapply(rsem.tpms.log,calculate_density)
density_data <- data.frame(bins=seq(-2,10.5,length.out=200),density)

dygraph(density_data, main="RNA-Seq log2(TPM) Density Plot",
		xlab="Log2(TPM)", ylab="Density") %>%
		dyHighlight(highlightSeriesOpts=list(strokeWidth=3)) %>%
		dyLegend(show="always",width=950) %>%
		dyRangeSelector() %>%
		dyCSS("dygraph.css")

g <- dygraph(density_data, main="RNA-Seq log2(TPM) Density Plot",
		xlab="Log2(TPM)", ylab="Density") %>%
		dyHighlight(highlightSeriesOpts=list(strokeWidth=3)) %>%
		dyLegend(show="always",width=950) %>%
		dyRangeSelector() %>%
		dyCSS("dygraph.css")

file_path <- file.path("./", "tpm_density_plot.html")
htmlwidgets::saveWidget(g, file_path)


# TPM Boxplot

row_set <- ceiling((ncol(rsem.tpms.log))/4)

data=stack(rsem.tpms.log)

dp <- ggplot(data,aes(x=ind,y=values,color=ind))+geom_boxplot()+
		labs(x="Samples", y="Log2(TPM)") +
        ggtitle("Boxplot of Log2(TPM)") +
        theme_x +
        theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
        guides(col=guide_legend(nrow=row_set,title="samples"))

suppressWarnings(print(dp))

file_path <- file.path("./", "tpm_boxplot.png")

ggsave(filename = file_path, plot = dp, width = 15, height = 15, units = "in")

# Spearman Correlations Coefficient Heatmap

cor.method <- "spearman"

rsem.tpms.exp              <- rsem.tpms.wide.orig

# convert to log2
rsem.tpms.exp[,-1]         <- log2(rsem.tpms.exp[,-1])
rsem.tpms.exp.m            <- as.matrix(rsem.tpms.exp[,-1])
rownames(rsem.tpms.exp.m)  <- rsem.tpms.exp.m[,1]

rsem.tpms.exp.m[is.infinite(rsem.tpms.exp.m)] <- NA
rsem.tpms.exp.m[is.nan(rsem.tpms.exp.m)]      <- NA

rsem.tpms.max.exp     <- rowMaxs(rsem.tpms.exp.m, na.rm=TRUE) > 1
rsem.tpms.exp.m       <- rsem.tpms.exp.m[rsem.tpms.max.exp,  ]

cor.exp.matrix         <- cor(rsem.tpms.exp.m, method=cor.method, use="complete.obs")

plot.title <- paste(cor.method,"correlation for ", nrow(rsem.tpms.exp.m), "expressed genes without clustering")

pheatmap(cor.exp.matrix, main=plot.title,border_color=NA, cluster_cols=F, cluster_rows=F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

g <- pheatmap(cor.exp.matrix, main=plot.title,border_color=NA, cluster_cols=F, cluster_rows=F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

file_path <- file.path("./", "Spearman_Correlations_Coefficient_Heatmap.png")

ggsave(filename = file_path, plot = g, width = 10, height = 10, units = "in")

plot.title <- paste(cor.method,"correlation for ", nrow(rsem.tpms.exp.m), "expressed genes with hierarchical clustering")

pheatmap(cor.exp.matrix, main=plot.title, border_color=NA,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

g <- pheatmap(cor.exp.matrix, main=plot.title, border_color=NA,color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

file_path <- file.path("./", "Spearman_Correlations_Coefficient_Heatmap_hier.png")

ggsave(filename = file_path, plot = g, width = 10, height = 10, units = "in")


#Count table

reformat_star <- function(df, label) {
    names(df)[ncol(df)] <- "counts"
    df <- df[,c("Geneid", "counts")]
    names(df) <- c("gene_id","counts")
    df$sample <- label
    df
}

load_star <- function(star.dir=getwd()) {
    star.files <- system(paste0("find ", getwd(), " -iname '*.featCount' | sort"), intern=TRUE)
    star.dfx.list  <- mclapply(star.files, function(star.file, star.dir) {
        star.data  <- read.delim(star.file, stringsAsFactors=F, header=T, comment.char = "#")
        star.label <- sub("^/", "",dirname(gsub(star.dir, "", star.file)))
        star.label <- tail(unlist(strsplit(star.label, "/")),1)
        reformat_star(star.data, star.label)
    }, star.dir)
    star.dfx.list
}

format_star <- function(){
	star.list <- load_star()
	star.tidy            <- do.call(rbind, star.list)
	star.tidy.c          <- star.tidy %>% group_by(gene_id, sample)
	# select <- ifelse(strand, "countNeg", "countAll") # What is this count positive and negative, this does not match stars output even!!
	star           <- dcast(star.tidy.c, gene_id ~ sample, value.var="counts")
	star
}

figure_path <- function(filename="") {
    # paste0(getOption("knitr.figure_dir"), "/", filename)
    paste0("./", filename)
}

star            <- format_star()
write.table(star, file=figure_path("star_count.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat('[STAR-Count Table](',figure_path("star_count.csv"),')\n\n')
