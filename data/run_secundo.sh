#!/bin/bash
#$ -N s_3681
#$ -V -j y
#$ -M bioinfo@stowers.org -m ae -b y
#$ -l mem_free=20G,h_vmem=40G

cd /n/core/Bioinformatics/secondary/Rohner/cc2673/MOLNG-3681.astMex_2.0.Ens_106/s_5dpfsurface1

STAR --version > star_version.log

# STAR aligner
STAR --readFilesIn /n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_1_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_2_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_3_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_4_1_ATCACG.fastq.gz \
		--genomeDir /n/analysis/indexes/astMex_2.0/annotation/Ens_106/STAR_76bp \
		--runThreadN 4 \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix s_5dpfsurface1. \
		--outSAMprimaryFlag OneBestScore \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --readFilesCommand zcat \
        --limitBAMsortRAM 10000000000 \
        --outSAMattributes NH HI MD AS nM \
        --quantMode TranscriptomeSAM GeneCounts

rsem-calculate-expression --version > rsem_version.log

#RSEM
rsem-calculate-expression --no-bam-output --estimate-rspd  --strandedness reverse \
--bam s_5dpfsurface1.Aligned.toTranscriptome.out.bam /n/analysis/indexes/astMex_2.0/annotation/Ens_106/RSEM/astMex_2.0.Ens_106.RSEM s_5dpfsurface1.RSEM
# output:
# rsem_version.log
# s_5dpfsurface1.RSEM.genes.results
# s_5dpfsurface1.RSEM.isoforms.results
# s_5dpfsurface1.RSEM.stat

# ERCC

bowtie2 -x /n/data1/genomes/indexes/ERCC/bowtie2/ERCC92 -p 4 \
		-U /n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_1_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_2_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_3_1_ATCACG.fastq.gz,/n/analysis/Rohner/cc2673/MOLNG-3681/HK5NCBGXT/n_4_1_ATCACG.fastq.gz \
		-S s_5dpfsurface1.sam 2>&1 | tee -a s_5dpfsurface1_ercc.log

samtools view -bS s_5dpfsurface1.sam > s_5dpfsurface1.bam

samtools sort s_5dpfsurface1.bam -o s_5dpfsurface1.sorted.ercc.bam

samtools index s_5dpfsurface1.sorted.ercc.bam

samtools idxstats s_5dpfsurface1.sorted.ercc.bam > s_5dpfsurface1.ercc.stats

rm s_5dpfsurface1.bam
rm s_5dpfsurface1.sam
rm s_5dpfsurface1.sorted.ercc.bam

Rscript /n/ngs/tools/secundo/scripts/bamTobw_rnaseq.r s_5dpfsurface1.Aligned.sortedByCoord.out.bam reverse

# ERCC outfiles: 
#  s_5dpfsurface1.ercc.stats
#  s_5dpfsurface1_ercc.log
#  s_5dpfsurface1.pos.rpm.bw
#  s_5dpfsurface1.neg.rpm.bw

# picard
java -Xmx10g -jar /n/apps/CentOS7/bin/picard.jar CollectRnaSeqMetrics input=s_5dpfsurface1.Aligned.sortedByCoord.out.bam \
  output=s_5dpfsurface1.rnaseq.stats strand_specificity=FIRST_READ_TRANSCRIPTION_STRAND \
  ref_flat=/n/analysis/indexes/astMex_2.0/annotation/Ens_106/tables/astMex_2.0.Ens_106.refFlat.txt ribosomal_intervals=/n/analysis/indexes/astMex_2.0/annotation/Ens_106/extras/astMex_2.0.Ens_106.riboList.default.txt 2>&1 | tee -a s_5dpfsurface1_picard.log

#picard outfiles:
# picard.log
# rnaseq.stats