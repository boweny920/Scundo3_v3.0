process bowtie2_samtools_n_mkdirs {
    
    label 'big_mem'

    input:
    val meta

    output:
    tuple val("${scundo_outDir}"), path("${meta.scundoname}.sorted.bam"), path("${meta.scundoname}.sorted.bam.bai"), path("bam_stats.txt"), path("${meta.scundoname}_bowtie2.stat"), path("${meta.scundoname}.chancedata"), path("${meta.scundoname}.rpm.bw"), path("${meta.scundoname}_cross_corr.data"), path("${meta.scundoname}_cross_corr.jpg")

    publishDir "${output_lib_folder}", mode: 'copy'

    script:

    output_lib_folder = file("${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${meta.genome_ver}.${meta.annotation}/${meta.ID}" )
    index_genome = "${meta.species}/${meta.genome_ver}"
    genomeVer = "${meta.genome_ver}"
    scundo_outDir ="${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${meta.genome_ver}.${meta.annotation}" 

    output_lib_folder.mkdirs()
    
    // Check if data is paired or single reads 
    if (meta.readType.contains("PairEnd")) {
        reads = meta.fastqs.split(' ')
        read1 = reads[0]
        read2 = reads[1]
        read_flag = "-1 ${read1} -2 ${read2}"
    } else {
        read_flag = "-U ${meta.fastqs}"
    }

    """
    ml bowtie2 samtools

    bowtie2 -x ${params.indexDir}/${index_genome}/bowtie2/${genomeVer} -p 4 \
        ${read_flag} -S ${meta.scundoname}.sam 2>&1 | tee -a ${meta.scundoname}_bowtie2.stat
    
    samtools view -bS ${meta.scundoname}.sam > ${meta.scundoname}.bam

    samtools sort ${meta.scundoname}.bam -o ${meta.scundoname}.sorted.bam

    samtools index ${meta.scundoname}.sorted.bam

    /opt/apps/dev/containers/clean_ngs/1.0/bin/bam_stats ${meta.scundoname}.sorted.bam

    bamTobw.r ${meta.scundoname}.sorted.bam pair ChIP-seq ${genomeVer} ${index_genome}
    """
}

process deep_tools_correlation {
    
    label 'deeptools'

    input:
    val scundo_outDir

    output:
    path "all_samples_bamsummary.npz" , emit: npz
    path "all_samples_bamsummary_heatmap.png" , emit: heatmap
    path "all_samples_bamsummary_pca.png" , emit: pca

    publishDir "${scundo_outDir}", mode: 'copy'

    script:
    """
    ml deeptools

    multiBamSummary bins --bamfiles ${scundo_outDir}/*/*.sorted.bam -o all_samples_bamsummary.npz -p 30 -e
    plotCorrelation --corData all_samples_bamsummary.npz --plotFile all_samples_bamsummary_heatmap.png --whatToPlot heatmap -c spearman
    plotPCA --corData all_samples_bamsummary.npz --plotFile all_samples_bamsummary_pca.png
    """
}

process deep_tools_coverage {

    label 'deeptools'

    input:
    val scundo_outDir

    output:
    path "all_samples_bamsummary_coverage.png", emit: png
    path "coverage.tsv" , emit: tsv

    publishDir "${scundo_outDir}", mode: 'copy'
    
    script:
    """
    ml deeptools
    
    plotCoverage -b ${scundo_outDir}/*/*.sorted.bam --plotFile all_samples_bamsummary_coverage --outRawCounts coverage.tsv
    """ 
}