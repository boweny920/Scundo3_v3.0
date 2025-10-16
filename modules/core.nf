process samplesheet_process {

    label 'small_mem'

    input:
    val fcid

    output:
    path "samplesheet.tsv"

    script:

    if (!params.samplereport) {

        if (!params.molng) {
        """
        samplesheet_Make_V5.py --fcid ${fcid} -i ${params.roboSamplesheet} -s ${params.species} -g ${params.genomeVer} -a ${params.annotation}
        """
        } else { // This is to account for the situation when 1 FCID can have 2 or more MOLNG-ID
        """
        samplesheet_Make_V5.py --fcid ${fcid} -i ${params.roboSamplesheet} --molng ${params.molng} -s ${params.species} -g ${params.genomeVer} -a ${params.annotation}
        """
        }
    } else {
        """
        samplesheet_Make_withSampleReport.py -s ${params.samplereport} -i ${params.roboSamplesheet} --requester ${params.requester} --lab ${params.lab}
        """
    }

}

process STAR_n_mkdirs {

    memory { 50.GB * task.attempt }
    cpus { 10 * task.attempt }
    // time { 1.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 4

    // conda "${projectDir}/assets/conda_star_2_7_3a.yml" // This is no longer necessary, as new indexes are not built with this older STAR version

    input:
    val meta

    output: 
    tuple val(meta), val(output_lib_folder), val("${meta.ID}.Aligned.sortedByCoord.out.bam"), val("${meta.ID}.Aligned.toTranscriptome.out.bam"), path("star_version.log"), path("${meta.ID}.*"), val(index_genome), val(genomeVer)

    publishDir "${output_lib_folder}", mode: 'copy'

    script:
    output_lib_folder = file("${params.outdir}/${meta.pi_name}/${meta.requester_name}/${meta.molngID}.${meta.genome_ver}.${meta.annotation}/${meta.ID}" )
    index_genome = "${meta.species}/${meta.genome_ver}"
    genomeVer = "${meta.genome_ver}"
    output_lib_folder.mkdirs()

    """
    ml STAR 

    STAR --version > star_version.log

    STAR --readFilesIn ${meta.fastqs} \
		--genomeDir ${params.indexDir}/${index_genome}/annotation/${meta.annotation}/STAR_76bp \
		--runThreadN 10 \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ${meta.ID}. \
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
        --limitBAMsortRAM 20000000000 \
        --outSAMattributes NH HI MD AS nM \
        --quantMode TranscriptomeSAM GeneCounts
    """
}

process RSEM_calculate {

    label 'big_mem'

    input:
    tuple val(meta), val(output_lib_folder), val(star_sortedByCoord_out_bam), val(star_bamtoTranscriptome_out_bam), path(star_version_logFile), path(star_all_outFiles), val(index_genome), val(genomeVer)

    publishDir "${output_lib_folder}", mode: "copy"

    output:
    tuple path("${meta.ID}.*"), path("rsem_version.log")

    script:

    if (meta.readType.contains("PairEnd"))
    """
    ml rsem

    rsem-calculate-expression --version > rsem_version.log

    rsem-calculate-expression --no-bam-output --estimate-rspd  --strandedness reverse --paired-end \
    --bam ${star_bamtoTranscriptome_out_bam} ${params.indexDir}/${index_genome}/annotation/${meta.annotation}/RSEM/${genomeVer}.${meta.annotation}.RSEM \
    ${meta.ID}.RSEM
    """
    else
    """
    ml rsem

    rsem-calculate-expression --version > rsem_version.log

    rsem-calculate-expression --no-bam-output --estimate-rspd  --strandedness reverse \
    --bam ${star_bamtoTranscriptome_out_bam} ${params.indexDir}/${index_genome}/annotation/${meta.annotation}/RSEM/${genomeVer}.${meta.annotation}.RSEM \
    ${meta.ID}.RSEM
    """
}

process ERCC_control {

    label 'big_mem'

    input:
    tuple val(meta), val(output_lib_folder), val(star_sortedByCoord_out_bam), val(star_bamtoTranscriptome_out_bam), path(star_version_logFile), path(star_all_outFiles), val(index_genome), val(genomeVer)

    publishDir "${output_lib_folder}", mode: "copy"

    output:
    path "${meta.ID}.*"

    script:
    """
    ml bowtie2 samtools
    
    bowtie2 -x /n/data1/genomes/indexes/ERCC/bowtie2/ERCC92 -p 4 \
		-U ${meta.fastqs} \
		-S ERCC_control.sam 2>&1 | tee -a ${meta.ID}.ercc.log

    samtools view -bS ERCC_control.sam > ERCC_control.bam

    samtools sort ERCC_control.bam -o ERCC_control.sorted.ercc.bam

    samtools index ERCC_control.sorted.ercc.bam

    samtools idxstats ERCC_control.sorted.ercc.bam > ${meta.ID}.ercc.stats

    rm ERCC_control.bam
    rm ERCC_control.sam
    rm ERCC_control.sorted.ercc.bam

    bamTobw_rnaseq.r ${star_sortedByCoord_out_bam} reverse
    """
}

process picard_run {

    label 'big_mem'

    input:
    tuple val(meta), val(output_lib_folder), val(star_sortedByCoord_out_bam), val(star_bamtoTranscriptome_out_bam), path(star_version_logFile), path(star_all_outFiles), val(index_genome), val(genomeVer)

    publishDir "${output_lib_folder}", mode: "copy"

    output:
    path "${meta.ID}.*"

    script:

    """
    ml picard

    picard CollectRnaSeqMetrics input=${star_sortedByCoord_out_bam} \
    output=${meta.ID}.rnaseq.stats strand_specificity=FIRST_READ_TRANSCRIPTION_STRAND \
    ref_flat=${params.indexDir}/${index_genome}/annotation/${meta.annotation}/tables/${genomeVer}.${meta.annotation}.refFlat.txt\
    ribosomal_intervals=${params.indexDir}/${index_genome}/annotation/${meta.annotation}/extras/${genomeVer}.${meta.annotation}.riboList.default.txt 2>&1 | tee -a ${meta.ID}.picard.log
    """
}
