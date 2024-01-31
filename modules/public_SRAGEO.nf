process samplesheet_process_publicSRAGEO { // Need to test!

    label 'lil_mem'

    input:
    val public_datacsv

    output:

    path "Secundo3_SampleSheet_PublicData.tsv"

    publishDir "${launchDir}" , mode: 'copy'

    script:

    """
    ml sratoolkit/3.0.1

    sra_download.dry_run.r -s ${public_datacsv} -l ${params.lab} -r ${params.requester} 

    cp ./*_samplesheet.csv ./*_cmds.sh /n/core/Bioinformatics/PublicData/analysis_sample_sheet/

    PublicData_samplesheet_Make.py -s ./*_samplesheet.csv -i ${params.roboSamplesheet} -l ${params.lab} -r ${params.requester}

    """
}