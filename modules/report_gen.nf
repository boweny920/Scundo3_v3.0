process copy_index_files {

    label 'small_mem'

    input:
    path samplesheet
    val meta
    //This is just to make the workflow wait until the SECUNDO workflow is finished then to execute 
    val picard_out 
    val rsem_out
    val ercc_out
    //

    output:
    stdout // this will have a \n at the end! Beware!

    script:
    """
    copy_indexFiles_V2.py --samplesheet ${samplesheet} --outdir_of_secundo ${params.outdir} --annotation_Version ${meta.annotation} --indexDir ${params.indexDir} 
    """
}

process target_tsv_lims_json {

    label 'small_mem'

    input:
    val secundo_folder 
    path samplesheet

    output: 
    val "${secundo_folder}"

    // publishDir "${secundo_folder}", mode: 'copy' 

    script:
    if (!params.samplereport) {

        if (!params.molng){

            """
            target_tsv_Make_V2.py --lims ${params.fcid}

            lims_json_Make.py --lims ${params.fcid}

            cp lims_order.json targets.tsv ${secundo_folder}

            cp ${projectDir}/assets/dygraph.css ${secundo_folder}

            """
        
        } else {
            """
            target_tsv_Make_V2.py --lims ${params.fcid} --molng ${params.molng}

            lims_json_Make.py --lims ${params.fcid} --molng ${params.molng}

            cp lims_order.json targets.tsv ${secundo_folder}

            cp ${projectDir}/assets/dygraph.css ${secundo_folder}

            """
        }

    } else {
        """
        target_tsv_Make_PublicData.py -s ${samplesheet}

        lims_json_Make_PublicData.py -s ${samplesheet}

        cp lims_order.json targets.tsv ${secundo_folder}

        cp ${projectDir}/assets/dygraph.css ${secundo_folder}
        """
    }

}

process target_tsv_lims_json_PublicData {

    label 'small_mem'

    input:
    val secundo_folder 
    path samplesheet

    output: 
    val "${secundo_folder}"

    // publishDir "${secundo_folder}", mode: 'copy' 

    script:

    """
    target_tsv_Make_PublicData.py -s ${samplesheet}

    lims_json_Make_PublicData.py -s ${samplesheet}

    cp lims_order.json targets.tsv ${secundo_folder}

    cp ${projectDir}/assets/dygraph.css ${secundo_folder}

    """
    
}

process run_RmarkDown_script {

    label 'small_mem'

    input:
    val secundo_folder
    val meta

    output:
    val "https://webfs${secundo_folder.replaceAll("\\n", "")}/analysis_rnaseq.html"

    script:
    """
    cp ${projectDir}/assets/analysis_rnaseq.rmd ${secundo_folder}

    rmarkdown_run.py --secundo_dir ${secundo_folder.replaceAll("\\n", "")} --Reference_Genome ${genome_ver} --Annotation_Version ${meta.annotation}
    """
}

process biotools_orderAppend {
    label 'small_mem'

    input:
    val report_link

    output:
    val "order appended to biotools"

    script:
    order_path = report_link.replaceAll("\\n", "")
    // order_path = report_link.replaceAll("\\n", "")

    """
    echo "${order_path}" >> /n/core/Bioinformatics/tools/biotools/db/orders.txt

    awk '!seen[\$0]++' /n/core/Bioinformatics/tools/biotools/db/orders.txt > db.txt

    cp db.txt /n/core/Bioinformatics/tools/biotools/db/orders.txt

    db_nfScundo3.py /n/core/Bioinformatics/tools/biotools/db/orders.txt

    ssh -o StrictHostKeyChecking=no compbiotools "conda activate /n/core/Bioinformatics/tools/biotools/conda_env; mongoimport --db RNAseqAnalysis --collection OrderInfo --jsonArray --drop --file /n/core/Bioinformatics/tools/biotools/db/datafile.json"
    """
}

process run_RmarkDown_script_chip {
    label 'lil_mem'

    input:
    val secundo_folder
    val meta

    output:
    val "https://webfs${secundo_folder.replaceAll("\\n", "")}/analysis_chipseq.html"

    script:
    """
    cp ${projectDir}/assets/analysis_chipseq.rmd ${secundo_folder}

    rmarkdown_run_chip.py --secundo_dir ${secundo_folder.replaceAll("\\n", "")}

    """
} 
