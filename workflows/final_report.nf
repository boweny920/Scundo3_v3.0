include {copy_index_files; target_tsv_lims_json; target_tsv_lims_json_PublicData;run_RmarkDown_script; run_RmarkDown_script_chip; biotools_orderAppend} from "../modules/report_gen.nf"

workflow FINAL_REPORT {
    take:
    secundo_output_rsem
    secundo_output_ercc
    secundo_output_picard
    samplesheet
    meta

    main:

    def single_meta_channel = meta
                                .collect()
                                .map{it[0]}

    copy_index_files(samplesheet, single_meta_channel, secundo_output_picard, secundo_output_rsem, secundo_output_ercc)   

    def secundoFolder= copy_index_files.out
                    .collect()
                    .map{it[0]}
    
    if (!params.public_dataxlsx) {
        target_tsv_lims_json(secundoFolder, samplesheet)

        if (!params.skip_biotools) { // If one does not wish data be published on biotools
            biotools_orderAppend(target_tsv_lims_json.out)
        }
        
        if (!params.chip_cut_atac) {
            run_RmarkDown_script(target_tsv_lims_json.out, single_meta_channel)
            report = run_RmarkDown_script.out 
        } else {
            run_RmarkDown_script_chip(target_tsv_lims_json.out, single_meta_channel)
            report = run_RmarkDown_script_chip.out
        }
        
    } else { // To account for Public SRA data analysis
        target_tsv_lims_json_PublicData(secundoFolder, samplesheet)
        biotools_orderAppend(target_tsv_lims_json_PublicData.out)
        run_RmarkDown_script(target_tsv_lims_json_PublicData.out, single_meta_channel)
        report = run_RmarkDown_script.out  
    }
    
    emit:
    report
}
