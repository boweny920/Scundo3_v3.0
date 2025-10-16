include {SAMPLESHEET_CHANNEL} from "./workflows/samplesheet_channel.nf"
include {SECUNDO3} from "./workflows/secundo3.nf"
include {FINAL_REPORT} from "./workflows/final_report.nf"
include {CHIP_RELATED} from "./workflows/cut_chip_atac.nf"
include {samplesheet_process} from "./modules/core.nf"
include {run_RmarkDown_script} from "./modules/report_gen.nf"
include {samplesheet_process_publicSRAGEO} from "./modules/public_SRAGEO.nf"

def report

workflow {
        if (!params.public_dataxlsx) {
            if (!params.fcid) {
                exit 1, "ERROR: Provide a Flowcell ID for Secundo3."
            } else { 
                samplesheet_process(params.fcid)
                SAMPLESHEET_CHANNEL(samplesheet_process.out)

                // RNA seq related Block
                if (!params.chip_cut_atac) {
                    SECUNDO3(SAMPLESHEET_CHANNEL.out.data_meta)
                    FINAL_REPORT(SECUNDO3.out.rsem, SECUNDO3.out.ercc, SECUNDO3.out.picard, samplesheet_process.out, SAMPLESHEET_CHANNEL.out.data_meta)
                    report = FINAL_REPORT.out

                // Chip Cut&RUN/TAG ATAC seq Block
                } else {
                    CHIP_RELATED(SAMPLESHEET_CHANNEL.out.data_meta)
                    FINAL_REPORT(CHIP_RELATED.out.scundo_dir, CHIP_RELATED.out.correlation, CHIP_RELATED.out.coverage, samplesheet_process.out, SAMPLESHEET_CHANNEL.out.data_meta )

                    // FINAL_REPORT(CHIP_RELATED.out.scundo_dir, CHIP_RELATED.out.correlation, CHIP_RELATED.out.coverage, samplesheet_process.out, SAMPLESHEET_CHANNEL.out.data_meta )
                    report = FINAL_REPORT.out
                }            
            }
            
        } else { // Pubulic SRA/GEO RNA-seq analysis block   
            if (!params.requester) {
                exit 1, "ERROR: Provide a requester ID for public data analysis."
            } else if (!params.lab) {
                exit 1, "ERROR: Provide a Lab ID for public data analysis."
            } else {
                samplesheet_process_publicSRAGEO(params.public_dataxlsx)
                SAMPLESHEET_CHANNEL(samplesheet_process_publicSRAGEO.out) 
                SECUNDO3(SAMPLESHEET_CHANNEL.out.data_meta)
                FINAL_REPORT(SECUNDO3.out.rsem, SECUNDO3.out.ercc, SECUNDO3.out.picard, samplesheet_process_publicSRAGEO.out, SAMPLESHEET_CHANNEL.out.data_meta)
                report = FINAL_REPORT.out 
            }
        
        }
    }


workflow.onComplete {

    def msg = """\
        Secundo3 Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        workDir     : ${workflow.workDir}
        Job status : ${ workflow.success ? "Success! SECUNDO run Report link : ${report.value}" : "failed, error message: ${workflow.errorMessage}" }
        """
        .stripIndent()
    
    println msg 
    sendMail(to: "${params.user_email}", subject: "${params.fcid} Secundo3 pipeline execution finished", body: msg)
}
