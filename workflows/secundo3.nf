include {STAR_n_mkdirs; RSEM_calculate; ERCC_control; picard_run} from "../modules/core.nf"

workflow SECUNDO3 {
    take:
    meta

    main:
    STAR_n_mkdirs(meta)
    RSEM_calculate(STAR_n_mkdirs.out)
    ERCC_control(STAR_n_mkdirs.out)
    picard_run(STAR_n_mkdirs.out)

    emit:
    rsem = RSEM_calculate.out.collect()
    ercc = ERCC_control.out.collect()
    picard = picard_run.out.collect()
    
}