include {bowtie2_samtools_n_mkdirs; deep_tools_correlation; deep_tools_coverage} from "../modules/chip_related.nf"


workflow CHIP_RELATED {
    take:
    meta

    main:
    bowtie2_samtools_n_mkdirs(meta)

    def scundo_outDir = bowtie2_samtools_n_mkdirs.out
                                    .collect()
                                    .map {it[0]}

    // deep_tools_correlation(scundo_outDir)
    // deep_tools_coverage(scundo_outDir)

    emit:
    scundo_dir = scundo_outDir
    // correlation = deep_tools_correlation.out.npz
    correlation = "PlaceHolder"
    // coverage = deep_tools_coverage.out.tsv
    coverage = "PlaceHolder"

    
}