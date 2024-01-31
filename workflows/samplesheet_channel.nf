workflow SAMPLESHEET_CHANNEL {
    take:
    samplesheet

    main:

    // Channel.from( samplesheet )
    samplesheet.splitCsv ( header:true, sep:'\t' )
        .map { create_metadata_channel(it) }
        .set { data_meta }

    emit:
    data_meta
}

def create_metadata_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.ID = row.secundoname
    meta.samplename = row.samplename
    meta.scundoname = row.secundoname
    meta.genome_ver = row.genome_ver
    meta.species = row.species
    meta.fastqs = row.fastqs
    meta.molngID = row.molngID
    meta.pi_name = row.pi_name
    meta.requester_name = row.requester_name
    meta.readType = row.readType
    meta.annotation = row.annotation

    return meta
}