
params.transcriptome_file = "$projectDir/data/transcriptome.fa.gz"

process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index -t $transcriptome -i salmon_index
    """
}

workflow {
    
    transcriptome_file = file(params.transcriptome_file, checkIfExists: true)

    index_ch = INDEX(transcriptome_file)
}