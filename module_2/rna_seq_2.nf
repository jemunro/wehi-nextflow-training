
params.transcriptome_file = "$projectDir/data/transcriptome.fa.gz"
params.reads = "$projectDir/data/*_{1,2}.fq.gz"

process INDEX {
    module 'salmon/1.9.0'
    memory '2 GB'
    cpus 1

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    module 'salmon/1.9.0'
    memory '2 GB'
    cpus 2
    tag "$sample_id"

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path output

    script:
    output = "${sample_id}.sf"
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o out
    mv out/quant.sf $output
    """
}

workflow {

    transcriptome_file = file(params.transcriptome_file, checkIfExists: true)
    
    index_ch = INDEX(transcriptome_file)

    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    quant_ch.view()
}
