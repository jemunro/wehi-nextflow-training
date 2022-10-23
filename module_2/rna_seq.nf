
params.transcriptome_file = "$projectDir/data/transcriptome.fa.gz"
// params.reads = "$projectDir/data/gut_{1,2}.fq.gz"
params.reads = "$projectDir/data/*_{1,2}.fq.gz"
params.outdir = "results"

process INDEX {
    module 'salmon/1.9.0'

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
    tag "$sample_id"
    // module 'salmon/1.9.0'
    container 'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1'


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

process PLOT_TPM {
    container 'rocker/tidyverse:4.1.3'

    input:
    path quant_results

    output:
    path 'TPM.png'

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    data.frame(filename = list.files(pattern='.sf')) %>% 
        mutate(tissue = str_remove(basename(filename), '.sf')) %>% 
        mutate(data = map(filename, read_tsv, col_types = cols())) %>% 
        unnest(data) %>% 
        select(tissue, transcript = Name, TPM) %>% 
        ggplot(aes(transcript, TPM, fill = tissue)) +
        geom_col(position = 'dodge') +
        coord_flip()

    ggsave('TPM.png', width = 6, height = 4)
    """
}


workflow {

    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // index_ch = INDEX(params.transcriptome_file)

    // quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    // plot_ch = PLOT_TPM(quant_ch.collect())

    // plot_ch.view()
}
