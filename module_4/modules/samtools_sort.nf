
process SAMTOOLS_SORT {
    // TODO: set cpus to 2
    // TODO: set memory to 2 Gigabytes
    // TODO: set time to 1 hours
    // TODO: tag with sample name
    // TODO: provide samtools through a module

    input:
    tuple val(sample), path(input_bam)

    // TODO: output files sample, sorted_bam and bam_index

    script:
    sorted_bam = sample + '.sorted.bam'
    bam_index = sorted_bam + '.bai'
    """
    samtools sort --threads $task.cpus $input_bam > $sorted_bam
    samtools index $sorted_bam
    """
}
