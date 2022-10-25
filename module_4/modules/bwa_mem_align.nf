
process BWA_MEM_ALIGN {
    // TODO: set cpus to 2
    // TODO: set memory to 4 Gigabytes
    // TODO: set time to 2 hours
    // TODO: tag with sample name
    // TODO: provide bwa and samtools through modules

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    tuple path(ref_fasta), path(ref_indices)

    output:
    tuple val(sample), path(bam)

    script:
    bam = sample + '.bam'
    """
    bwa mem -M -t $task.cpus -R '@RG\\tID:$sample\\tSM:$sample' $ref_fasta $fastq1 $fastq2 |
        samtools view -b > $bam
    """
}
