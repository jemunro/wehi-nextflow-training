
process INDEX_REF {
    // TODO: set cpus to 1 - https://www.nextflow.io/docs/latest/process.html#cpus
    // TODO: set memory to 1 Gigabyte - https://www.nextflow.io/docs/latest/process.html#memory
    // TODO: set time to 1 hours - https://www.nextflow.io/docs/latest/process.html#time
    // TODO: provide bwa and samtools through modules - https://www.nextflow.io/docs/latest/process.html#module

    input:
    path(ref_fasta_gz)

    output:
    tuple path('ref.fasta'), path("ref.fasta.*")

    script:
    """
    gzip -cd $ref_fasta_gz > ref.fasta
    samtools faidx ref.fasta
    bwa index ref.fasta
    """
}
