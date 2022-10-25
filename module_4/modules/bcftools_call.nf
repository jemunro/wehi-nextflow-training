
process BCFTOOLS_CALL {
    // TODO: set cpus to 2
    // TODO: set memory to 2 Gigabytes
    // TODO: set time to 1 hours
    // TODO: tag with sample name
    // TODO: provide bcftools through container using tag '1.16--hfe4b78e_1'
    //       see https://bioconda.github.io/recipes/bcftools/README.html?highlight=bcftools

    input:
    tuple val(sample), path(sorted_bam), path(bam_index)
    tuple path(ref_fasta), path(ref_indices)

    // TODO: output bcf file

    script:
    bcf =  sample + '.bcf'
    //TODO: write bash script for variant calling
    // see https://samtools.github.io/bcftools/howtos/variant-calling.html
}
