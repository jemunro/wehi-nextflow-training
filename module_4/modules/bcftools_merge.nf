
process BCFTOOLS_MERGE {
    // TODO: set cpus to 2
    // TODO: set memory to 2 Gigabytes
    // TODO: set time to 1 hours
    // TODO: provide bcftools through container using tag '1.16--hfe4b78e_1'
    //       see https://bioconda.github.io/recipes/bcftools/README.html?highlight=bcftools

    input:
    path(bcfs)

    //TODO: output merged vcf file
    
    script:
    merged_vcf = 'merged.vcf.gz'
    """
    bcftools merge --threads $task.cpus --no-index --missing-to-ref -Oz $bcfs > $merged_vcf
    """
}
