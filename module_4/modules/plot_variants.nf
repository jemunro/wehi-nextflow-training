
process PLOT_VARIANTS {
    // TODO: set cpus to 1
    // TODO: set memory to 2 Gigabytes
    // TODO: set time to 1 hours
    // TODO: provide R dependencies though container 'library://jemunro/training/tidyverse-pheatmap'
    // TODO publish output to "results" using publishDir directive - see https://www.nextflow.io/docs/latest/process.html#publishdir

    input:
    path(vcf)
    path(metadata)

    output:
    path(plot)

    script:
    plot = 'plot.png'
    """
    plot_variants.R $vcf $metadata $plot
    """
}
