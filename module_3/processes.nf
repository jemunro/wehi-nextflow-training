
process index_ref {
    cpus 1
    memory '1 GB'
    time '1 h'
    module 'bwa/0.7.17'
    module 'samtools/1.16.1'

    input:
    path(ref_fasta_gz)

    output:
    tuple path(ref_fasta), path("${ref_fasta}.*")

    script:
    ref_fasta = 'ASM985889v3_genomic.fasta'
    """
    gzip -cd $ref_fasta_gz > $ref_fasta
    samtools faidx $ref_fasta
    bwa index $ref_fasta
    """
}

process bwa_mem_align {
    cpus 2
    memory '2 GB'
    time '2 h'
    module 'bwa/0.7.17'
    module 'samtools/1.16.1'
    tag { sample }

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    tuple path(ref_fasta), path(ref_indices)

    output:
    tuple val(sample), path(bam)

    script:
    bam = sample + '.bam'
    """
    bwa mem -M -t 2 -R '@RG\\tID:$sample\\tSM:$sample' $ref_fasta $fastq1 $fastq2 |
        samtools view -b > $bam
    """
}

/*
TODO: process samtools_sort
sort and index bam file
http://www.htslib.org/doc/samtools-sort.html
e.g. samtools sort -t 2 -b $input_bam > $sorted_bam
inputs from bwa_mem_align
output: sample, sorted_bam, bam_index
*/

/*
TODO: process bcftools_call
call variants using bcftools
see https://samtools.github.io/bcftools/howtos/variant-calling.html
input from samtools_sort
input from index_ref
output: vcf
*/

/*
TODO: process bcftools merge
http://samtools.github.io/bcftools/bcftools.html#merge
merge variant calls
*/

process plot_variants {
    cpus 2
    memory '2 GB'
    time '2 h'
    container 'rocker/tidyverse:4.2.1'
    publishDir "output", mode: 'copy'
    tag { sample }

    input:
    path(vcf)

    output:
    path(plot)

    script:
    plot = 'plot.png'
    """
    plot_variants.R $vcf $plot
    """
}
