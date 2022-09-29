
process index_ref {
    cpus 1
    memory '1 GB'
    time '1 h'
    // TODO: provide bwa and samtools through modules
    // see https://www.nextflow.io/docs/latest/process.html#module

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
    // TODO: set cpus    https://www.nextflow.io/docs/latest/process.html#cpus
    // TODO: set memory  https://www.nextflow.io/docs/latest/process.html#memory
    // TODO: set time    https://www.nextflow.io/docs/latest/process.html#time
    // TODO: provide bwa and samtools through modules
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

process samtools_sort {
    cpus 2
    memory '2 GB'
    time '2 h'
    //TODO: provide samtools through a modules
    tag { sample }

    input:
    tuple val(sample), path(input_bam)

    output:
    // TODO: output sample, sorted_bam and bam_index

    script:
    sorted_bam = sample + '.sorted.bam'
    bam_index = sorted_bam + '.bai'
    """
    samtools sort -t 2 $input_bam > $sorted_bam
    samtools index $sorted_bam
    """
}


process bcftools_call {
    cpus 2
    memory '2 GB'
    time '2 h'
    //TODO: provide bcftools through a module
    tag { sample }

    input:
    tuple val(sample), path(sorted_bam), path(bam_index)
    tuple path(ref_fasta), path(ref_indices)

    output:
    path(bcf)

    script:
    bcf =  sample + '.bcf'
    //TODO: write bash script for variant calling
    // see https://samtools.github.io/bcftools/howtos/variant-calling.html
}


process bcftools_merge {
    cpus 2
    memory '2 GB'
    time '2 h'
    //TODO: provide bcftools through a module or conda

    input:
    path(bcf)

    output:
    //TODO: add output
    
    script:
    merged_vcf = 'merged.vcf.gz'
    """
    bcftools merge --no-index --missing-to-ref -Oz $bcf > $merged_vcf
    """
}

process plot_variants {
    cpus 1
    memory '2 GB'
    time '1 h'
    container 'library://jemunro/training/tidyverse-pheatmap'
    // TODO publish output using publishDir directive
    // see https://www.nextflow.io/docs/latest/process.html#publishdir

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
