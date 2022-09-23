
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
e.g. samtools sort --threads 2 $input_bam > $sorted_bam
     samtools index $sorted_bam
inputs from bwa_mem_align
output: sample, sorted_bam, bam_index
*/

process samtools_sort {
    cpus 2
    memory '2 GB'
    time '2 h'
    module 'samtools/1.16.1'
    tag { sample }

    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path(sorted_bam), path(bam_index)

    script:
    sorted_bam = sample + '.sorted.bam'
    bam_index = sorted_bam + '.bai'
    """
    samtools sort --threads 2 $input_bam > $sorted_bam
    samtools index $sorted_bam
    """
}

/*
TODO: process bcftools_call
call variants using bcftools
find an apropriate module or bioconda pacakge
see https://samtools.github.io/bcftools/howtos/variant-calling.html
input from samtools_sort
input from index_ref
output: vcf
*/

process bcftools_call {
    cpus 2
    memory '2 GB'
    time '2 h'
    module 'bcftools/1.16'
    tag { sample }

    input:
    tuple val(sample), path(sorted_bam), path(bam_index)
    tuple path(ref_fasta), path(ref_indices)

    output:
    path(vcf)

    script:
    vcf =  sample + 'vcf.gz'
    """
    bcftools mpileup -Ou -f $ref_fasta $sorted_bam | bcftools call -mv -Oz -o $vcf
    """
}

/*
TODO: process bcftools merge
http://samtools.github.io/bcftools/bcftools.html#merge
merge variant calls
*/

process bcftools_merge {
    cpus 2
    memory '2 GB'
    time '2 h'
    module 'bcftools/1.16'

    input:
    path(vcfs)

    output:
    path(merged_vcf)
    
    script:
    merged_vcf = 'merged.vcf.gz'
    """
    bcftools merge --no-index --missing-to-ref -Oz $vcfs > $merged_vcf
    """
}

process plot_variants {
    cpus 1
    memory '2 GB'
    time '1 h'
    container 'library://jemunro/training/tidyverse-pheatmap'
    publishDir "output", mode: 'copy'

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
