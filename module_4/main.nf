
params.ref_fasta_gz    = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz'
params.sample_manifest = "$projectDir/vic_sarscov2_2020.tsv"

include { DOWNLOAD_REF    } from './modules/download'
include { DOWNLOAD_FASTQS } from './modules/download'
// include { INDEX_REF       } from './modules/index_ref'
// include { BWA_MEM_ALIGN   } from './modules/bwa_mem_align'
// include { SAMTOOLS_SORT   } from './modules/samtools_sort'
// include { BCFTOOLS_CALL   } from './modules/bcftools_call'
// include { BCFTOOLS_MERGE  } from './modules/bcftools_merge'
// include { PLOT_VARIANTS   } from './modules/plot_variants'

workflow {

    ref_files_ch = DOWNLOAD_REF(params.ref_fasta_gz)

    fastq_url_ch = Channel
        .fromPath(params.sample_manifest, checkIfExists: true)
        .splitCsv(sep: '\t', header: true)
        .map { [it.accession, it.fastq1, it.fastq2] }
    
//     fastq_ch = DOWNLOAD_FASTQS(fastq_url_ch)

//     index_ch = INDEX_REF(ref_files_ch)

//     aligned_ch = BWA_MEM_ALIGN(fastq_ch, index_ch)

//     sorted_ch = SAMTOOLS_SORT(aligned_ch)

//     called_ch = BCFTOOLS_CALL(sorted_ch, index_ch)

//     merged_ch = BCFTOOLS_MERGE(called_ch.collect())

//     PLOT_VARIANTS(merged_ch, file(params.sample_manifest, checkIfExists: true))
}

// workflow.onComplete {
//     log.info( 
//         workflow.success ? 
//         "Pipeline Complete!\nOutput: $launchDir/results/plot.png\n" : 
//         "Pipeline Failed.\n" )
// }