
params.ref_fasta_gz    = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz'
params.sample_manifest = "$projectDir/vic_sarscov2_2020.tsv"

include { DOWNLOAD        } from './download'
include { DOWNLOAD_FASTQS } from './download'
// include { INDEX_REF       } from './processes'
// include { BWA_MEM_ALIGN   } from './processes'
// include { SAMTOOLS_SORT   } from './processes'
// include { BCFTOOLS_CALL   } from './processes'
// include { BCFTOOLS_MERGE  } from './processes'
// include { PLOT_VARIANTS   } from './processes'

workflow {

    ref_files_ch = DOWNLOAD(params.ref_fasta_gz)

    fastq_url_ch = Channel
        .fromPath(params.sample_manifest)
        .splitCsv(sep: '\t', header: true)
        .map { [it.accession, it.fastq1, it.fastq2] }
    
    fastq_ch = DOWNLOAD_FASTQS(fastq_url_ch)

    // index_ch = INDEX_REF(ref_files_ch)
    // aligned_ch = BWA_MEM_ALIGN(fastq_ch, index_ch)
    // sorted_ch = SAMTOOLS_SORT(aligned_ch)
    // called_ch = BCFTOOLS_CALL(sorted_ch, index_ch)
    // merged_ch = BCFTOOLS_MERGE(called.collect())
    // plot_ch = PLOT_VARIANTS(merged, file(params.sample_manifest))
    
    // plot_ch.view()
}
