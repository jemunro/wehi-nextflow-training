
params.ref_fasta_gz    = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz'
params.sample_manifest = "$projectDir/vic_sarscov2_2020.tsv"

include { index_ref      } from './processes'
include { bwa_mem_align  } from './processes'
// include { samtools_sort  } from './processes'
// include { bcftools_call  } from './processes'
// include { bcftools_merge } from './processes'
include { plot_variants  } from './processes'


workflow {
    
    ref_files = index_ref(file(params.ref_fasta_gz))

    sample_fastqs = Channel.fromPath(params.sample_manifest) |
        splitCsv(sep: '\t', header: true) |
        map { [it.accession, file(it.fastq1), file(it.fastq2)] }

    aligned = bwa_mem_align(sample_fastqs, ref_files)

    // sorted = samtools_sort(aligned)
    // called = bcftools_call(sorted)
    // merged = bcftools_merge(called.toList())


}