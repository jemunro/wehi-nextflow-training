
process download {
    cpus 1
    memory '1 GB'
    time '1 h'
    executor 'local'

    input:
    val(url)

    output:
    path(output)

    script:
    output = file(url).name
    """
    wget $url
    """
}

process download_fastqs {
    cpus 1
    memory '1 GB'
    time '1 h'
    executor 'local'

    input:
    tuple val(sample), val(url1), val(url2)

    output:
    tuple val(sample), path(fastq1), path(fastq2)

    script:
    fastq1 = file(url1).name
    fastq2 = file(url2).name
    """
    wget $url1
    wget $url2
    """
}
