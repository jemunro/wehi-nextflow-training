process LocalHostName {
    executor 'local'
    output: stdout
    script: "hostname"
}

process SlurmHostName {
    executor 'slurm'
    output: stdout
    script: "hostname"
}

workflow {
    LocalHostName().view { "local: ${it.trim()}" }
    SlurmHostName().view { "slurm: ${it.trim()}" }
}
