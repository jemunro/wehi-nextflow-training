
workflow {
    Channel.fromPath("$projectDir/logos.csv") |
        splitCsv(header: true) |
        view
}





