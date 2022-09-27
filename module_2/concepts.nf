
workflow {
    languages = Channel.fromPath("$projectDir/languages.csv") |
        splitCsv(skip: 1) |
        view`
}


