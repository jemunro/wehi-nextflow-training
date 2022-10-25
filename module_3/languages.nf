
workflow {

    year_created_ch = Channel.fromPath("$projectDir/data/year_created.csv", checkIfExists: true)
    year_created_ch.view()
}
