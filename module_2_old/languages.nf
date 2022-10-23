
workflow {
    details = Channel.fromPath("$projectDir/details.csv", checkIfExists: true)
    details.view()
}
