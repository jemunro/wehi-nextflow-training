
greeting = 'Hello'

process greet {
    output: stdout
    script: "echo -n $greeting world"
}

workflow {
    greet | view
}
