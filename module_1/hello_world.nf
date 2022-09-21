
greeting = 'Hello world'

process greet {
    output: stdout
    script: "echo -n $greeting"
}

workflow {
    greet() | view
}
