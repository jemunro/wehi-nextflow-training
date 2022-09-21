
message = 'Hello world'

process greet {
    output: stdout
    script: "echo -n $message"
}

workflow {
    greet() | view
}
