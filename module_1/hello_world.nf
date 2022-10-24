
process GREET {
    input: val(x)
    output: stdout
    script: "echo -n Hello $x"
}

workflow {
    input_ch = Channel.of('world')
    greet_ch = GREET(input_ch)
    greet_ch.view()
}
