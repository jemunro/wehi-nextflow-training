
process Greet {
    input: val(x)
    output: stdout
    script: "echo -n Hello $x"
}

workflow {
    input_ch = Channel.of('world')
    Greet(input_ch)
    Greet.out.view()
}
