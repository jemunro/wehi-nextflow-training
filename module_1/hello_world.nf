
audience = ['world']

process Greet {
    input: val(x)
    output: stdout
    script: "echo -n Hello $x!"
}

workflow {
    channel.fromList(audience) |
        Greet |
        view
}