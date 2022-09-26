
audience = ['world']

process greet {
    input: val(x)
    output: stdout
    script: "echo -n Hello $x!"
}

workflow {
    channel.fromList(audience) |
        greet |
        view
}