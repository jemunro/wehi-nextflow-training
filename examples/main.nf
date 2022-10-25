process TOUPPER {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | tr a-z A-Z"
}

process REVERSE {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo $x | rev | tr -d '\\n'"
}

workflow {
    input_ch = channel.of('foo')
    upper_ch = TOUPPER(input_ch)
    reverese_ch = REVERSE(input_ch)
    output_ch = upper_ch.join(reverese_ch)
    output_ch.view()
}




