process TOUPPER {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | tr a-z A-Z"
}

process REVERSE {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | rev"
}

process CONCAT {
    input: tuple val(x), val(y), val(z)
    output: stdout
    script: "echo -n $x $y $z"
}

workflow {
    input = channel.of('foo', 'bar', 'baz')
    upper = TOUPPER(input)
    reveresed = REVERSE(input)
    output = CONCAT(upper.join(reveresed))
    output.view()
}

