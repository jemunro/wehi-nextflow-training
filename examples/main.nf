process ToUpper {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | tr a-z A-Z"
}

process Reverse {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | rev"
}

process Concat {
    input: tuple val(x), val(y), val(z)
    output: stdout
    script: "echo -n $x $y $z"
}

workflow {
    input = channel.of('foo', 'bar', 'baz')
    upper = ToUpper(input)
    reveresed = Reverse(input)
    output = Concat(upper.join(reveresed))
    output.view()
}

