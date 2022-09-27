
process stringUpper {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | tr a-z A-Z"
}

process stringReverse {
    input: val(x)
    output: tuple val(x), stdout
    script: "echo -n $x | rev"
}

process stringConcat {
    input: tuple val(x), val(y), val(z)
    output: stdout
    script: "echo -n $x-$y-$z"
}

workflow {
    input = channel.of('foo', 'bar', 'baz')
    
    input |
        toList() | 
        view { "input: $it" }

    stringUpper(input) |
        join(stringReverse(input)) |
        stringConcat |
        toSortedList() |
        view { "output: $it" }
}


