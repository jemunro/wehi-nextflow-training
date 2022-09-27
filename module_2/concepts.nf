

workflow {

    languages = Channel.fromPath("$projectDir/languages.csv") |
        splitCsv(skip: 1) |
        view


    logos = Channel.fromPath("$projectDir/logos.csv") |
        splitCsv(skip: 1) |
        download |
        jp2a

    languages |
        join(logos) |
        view { name, year, author, logo -> 
            "$name was created in $year by $author:\n$logo"
        }
}

process download {
    input: tuple val(name), val(url)

    output: tuple val(name), path(output)

    script:
    output = file(url).name
    """
    wget $url
    """
}

process jp2a {
   container 'talinx/jp2a' 
   input: tuple val(name), path(image)

   output: tuple val(name), stdout

   script:
   """
   jp2a $image --colors --width=35 --color-depth=24
   """
}

