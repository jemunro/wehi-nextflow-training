
process jp2a {
   cpus 1
   memory '1 GB'
   time '1 h'
   container 'talinx/jp2a' 

   input:
   tuple val(name), val(image)

   output:
   tuple val(name), stdout

   script:
   """
   jp2a $image --colors --width=35 --color-depth=24
   """
}

workflow {

    languages = Channel.fromPath("$projectDir/languages.csv") |
        splitCsv(skip: 1)

    // languages | view

    logos = Channel.fromPath("$projectDir/logos.csv") |
        splitCsv(skip: 1) |
        jp2a

    languages |
        join(logos) |
        view { name, year, author, logo -> 
            "$name was created in $year by $author:\n$logo"
        }

}

