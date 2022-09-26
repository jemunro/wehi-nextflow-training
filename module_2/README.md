# Module 2: Nextflow Concepts

### Learning Objectives
1. Nextflow/Groovy scripting basics
2. Understand channel creation
3. Use operators to transform channels

## 1. Scripting
### 1.1 Groovy Scripting
* Much of Nextflow scripting is done using Groovy
   ```groovy
   x = 1                    // integer
   y = 1.5                  // float
   println(x + y)           // print statement, addition
   println x * y            // print statement, optional parenthesis, multiplication
   s1 = 'foo'               // String
   s2 = 'bar'               // String
   println s1 + '-' + s2    // String concatenation
   println "$s1-$s2"        // Variable interpolation
   b = true                 // boolean
   println( b ? 1 : 0 )     // ternary operator (if_else)
   ```
* **Lists**
   * Nextflow/groovy:
      ```groovy
      l = [1, 2, 3]
      l = l + 4
      l = l.collect { it * 2 } // <-- closure
      println(l)
      ```
   * R:
      ```R
      l = list(1, 2, 3)
      l = c(l, 4)
      l = lapply(l, function(x) x * 2)
      print(l)
      ```
   * python:
      ```python
      l = [1, 2, 3]
      l.append(4)
      l = map(lambda x: x * 2, l)
      print(l)
      ```
* **Maps**
   * Equivalent to named lists in `R` or dicts in `python`
   * Nextflow/groovy:
      ```groovy
      x = [foo: 1, bar: 2]
      println(x.bar)
      x.baz = 3
      ```
   * R:
      ```R
      x = list(foo = 1, bar = 2)
      print(x$bar)
      x$baz = 3
      ```
   * python:
      ```python
      x = dict(foo = 1, bar = 2)
      print(x['bar'])
      x['baz'] = 3
      ```

* **Optional Function Arguments**
   * In `R` and `python`, optional arguments are given to a function using `=`, e.g.:
   ```R
   myFunction(foo = bar)
   ```
   * In Nextflow/Groovy, optional arguments are provided using `:`, e.g.:
   ```groovy
   myFunction(foo: bar)
   ```
### **Exercise 2.1**
1. Try out some of above groovy examples using the groovy shell:
   ```
   module load java/1.8.0_92 groovy/4.0.0
   groovysh
   ```
   or at https://groovyconsole.appspot.com/


### 1.2 Nextflow Scripting
**Implicit Variables**:
* A number of variables are available in all nextflow scripts:
* `params`: map storing workflow parameters
* `launchDir`: The directory the nextflow script was run from
* `projectDir`: The directory containing the main nextflow script
* See https://www.nextflow.io/docs/latest/script.html#implicit-variables

**Files**:
   * For creating `file` type variables, nextflow defines the function `file()`
   * This must be used to exlicity convert a filepath string into a file object before passing to a process
   * We can use the `checkIfExists` optional argument to make sure the file exists:
      ```groovy
      input = file('/path/to/my_file.txt', checkIfExists: true)
      ```
   * Nextflow files also support the  FTP/HTTP protocols. These will be downloaded automatically when used by a process:
      ```groovy
      input = file('https://www.wehi.edu.au/sites/default/files/wehi-logo-2020.png')
      ```



## 2. Channels & Operators
### **Channel Creation**:
* Nextflow includes a number of ways to create channels
* `channel.of(...)`: create a channel that emits each of the arguments one at a time.
   ```groovy
   channel.of('A', 'B') | view
   ```
* `channel.fromList(list)`: given a list, create a channel the emits each element of the list one at a time
   ```groovy
   channel.fromList(['A', 'B']) | view
   ```
* `channel.fromPath(path)`: Create a channel from a file path or glob pattern, emitting a `file` object.
   ```groovy
   channel.fromPath('/path/to/sample.bam') | view
   ```
* See https://www.nextflow.io/docs/latest/channel.html#channel-factory
### **Operators**
### Map
* `map` is the most commonly used nextflow operater
* Functionally similar to R's `lapply()` or Python's `map()`functions in R
* The default variable name for a closure is `it`
* for example:
   ```groovy
   channel.of(1, 2, 3) |
      map { it * 2 } |  
      view()
   ```
   will print:
   ```
   2
   4
   6
   ```
### splitCsv
* Convert a CSV (or TSV) file into a nextflow channel
* Optional argument `skip` may be used to skip a number of lines (e.g. the header)
* Most often used on the output of `channel.fromPath()`, e.g.
   ```groovy
   channel.fromPath("$projectDir/logos.csv") |
        splitCsv(skip: 1) |
        view
   ```
### join
* Join two channels by a matching key
   ```groovy
   left  = channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
   right = channel.of(['Z', 6], ['Y', 5], ['X', 4])
   left | join(right) | view()
   ```
   ```
   [Z, 3, 6]
   [Y, 2, 5]
   [X, 1, 4]
   ```
* This operation is roughly equivalent to Rs `dplyr::inner_join()` and Python's pandas `pd.merge()`. For example, in R:
   ```R
   left  = data.frame(V1 = c('X', 'Y', 'Z', 'P'),
                      V2 = c( 1,   2,   3,   7))
   right = data.frame(V1 = c('Z', 'Y', 'X'),
                      V3 = c( 6,   5,   4))
   left %>% inner_join(right)
   ```
   ```
     V1 V2 V3
   1  X  1  4
   2  Y  2  5
   3  Z  3  6
   ```
### More
* Many more operators are availabe, see https://www.nextflow.io/docs/latest/operator.html
### **Exercise 2.2**

1. Open [concepts.nf](concepts.nf), [languages.csv](languages.csv) and [logos.csv](logos.csv)
1. Run `concepts.nf` and observe the output
   ```
   nextflow run ~/wehi-nextflow-training/module_2/concepts.nf
   ```
1. Create a channel named `logos` from [logos.csv](logos.csv) as is done with `languages`
1. Use the `join()` operator to join `languages` with `logos`, and print the result with `view()`

## 3. Processes
### Inputs
* **val**
* **path**
* **tuple**
* **multiple inputs**

```groovy
process jp2a {
   cpus 1                  //
   memory '1 GB'           // (1) process directives
   time '1 h'              // 
   container 'talinx/jp2a' //https://hub.docker.com/r/talinx/jp2a

   input:
   tuple val(name), path(image) // (2) tuple input

   output:
   tuple val(name), stdout

   script:
   """ // (5) Multiline string script
   jp2a https://www.wehi.edu.au/sites/default/files/wehi-logo-2020.png --colors --width=100
   """ 
}
```


## 5. Configutation
* look at `nextflow.config` and ~/.nextflow/config


## 2. Processes
```nextflow
process sort_bam {
   cpus 2                  //
   memory '2 GB'           // (1) process directives
   time '1 h'              // 
   module 'samtools/1.15'  //

   input:
   tuple val(sample), path(unsorted_bam) // (2) tuple input

   output:
   tuple val(sample), path(sorted_bam), path(index) // (3) tuple output

   script:
   sorted_bam = sample + '.sorted.bam'  // (4) defining task variables
   index = sorted_bam + '.bai'
   """ // (5) Multiline string script
   samtools sort $unsorted_bam -o $sorted_bam \\
      --output-fmt BAM \\
      --threads 2 \\
      --write-index
   """ 
}
```

### **Exercise XX**
* Write a process using a different aligner