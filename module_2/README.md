# Module 2: Nextflow Concepts

### Learning Objectives
1. Understand nextflow process directives
1. Manage software using `module`
2. Manage software using `conda`
3. Managing software using `container`

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
      l = lapply(l, function(x) x + 1)
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

### 1.2 Nextflow Scripting
* **Implicit Variables**:
* **Files**:
   * For creating `file` type variables, nextflow defines the function `file()`

### **Exercise 1.1** ?
* try out the above  groovy examples at https://groovyconsole.appspot.com/
* ..?

## 2. Channels & Operators
### Channel Creation:
* Nextflow includes a number of ways to create channels
* The most importand methods are:
   * `channel.of(...)`: create a channel that emits each of the arguments one at a time.
      ```nextflow
      channel.of('A', 'B') | view
      ```
   * `channel.fromList(list)`: given a list, create a channel the emits each element of the list one at a time
      ```nextflow
      channel.fromList(['A', 'B']) | view
      ```
   * `channel.fromPath(path)`: Create a channel from a file path or glob pattern, emitting a `file` object.
      ```nextflow
      channel.fromPath('/path/to/sample.bam') | view
      ```
### **Operators**
### Map
* `map` is the most commonly used nextflow operater
* Functionally similar to R's `lapply()` or Python's `map()`functions in R
* The default variable name for a closure is `it`
* for example:
   ```nextflow
   channel.of(1, 2, 3) |
      map { it * 2 } | // "it" is the 
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
* Most often used on the output of `channel.fromPath()`, e.g.
   ```nextflow
   channel.fromPath("$projectDir/logos.csv") |
        splitCsv(header: true) |
        view
   ```
### More
* Many more operators are availabe, see https://www.nextflow.io/docs/latest/operator.html
### **Exercise 2.X**
1. Look at [`concepts.nf`](concepts.nf) and [logos.csv](`logos.csv`)
1. Run `concepts.nf` and observe the output
   ```
   nextflow run ~/wehi-nextflow-training/module_2/concepts.nf
   ```
1. Add the following map statement between `splitCsv` and `view` in [`concepts.nf`](concepts.nf) and run it again:
   ```nextflow
   map { [it.name, file(it.url)] } 
   ```


 see the [channel factory documentation](https://www.nextflow.io/docs/latest/channel.html#channel-factory)


## 4. Processes

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