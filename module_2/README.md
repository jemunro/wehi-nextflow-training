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
   ```
* **Lists**
   * groovy:
      ```groovy
      l = [1, 2, 3]
      l = l + 4
      l = l.collect { it * 2 }
      print(l)
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
   * groovy:
      ```groovy
      x = [foo: 1, bar: 2]
      println x.foo
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
### **Exercise 1.1**
* try out the above  groovy examples at https://groovyconsole.appspot.com/

## 2. Channels

## 3. Operators

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