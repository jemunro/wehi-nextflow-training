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
      println(l[0]) // (0-based index)
      ```
   * R :
      ```R 
      l = list(1, 2, 3)
      l = c(l, 4)
      l = lapply(l, function(x) x * 2)
      print(l)
      print(l[1]) # (1-based index)
      ```
   * python (0-based):
      ```python 
      l = [1, 2, 3]
      l.append(4)
      l = map(lambda x: x * 2, l)
      print(l)
      print(l[0]) # (0-based index)
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

### **Exercise 2.1**
1. Try out some of above groovy examples using the groovy shell:
   ```
   module load java/1.8.0_92 groovy/4.0.0
   groovysh
   ```

### 1.2 Groovy Closures
* Closures in groovy act as functions that can be passed to other functions
* For example the `.collect()` function which is a property of lists in groovy
   ```groovy
   l = [1, 2, 3]
   println( l.collect { it + 1 } )
   ```
* The default variable name in a closure is `it`, but this can be overridden:
   ```groovy
   l = [1, 2, 3]
   l.collect { x -> x + 1 }
   ```
* In Groovy and Nextflow we often deal with nested lists
   ```groovy
   nested = [[1, 'A'], [2, 'B'], [3, 'C']]
   ```
* When we call `collect` on a nested list such as this, the variable `it` will be a list. For example here we are concatenating the elements of each inner list togeher:
   ```groovy
   nested.collect { "${it[0]}-${it[1]}" }
   ```
   ```
   [1-A, 2-B, 3-C]
   ```
* We can also unpack the values of the list and assign them their own variables, which makes writing closures clearer
   ```groovy
   nested.collect { number, letter -> "$number-$letter" }
   ```
* If we wanted to keep the `number` and `letter` elements, we could instead do this:
   ```groovy
   nested.collect { number, letter -> [number, letter, "$number-$letter"] }
   ```
   ```
   [[1, A, 1-A], [2, B, 2-B], [3, C, 3-C]]
   ```
* See https://www.nextflow.io/docs/latest/script.html#closures

### **Exercise 2.2**
1. In the Groovy shell, define the variables `data` as below
   ```groovy
   data = [['foo', 1, 2], ['bar', 3, 4], ['baz', 5, 6]]
   ```
1. Using `collect`, transform data to have a new value which is the product of the first and second numeric values in each list, e.g. `['bar', 3, 4]` -> `['bar', 3, 4, 12]`

### 1.3 Nextflow Scripting
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
* `map` is the most commonly used nextflow operater, and works the same as Groovy's `collect()` but applied to channels instead of lists.
* Functionally similar to R's `lapply()` or Python's `map()` 
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
### View
* `view` takes an input channel and prints the contents to the terminal
* We can optionally provide a closure, similarly to `map {}`, which will instead print the value of applying the closure
   ```groovy
   channel.of(1, 2, 3) |  view { it * 2 } 
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
   channel.fromPath("$projectDir/languages.csv") |
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

* Many more operators are availabe, see https://www.nextflow.io/docs/latest/operator.html

### **Exercise 2.3**
1. Open [concepts.nf](concepts.nf), [languages.csv](languages.csv) and [logos.csv](logos.csv)
1. Run `concepts.nf` and observe the output
   ```
   nextflow run ~/wehi-nextflow-training/module_2/concepts.nf
   ```
1. Create a channel named `logos` from [logos.csv](logos.csv) as is done with `languages`
1. Use the `join()` operator to join `languages` with `logos`, and print the result with `view()`

## 3. Processes
### **Inputs**
* `val()` - A val type input denotes a regular groovy variable. It could be a String, Integer, Boolean, double etc.
* `path()` - A path represents an input file.
* `tuple` - A tuple represents a list of inputs. There may be of either `val` or `path` types
### **Outputs**
* similarly to inputs, outputs may be of type `val`, `path` or `tuple`
* `path()` - when a `path` is declared as an output, after the process has run sucessfully nextflow will check that the path exists, and throw an error if not.
* `stdout` - stdout is a special output type that will return the standard output of the process run

### **Exercise 2.4**
1. Add the following processes to [concepts.nf](concepts.nf):
   ```groovy
   process download {
      input:
      tuple val(name), val(url)
      output:
      tuple val(name), path(output)
      script:
      output = file(url).name
      """
      wget $url
      """
   }
   ```
   ```groovy
   process jp2a {
      container 'talinx/jp2a' 
      input:
      tuple val(name), path(image)
      output:
      tuple val(name), stdout
      script:
      """
      printf '\n'
      jp2a $image --width=40 --color-depth=24 --fill --border
      """
   }
   ```
   `jp2a` (jpeg to ASCII) is a program that converts an image file into an ASCII (text) representation.
2. Modify the `logos` channel created in 2.3.3 by piping into the processes `download` and then `jp2a` and then `view()` the result
3. Using the joined `languages` and `logos` channels (as in 2.3.4), use the view operator with a closure to print the following text followed by the ASCII logo:
   ```
   Python was created in 1991 by Guido van Rossum:
   <ASCII logo here>
   ```


## 4. Configutation
* Look at [nextflow.config](nextflow.config)
* When a file named `nextflow.config` is present in the same directory as a nextflow script, it provides project-level configuration to be used when running that script. 
* Look at `~/.nextflow/config`. This provides system wide nextflow configuration, and is tailored to Milton/SLURM (it was created when you loaded the nextflow module).
* Any settings provided by both the system wide `~/.nextflow/config` and project `nextflow.config` are overridden by the project `nextflow.config` 
* see https://www.nextflow.io/docs/latest/config.html
