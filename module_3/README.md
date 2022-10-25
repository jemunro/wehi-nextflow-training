# Module 3: Nextflow Scripting

### Learning Objectives
1. Groovy scripting basics
2. Understand channel creation
3. Use operators to transform channels

## 3.1 Groovy Scripting
* Much of Nextflow scripting is done using the Groovy language
* **Numeric Operations**
   ```groovy
   x = 2                      // integer variable
   y = 2.5                    // float variable 
   z = 0.99                   // float variable 
   assert x * y == 5          // assertion, throw error if not true
   assert Math.round(z) == 1  // round to nearest integer
   assert [x, y].max() == 2.5 // maximum 
   ```
* **String Operations**
   ```groovy
   s1 = 'foo'                 // String variable
   s2 = 'bar'                 // String variable
   println(s1)                // print value to standard output
   c1 = s1 + '-' + s2         // String concatenation
   c2 = "$s1-$s2"             // String interpolation of variables
   assert c1 == c2  
   s3 = "1 + 2 = ${1 + 2}"    // Strining interpolation of closure
   assert s3 == "1 + 2 = 3"
   ```

* **Logic Operations**
   ```groovy
   isValid = true // boolean variable
   // ------------ if else statement ------------ //
   if (isValid) {           
      println("valid")
   } else {
      println("not valid")
   }
   // ------------ ternary operator ------------ //
   isValid ? println("valid") : println("not valid")
   ```
* See example code snippits in Groovy and Python: https://programming-idioms.org/cheatsheet/Groovy/Python

## 3.2 Lists & Maps
* **Lists**
   * Nextflow/Groovy:
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
### **Exercise 3.2**
1. Try out some of above Groovy examples using the groovy shell:
   ```
   module load java/1.8.0_92 groovy/4.0.0
   groovysh
   ```

## 3.3 Closures
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
* If we wanted to add to the inner list, we could instead do this:
   ```groovy
   nested.collect { number, letter -> [number, letter, "$number-$letter"] }
   ```
   ```
   [[1, A, 1-A], [2, B, 2-B], [3, C, 3-C]]
   ```
* See https://www.nextflow.io/docs/latest/script.html#closures

### **Exercise 3.3**
1. In the Groovy shell, define the variable `data` as below
   ```groovy
   data = [['foo', 1, 2], ['bar', 3, 4], ['baz', 5, 6]]
   ```
   Using `collect`, transform the data to the following list:  
   `["foo: 3", "bar: 12", "baz: 30"]`  
   (the number is equal to the product of the numbers in the list)
   <details>
   <summary>Solution</summary>

   ```groovy
   data.collect { s, n1, n2 -> "$s: ${n1*n2}" } 
   ```
   </details>
1. In the Groovy shell, define the variable `data` as below
   ```groovy
   data = [['foo', 1, 2], ['bar', 3, 4], ['baz', 5, 6]]
   ```
   Using `collect`, transform data to the following nested list:  
   `[['foo', 1, 2, 2], ['bar', 3, 4, 12], ['baz', 5, 6, 30]]`
   <details>
   <summary>Solution</summary>

   ```groovy
   data.collect { s, n1, n2 -> [s, n1, n2, n1 * n2] }
   ```
   </details>


## 3.4 Nextflow Scripting
**Implicit Variables**:
* A number of variables are available in all nextflow scripts:
* `params`: map storing workflow parameters
* `projectDir`: A string variable of the directory containing the nextflow script being run. Useful for accessing workflow resource files.
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

## 3.5 Channels
* Nextflow includes a number of ways to create channels
* `channel.of(...)` 
  * create a channel that emits each of the arguments one at a time.
      ```groovy
      channel.of('A', 'B', 'C')
      ```
* `channel.fromList(list)` 
   * given a list, create a channel that emits each element of the list one at a time
      ```groovy
      list = ['A', 'B', 'C']
      channel.fromList(list)
      ```
* `channel.fromPath(path)`: 
  * Create a channel from a file path, emitting a `file` variable. Functions similarly to the `file()` method but returns a channel.
      ```groovy
      channel.fromPath('/path/to/sample.bam')
      ```
  * Optional argument `checkIfExists` will throw an error if file does not exists
      ```groovy
      channel.fromPath('/path/to/sample.bam', checkIfExists: true)
/home/users/allstaff/munro.j/pipelines/wehi-nextflow-training/examples      ```
   * Files may be located on the web and will downloaded when needed:
       ```groovy
      channel.fromPath('https://www.wehi.edu.au/sites/default/files/wehi-logo-2020.png')
      ```
   * See https://www.nextflow.io/docs/latest/channel.html#channel-factory for more ways to create channels
 ### **Exercise 3.5**
1. Open [languages.nf](languages.nf) and look at the CSV files in [data](data)
1. Create and `view()` channels `authors_ch` and `homepage_ch` in the same way as `year_created_ch`
1. Run [languages.nf](languages.nf)
   ```
   nextflow run ~/wehi-nextflow-training/module_3/languages.nf
   ```
   <details>
   <summary>Solution</summary>

   ```nextflow
   workflow {

      year_created_ch = Channel.fromPath("$projectDir/data/year_created.csv", checkIfExists: true)
      year_created_ch.view()

      authors_ch = Channel.fromPath("$projectDir/data/authors.csv", checkIfExists: true)
      authors_ch.view()

      homepage_ch = Channel.fromPath("$projectDir/data/homepage.csv", checkIfExists: true)
      homepage_ch.view()
   }
   ```
   </details>

## 3.6 Operators
### Map
* `map` is the most commonly used nextflow operater, and works the same as Groovy's `collect()` but applied to channels instead of lists.
* Functionally similar to R's `lapply()` or Python's `map()` 
* The default variable name for a closure is `it`
* for example:
   ```groovy
   channel.of(1, 2, 3).map { it * 2 }.view()
   ```
   will print:
   ```
   2
   4
   6
   ```
### View
* `view` takes an input channel and prints the contents to the terminal
* We can optionally provide a closure, similarly to `map`, which will instead print the value of applying the closure
   ```groovy
   channel.of(1, 2, 3).view { x -> "$x * 2 = ${x * 2}" } 
   ```
   will print:
   ```
   1 * 2 = 2
   2 * 2 = 4
   3 * 2 = 6
   ```
### splitCsv
* Convert a CSV (or TSV) file into a nextflow channel
* Optional argument `skip` may be used to skip a number of lines (e.g. the header)
* Most often used on the output of `channel.fromPath()`, e.g.
   ```groovy
   channel.fromPath("$projectDir/data/year_created.csv")
      .splitCsv(skip: 1)
      .view()
   ```
   will print:
   ```
   [Python, 1991]
   [R, 1993]
   [Nextflow, 2013]
   ```
 ### **Exercise 3.6.1**
1. Update all channels in [languages.nf](languages.nf) with splitCsv as follows:
   ```nextflow
      workflow {
         year_created_ch = Channel
            .fromPath("$projectDir/data/year_created.csv", checkIfExists: true)
            .splitCsv(skip: 1)
         year_created_ch.view()
      }
   ```
1. Run [languages.nf](languages.nf)
   <details>
   <summary>Solution</summary>

   ```nextflow
   workflow {
      year_created_ch = Channel
         .fromPath("$projectDir/data/year_created.csv", checkIfExists: true)
         .splitCsv(skip: 1)
      year_created_ch.view()

      authors_ch = Channel
         .fromPath("$projectDir/data/authors.csv", checkIfExists: true)
         .splitCsv(skip: 1)
      authors_ch.view()

      homepage_ch = Channel
         .fromPath("$projectDir/data/homepage.csv", checkIfExists: true)
         .splitCsv(skip: 1)
      homepage_ch.view()
   }
   ```
   </details>

### join
* Join two channels by a matching key
   ```groovy
   left  = channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
   right = channel.of(['Z', 6], ['Y', 5], ['X', 4])
   left.join(right).view()
   ```
   will print:
   ```
   [Z, 3, 6]
   [Y, 2, 5]
   [X, 1, 4]
   ```
* This operation is roughly equivalent to R's `dplyr::inner_join()` and Python's pandas `pd.merge()`. For example, in R:
   ```R
   left  = data.frame(V1 = c('X', 'Y', 'Z', 'P'),
                      V2 = c( 1,   2,   3,   7))
   right = data.frame(V1 = c('Z', 'Y', 'X'),
                      V3 = c( 6,   5,   4))
   left %>% inner_join(right) %>% View()
   ```
   ```
     V1 V2 V3
   1  X  1  4
   2  Y  2  5
   3  Z  3  6
   ```

* Many more operators are availabe, see https://www.nextflow.io/docs/latest/operator.html

### **Exercise 3.6**
1. Remove calls to the `view()` operator for channeles `year_created_ch`, `authors_ch` and `homepage_ch`
1. Using the `join()` operator, created a new channel `joined_ch` that joins `year_created_ch`, `authors_ch` and `homepage_ch` and `view()` the output
   <details>
   <summary>Solution</summary>

   ```nextflow
   workflow {

      year_created_ch = Channel
         .fromPath("$projectDir/data/year_created.csv", checkIfExists: true)
         .splitCsv(skip: 1)

      authors_ch = Channel
         .fromPath("$projectDir/data/authors.csv", checkIfExists: true)
         .splitCsv(skip: 1)

      homepage_ch = Channel
         .fromPath("$projectDir/data/homepage.csv", checkIfExists: true)
         .splitCsv(skip: 1)

      joined_ch = year_created_ch
         .join(authors_ch)
         .join(homepage_ch)
         .view()
   }
   ```
   </details>

1. Using the `view {}` operator on `joined_ch`, print the following:
   ```
   R was created in 1993 by Ross Ihaka and Robert Gentleman. To learn more vist https://www.r-project.org/
   Python ...
   Nextflow ...
   ```
   <details>
   <summary>Solution</summary>

   ```nextflow
      joined_ch = year_created_ch
         .join(authors_ch)
         .join(homepage_ch)
         .view { lang, year, auth, url -> 
            "$lang was created in $year by $auth. To learn more vist $url" }
   ```
   </details>


