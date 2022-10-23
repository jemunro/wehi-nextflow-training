# Module 1: Hello World

## Learning Objectives
1. Run and modify a Nextflow script
1. Recognise `process` definition
1. Recognise `workflow` definition
1. Understand and use pipeline `params`
1. Understand and use workflow caching (`-resume`)

## 1.1 hello_world.nf

* [hello_world.nf](hello_world.nf) is a simple Nextflow script.
   ```nextflow
   // -------- (1) process definition -------- //
   process Greet {
      input: val(x) 
      output: stdout
      script: "echo -n Hello $x"  // (2) String interpolation
   }

   // -------- (3) workflow definition -------- //
   workflow {
      input_ch = Channel.of('world') // create input channel
      Greet(input_ch)  // run process create on input channel
      Greet.out.view() // view the output channel of process Greet
   }
   ```
1. Nextflow processes define units of computations carried out by other software. The `script:` section defines a bash script to be run.  
2. String interpolation: `$x` is replaced with the value of the variable `x`
3. Workflows define connections bewteen `channels`, `processes` and `operators`.

### **Exercise 1.1**
1. Run `hello_world.nf`
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf
   ```
2. Add additional items to the input channel (e.g. `input_ch = Channel.of('world', 'WEHI')`) and run `hello_world.nf`. This will create additional inputs for the processes to be run on.

## 1.2 Workflow Parameters
* Nextflow workflows contain a special variable `params` which is a 'map' (equivalent to a named list in R or a dict in python)
* We can set default values for `params` in the Nextflow script, e.g.:
   ```nextflow
   params.greeting = 'Hello'

   process Greet {
      input: val(x)
      output: stdout
      script: "echo -n $params.greeting $x"
   }
   ```
* We can then override parameters with command line arguments, e.g.:
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf --greeting 'Hey'
   ```
* See https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters
### **Exercise 1.2**
1. Convert `greeting` to an input paramter in `hello_world.nf` as in the example above. 
2.  Run `hello_world.nf`, providing `--greeting` as a command line arguemnt.

## 1.3 Adding Processes
Here will add another process to run on the output of `Greet`
### **Exercise 1.3**
1. Add the parameter and process definition below to `hello_world.nf`
   ```nextflow
   params.question = 'how are you?'

   process AskQuestion {
      input: val(x)
      output: stdout
      script: "echo -n $x, $question"
   }
   ```
1. Add call to process `AskQuestion` on the output of Greet, e.g. `AskQuestion(Greet.out)`
1. Remove `Greet.out.view()` and add `AskQuestion.out.view()` to the end of the workflow definition to view the output of `AskQuestion`
1. Run `hello_world.nf` with the parameter `--question`, e.g.:
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf --question 'what time is it?'
   ```

   <details>
   <summary>Solution</summary>

   ```nextflow
   params.greeting = 'Hello'
   params.question = 'how are you?'

   process Greet {
      input: val(x)
      output: stdout
      script: "echo -n $params.greeting $x"
   }

   process AskQuestion {
      input: val(x)
      output: stdout
      script: "echo -n $x, $params.question"
   }

   workflow {
      input_ch = Channel.of('world', 'WEHI')
      Greet(input_ch)
      AskQuestion(Greet.out)
      AskQuestion.out.view()
   }
   ```
   </details>

## 1.4 Workflow Caching
* Nextflow provides a mechanism to reuse results from previously run workflows, potentially saving costly processes from being recomputed.
* This works by creating a unique hash from all process inputs, and only running the process if a matching hash is not present in the work directory. If a match if found the existing outputs are used.
* To use this feature, provide the `-resume` argument at the command line:
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf -resume
   ```
* In general, we should always use this option

### **Exercise 1.4**
1. Experiment by running `hello_world.nf` with and without `-resume`, and with different parameters `--greeting` and `--question`.
