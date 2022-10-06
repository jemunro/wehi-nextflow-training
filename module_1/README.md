# Module 1: Hello World

## Learning Objectives
1. Run and modify a Nextflow script
1. Understand `process` definition
1. Understand `workflow` definition
1. Understand piping between channels and processes
1. Understand and use pipeline `params`
1. Understand pipeline results caching

## 1. hello_world.nf

* [hello_world.nf](hello_world.nf) is a simple Nextflow script.
   ```nextflow
   audience = ['world']            // (1) variable assignment

   // -------- (2) process definition --------
   process greet {
      input: val(x) 
      output: stdout
      script: "echo -n Hello $x!"  // (3) variable interpolation
   }

   // -------- (4) workflow definition --------
   workflow {
      channel.fromList(audience) | // (5) Channel creation
        greet | 
        view  
   }
   ```
1. We can assign variables as usual
2. Nextflow processes define units of computations carried out by other software. The `script:` section defines a bash script to be run.
3. Variable interpolation: `$x` is replaced with the value ("world")
4. Workflows define connections (Channels) bewteen `processes` and `operators`.
### **Exercise 1.1**
1. Run `hello_world.nf`
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf
   ```
2. Add additional items to the list `audience` to be a list of values (e.g. `audience = ['world',  'WEHI']`) and run `hello_world.nf`? What happens?


## 2. Connecting processes
* We can connect processes together using the pipe operator. Here we have added a new process `ask_question` that takes the output of `greet` and adds a question to it.

   ```nextflow
   question = 'how are you?'

   process ask_question {
      input: val(x)
      output: stdout
      script: "echo -n $x, $question"
   }
   ```
### **Exercise 1.2**
1. Add the above process definition to `hello_world.nf`
2. Add `ask_question |` to the workflow definition between `greet |` and `view`
3. Run `hello_world.nf`
   <details>
   <summary>Solution</summary>

   ```nextflow
   audience = ['world', 'WEHI']

   process greet {
      input: val(x)
      output: stdout
      script: "echo -n Hello $x!"
   }

   question = 'how are you?'

   process ask_question {
      input: val(x)
      output: stdout
      script: "echo -n $x, $question"
   }

   workflow {
      channel.fromList(audience) |
         greet |
         ask_question |
         view
   }
   ```
   </details>

## 3. Workflow Parameters
* Nextflow workflows contain a special variable `params` which is a 'map' (equivalent to a named list in `R` or a dict in `python`)
* We can set default values for `params` in the Nextflow script, e.g.:
   ```nextflow
   params.greeting = 'Hello'

   process greet {
      input: val(x)
      output: stdout
      script: "echo -n $params.greeting $x!"
   }
   ```
* We can then override parameters with command line arguments, e.g.:
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf --greeting 'Hey'
   ```
* See https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters
### **Exercise 1.3**
1. Convert `greeting` to an input paramter in `hello_world.nf` as in the example above. Run `hello world.nf`, providing `--greeting` as a command line arguemnt.
2. Also convert `question` to an input paramter. Run `hello world.nf`, providing both `--greeting` and `--question` as a command line arguments.
   <details>
   <summary>Solution</summary>

   ```nextflow
   params.greeting = 'Hello'

   process greet {
      input: val(x)
      output: stdout
      script: "echo -n $params.greeting $x!"
   }

   params.question = 'how are you?'

   process ask_question {
      input: val(x)
      output: stdout
      script: "echo -n $x, $params.question"
   }

   workflow {
      channel.fromList(audience) |
         greet |
         ask_question |
         view
   }
   ```
   </details>

