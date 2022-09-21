# Module 1: Hello World

## Learning Objectives
1. run a Nextflow script
1. assign variables
2. understand `process` definition
   * `output: stdout`
   * `script:`
3. understand `workflow` definition
5. use the `view` operator

## 1. hello_world.nf

* [hello_world.nf](hello_world.nf) is a simple Nextflow script.
   ```nextflow
   message = 'Hello world' // <-- (1) variable assignment

   // (2) process definition
   process greet {
      output: stdout
      script: "echo -n $message" // <-- (3) variable interpolation
   }

   // 4) workflow definition
   workflow {
      greet() | view  // <-- (5) piping
   }
   ```
1. We can assign variables as in most scripting languages
2. Nextflow processes define units of computations carried out by other software. The `script:` section defines a bash script to be run.
3. Variable interpolation: `$message` is replace with `Hello world`
4. Workflow definition. Workflows define connections (Channels) bewteen `processes` and `operators`

## 2. Connecting processes
* We can connect processes together using the pipe operator.
   ```nextflow
   greeting = 'Hello world'
   question = 'how are you'

   process greet {
      output: stdout
      script: "echo -n $greeting"
   }

   process ask_question {
      input: val(input)
      output: stdout
      script: "echo -n $input, $question?"
   }

   workflow {
      greet() | 
         ask_question |
         view  
   }
   ```
* **Exercise:** Add another process adding additional text to the output of `ask_question` (e.g. "Some weather we've had lately!")