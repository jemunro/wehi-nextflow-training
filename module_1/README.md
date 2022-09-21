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
   greeting = 'Hello world' // <-- (1) variable assignment

   // (2) process definition
   process greet {
      output: stdout
      script: "echo -n $greeting" // <-- (3) variable interpolation
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
### **Exercise 1.1**
* Run `hello_world.nf`
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf
   ```

## 2. Connecting processes
* We can connect processes together using the pipe operator. Here we have added a new process `ask_question` that takes the output of `greet` and adds a question to it.
   ```nextflow
   greeting = 'Hello world'
   question = 'how are you'

   process greet {
      output: stdout
      script: "echo -n $greeting"
   }

   process ask_question {
      input: val(greet_output)
      output: stdout
      script: "echo -n $greet_output, $question?"
   }

   workflow {
      greet() | 
         ask_question |
         view  
   }
   ```

### **Exercise 1.2**
*  Add another process that appends additional text to the output of `ask_question` (e.g. "Some weather we've had lately!")

## 3. Workflow Parameters
* Nextflow workflows container a special variable `params` which is a Map (equivalent to a named list in `R` or a dict in `python`)
* We can set default values for `params` in the Nextflow script
   ```nextflow
   params.greeting = 'Hello world'
   question = 'how are you'

   process greet {
      output: stdout
      script: "echo -n $params.greeting"
   }

   process ask_question {
      input: val(greet_output)
      output: stdout
      script: "echo -n $greet_output, question?"
   }

   workflow {
      greet() | 
         ask_question |
         view  
   }
   ```
* We can override parameters with command line arguments, e.g.
   ```
   nextflow run ~/wehi-nextflow-training/module_1/hello_world.nf --greeting 'Hey'
   ```
### **Exercise 1.3**
*  Move the variable `question` to a parameter and run the workflow overriding this value

### 5. Under the Hood
