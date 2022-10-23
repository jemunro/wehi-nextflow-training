# Module 3: Processes and Configuration

### Learning Objectives
1. TBD

## 1. Configutation
* Look at [nextflow.config](nextflow.config)
* When a file named `nextflow.config` is present in the same directory as a nextflow script, it provides project-level configuration to be used when running that script. 
* Look at `~/.nextflow/config`. This provides system wide nextflow configuration, and is tailored to Milton/SLURM (it was created when you loaded the nextflow module).
* Any settings provided by both the system wide `~/.nextflow/config` and project `nextflow.config` are overridden by the project `nextflow.config` 
* see https://www.nextflow.io/docs/latest/config.html

## 2. Scripts (bin) (TBD)

## 2.2 Processes
### **Inputs**
* `val()` - A val type input denotes a regular groovy variable. It could be a String, Integer, Boolean, double etc.
* `path()` - A path represents an input file.
* `tuple` - A tuple represents a list of inputs. There may be of either `val` or `path` types
### **Outputs**
* similarly to inputs, outputs may be of type `val`, `path` or `tuple`
* `path()` - when a `path` is declared as an output, after the process has run sucessfully nextflow will check that the path exists, and throw an error if not.
* `stdout` - stdout is a special output type that will return the standard output of the process run