## Module 2: Processes and Software Dependencies
* Here we will work though an exmple RNA-seq quantification using [Salmon](https://combine-lab.github.io/salmon/)  
* This is module is adapted from https://training.seqera.io/#_simple_rna_seq_pipeline

### Learning Objectives
1. Work on a practical example of a nextflow pipeline
1. Understand process inputs and outputs
1. Understand process directives
1. Use directives `module` and `container` to specify process dependencies
1. Understand nextflow configuration

## 2.1 Pipeline Processes
* Open [rna_seq_1.nf]('rna_seq_1.nf') and take a look.

*  Here we have a single process `INDEX` defined. This process builds the Salmon Index of the transcriptome provided.
    ```nextflow
    process INDEX {
        input:
        path transcriptome

        output:
        path 'salmon_index'

        script:
        """
        salmon index -t $transcriptome -i salmon_index
        """
    }
    ```
* We can see this process definies an input of type `path` (i.e. a file)
* The output is also of type `path`, and declares that a file with name 'salmon_index' should be created. 
* If the output file does not exist after the process has run, Nextflow will throw an error.

### **Exercise 2.1**
1. Run [rna_seq_1.nf]('rna_seq_1.nf')
   ```
   nextflow run ~/wehi-nextflow-training/module_2/rna_seq_1.nf
   ```
   This will fail with error message: `line 2: salmon: command not found`. This is because we haven't provided a specification for the software required.
2. Check for salmon module on milton
    ```
    module avail salmon
    ```
    We can see that salmon is indeed available.
3. Add the directive `module salmon/1.9.0` to the `INDEX` process as follows and run the pipeline again.
    ```nextflow
    process INDEX {
        module 'salmon/1.9.0'

        input:
        path transcriptome

        output:
        path 'salmon_index'

        script:
        """
        salmon index -t $transcriptome -i salmon_index
        """
    }
    ```
    With this, Nextflow will load the appropriate module prior to running the process script.

## 2.2 Process Directives
* Directives specify the execution environment of a nextflow process. For example the `module` directive above specifies the software modules to be used.
* Directives are placed at the top of a process definition.
* See https://www.nextflow.io/docs/latest/process.html#directives for all available directives

### **Exercise 2.2**
1. Add the following directives to `INDEX` to specify the cpu and memory resources required by the process.
    ```nextflow
    memory '2 GB'
    cpus 1
    ```
2. Run the workflow again and confirm it runs successfully.
    <details>
    <summary>Solution</summary>

    ```nextflow
    process INDEX {
        module 'salmon/1.9.0'
        memory '2 GB'
        cpus 1

        input:
        path transcriptome

        output:
        path 'salmon_index'

        script:
        """
        salmon index -t $transcriptome -i salmon_index
        """
    }
    ```
    </details>


## 2.3 Procecss inputs and outputs
* Open [rna_seq_2.nf]('rna_seq_2.nf') and take a look.
* Here we have added a process `QUANTIFICATION`. This process takes the RNA-seq data and counts the reads originating from each transcrpit in the transcriptome:
    ```nextflow
    process QUANTIFICATION {
        module 'salmon/1.9.0'
        memory '2 GB'
        cpus 2
        tag "$sample_id"

        input:
        path salmon_index
        tuple val(sample_id), path(reads)

        output:
        path output

        script:
        output = "${sample_id}.sf"
        """
        salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o out
        mv out/quant.sf $output
        """
    }
    ```
* Process inputs may be of the following types:
    * `val` - A val type denotes a regular groovy variable. It could be a String, Integer, Boolean, double etc.
    * `path` - A path represents an input/output file.
    * `stdout` - stdout is a special output type that will return the standard output of the process run
    * `tuple` - A tuple represents a collection of inputs/out. These may be of either `val`, `path` or `stdout` types
* Note that `QUANTIFICATION` defines two inputs channels, one for the salmon index and one for the reads for each sample
* `Channel.fromFilePairs()` is a special method designed to handle paired file inputs from NGS sequencing

## 2.4 Software Containers
* One of the most powerful features of Nextflow is it's support for software containers (Docker, Singularity, etc.).
* Using containers will improve the reproducibility and portability of your pipelines.
* Containers can be specified using the `container` directive.
* You can find pre-made containers for popular bioinformatics software through [Bioconda](https://bioconda.github.io/)

### **Exercise 2.4**
1. Visit https://bioconda.github.io/recipes/salmon/README.html. Here we see Salmon is available in a Docker container at "quay.io/biocontainers/salmon:<tag>". If we visit the "salmon/tags" link we can find that the latest available tag is "1.9.0--h7e5ed60_1"
1. Replace the directive `module 'salmon/1.9.0'` with `container 'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1'` in the processes `INDEX` and `QUANTIFICATION` in [rna_seq_2.nf]('rna_seq_2.nf')
1. Run [rna_seq_2.nf]('rna_seq_2.nf')
    ```
    nextflow run ~/wehi-nextflow-training/module_2/rna_seq_2.nf
    ```

## 2.5 Configurtaion
* Look at `~/.nextflow/config`. This provides system wide nextflow configuration, and is tailored to Milton (it was created when you first loaded the nextflow module).
    ```nextflow
    process {
        executor = 'slurm'
        cache = 'lenient'
    }

    executor {
        name = 'slurm'
        queueSize = 100
        queueStatInterval = '10 sec'
        pollInterval = '10 sec'
        submitRateLimit = '10sec'
    }

    singularity {
        enabled = true
        autoMounts = true
        runOptions = '-B /vast -B /stornext -B /wehisan'
    }

    docker.enabled = false
    ```
* When a file named `nextflow.config` is present in the same directory as a nextflow script, it provides project-level configuration to be used when running that script. 
* Any settings provided by both the system wide `~/.nextflow/config` and project `nextflow.config` are overridden by the project `nextflow.config`
* Open at [nextflow.config](nextflow.config) and take a configuration at the settings. The default 'slurm' executor is overwritten to use the 'local' executor.
* see https://www.nextflow.io/docs/latest/config.html


## 2.6 Publishing Outputs
* Open [rna_seq_3.nf]('rna_seq_3.nf') and take a look.
* Here we have added a process `PLOT_TPM`. This process is an R script that takes RNA seq quantification results and creates a plot:
    ```nextflow
    process PLOT_TPM {
        container 'rocker/tidyverse:4.1.3'
        publishDir "results", mode: 'copy'

        input:
        path quant_results

        output:
        path 'TPM.png'

        script:
        """
        #!/usr/bin/env Rscript

        library(tidyverse)

        data.frame(filename = list.files(pattern='.sf')) %>% 
            mutate(tissue = str_remove(basename(filename), '.sf')) %>% 
            mutate(data = map(filename, read_tsv, col_types = cols())) %>% 
            unnest(data) %>% 
            select(tissue, transcript = Name, TPM) %>% 
            ggplot(aes(transcript, TPM, fill = tissue)) +
            geom_col(position = 'dodge') +
            coord_flip()

        ggsave('TPM.png', width = 6, height = 4)
        """
    }
    ```
* The `publishDir` directives specifies that output files from this process should be copied to the folder 'results'
* The operator `collect()` is used to combine all the outputs in `quant_ch` into a single input for `PLOT_TPM`
    ```nextflow
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    
    PLOT_TPM(quant_ch.collect())
    ```

### **Exercise 2.6**
1. Open at [nextflow.config](nextflow.config) and change `process.executor` from 'local' to 'slurm'. This will direct jobs to be submitted to the slurm queue.
1. Run [rna_seq_3.nf]('rna_seq_3.nf') and observe the output
    ```
    nextflow run ~/wehi-nextflow-training/module_2/rna_seq_3.nf -resume
    ```

## 2.7 Execution Log
* Nextflow executions logs can be generated after a workflow has run, see https://www.nextflow.io/docs/latest/tracing.html#execution-log

### **Exercise 2.7**
1. Run `nextflow log` to list all previous executions, and note the `RUN NAME` of the most recent execution
2. Run `nextflow log <RUN NAME> -f name,workdir,native_id,status,exit` replacing `<RUN NAME>` with the name from 1. This will list all jobs run in the previous execution.