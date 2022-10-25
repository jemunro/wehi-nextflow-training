## Module 4: NGS Variant Calling Pipeline

### Learning Objectives
1. Understand workflow modularisation
1. Complete an NGS variant calling pipeline

## Modularisation
* Nextflow DSL2 pipelines can be separated into modules and included in other scripts
* This makes source code easier to maintain and share between pipelines
* Take a look at [main.nf](main.nf), we can see the processes `DOWNLOAD` and `DOWNLOAD_FASTQS` are included from the source file [modules/download.nf](modules/download.nf)

### Configuration
* Look at [nextflow.config](nextflow.config)
* This pipeline has been configured to run on the slurm queue
* The Singularity `cacheDir` is set so you can use pre-downloaded container images, otherwise Nextflow would download them automatically


## **Exercise 4**
* The pipeline in [main.nf](main.nf) is incomplete. 
* We will work through completing the "TODO" sections of each process until the pipeline is complete
* Note the use of the operator `.take(3)` on the channel `fastq_url_ch`.
    * This will take just the first 3 items from the channel
    * This is very useful to the purpose of testing and developing a pipeline

### **Exercise 4.1**
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf
    ```

### **Exercise 4.2**
1. Complete process `INDEX_REF` in [modules/index_ref.nf](modules/index_ref.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `INDEX_REF`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```

### **Exercise 4.3**
1. Complete process `BWA_MEM_ALIGN` in [modules/bwa_mem_align.nf](modules/bwa_mem_align.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BWA_MEM_ALIGN`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```

### **Exercise 4.4**
1. Complete process `SAMTOOLS_SORT` in [modules/samtools_sort.nf](modules/samtools_sort.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `SAMTOOLS_SORT`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```

### **Exercise 4.5**
1. Complete process `BCFTOOLS_CALL` in [modules/bcftools_call.nf](modules/bcftools_call.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BCFTOOLS_CALL`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```

### **Exercise 4.6**
1. Complete process `BCFTOOLS_MERGE` in [modules/bcftools_merge.nf](modules/bcftools_merge.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BCFTOOLS_MERGE`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```

## Project scripts
* Scripts placed in the bin directory may be called from a process
* take a look at [bin/plot_variants.R](bin/plot_variants.R), which is used by the process `PLOT_VARIANTS` in [modules/plot_variants.nf](modules/plot_variants.nf)

### **Exercise 4.7**
1. Complete process `PLOT_VARIANTS` in [modules/plot_variants.nf](modules/plot_variants.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `PLOT_VARIANTS`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
