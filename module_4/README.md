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
    <details>
    <summary>Solution</summary>

    ```nextflow
    process INDEX_REF {
        cpus 1
        memory '1 GB'
        time '1 h'
        module 'bwa'
        module 'samtools'

        input:
        path(ref_fasta_gz)

        output:
        tuple path('ref.fasta'), path("ref.fasta.*")

        script:
        """
        gzip -cd $ref_fasta_gz > ref.fasta
        samtools faidx ref.fasta
        bwa index ref.fasta
        """
    }
    ```
    </details>
    

### **Exercise 4.3**
1. Complete process `BWA_MEM_ALIGN` in [modules/bwa_mem_align.nf](modules/bwa_mem_align.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BWA_MEM_ALIGN`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
    <details>
    <summary>Solution</summary>

    ```nextflow
    process BWA_MEM_ALIGN {
        cpus 2
        memory '2 GB'
        time '2 h'
        module 'bwa'
        tag "$sample"
        module 'samtools'

        input:
        tuple val(sample), path(fastq1), path(fastq2)
        tuple path(ref_fasta), path(ref_indices)

        output:
        tuple val(sample), path(bam)

        script:
        bam = sample + '.bam'
        """
        bwa mem -M -t $task.cpus -R '@RG\\tID:$sample\\tSM:$sample' $ref_fasta $fastq1 $fastq2 |
            samtools view -b > $bam
        """
    }
    ```
    </details>

### **Exercise 4.4**
1. Complete process `SAMTOOLS_SORT` in [modules/samtools_sort.nf](modules/samtools_sort.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `SAMTOOLS_SORT`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
    <details>
    <summary>Solution</summary>

    ```nextflow
    process SAMTOOLS_SORT {
        cpus 2
        memory '2 GB'
        time '1 h'
        module 'samtools'
        tag "$sample"

        input:
        tuple val(sample), path(input_bam)

        output:
        tuple val(sample), path(sorted_bam), path(bam_index)

        script:
        sorted_bam = sample + '.sorted.bam'
        bam_index = sorted_bam + '.bai'
        """
        samtools sort --threads $task.cpus $input_bam > $sorted_bam
        samtools index $sorted_bam
        """
    }
    ```
    </details>

### **Exercise 4.5**
1. Complete process `BCFTOOLS_CALL` in [modules/bcftools_call.nf](modules/bcftools_call.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BCFTOOLS_CALL`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
    <details>
    <summary>Solution</summary>

    ```nextflow
    process BCFTOOLS_CALL {
        cpus 2
        memory '2 GB'
        time '1 h'
        container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        tag "$sample"

        input:
        tuple val(sample), path(sorted_bam), path(bam_index)
        tuple path(ref_fasta), path(ref_indices)

        output:
        path bcf

        script:
        bcf =  sample + '.bcf'
        """
        bcftools mpileup -Ou -f $ref_fasta $sorted_bam | bcftools call -mv -Ob -o $bcf
        """
    }
    ```
    </details>


### **Exercise 4.6**
1. Complete process `BCFTOOLS_MERGE` in [modules/bcftools_merge.nf](modules/bcftools_merge.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `BCFTOOLS_MERGE`
1. Run [main.nf](main.nf)
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
    <details>
    <summary>Solution</summary>

    ```nextflow
    process BCFTOOLS_MERGE {
        cpus 2
        memory '2 GB'
        time '1 h'
        container "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"

        input:
        path(bcfs)

        output:
        path merged_vcf
        
        script:
        merged_vcf = 'merged.vcf.gz'
        """
        bcftools merge --threads $task.cpus --no-index --missing-to-ref -Oz $bcfs > $merged_vcf
        """
    }
    ```
    </details>

## Project scripts
* Scripts placed in the bin directory may be called from a process
* take a look at [bin/plot_variants.R](bin/plot_variants.R), which is used by the process `PLOT_VARIANTS` in [modules/plot_variants.nf](modules/plot_variants.nf)

### **Exercise 4.7**
1. Complete process `PLOT_VARIANTS` in [modules/plot_variants.nf](modules/plot_variants.nf)
1. Uncomment the corresponding lines in [main.nf](main.nf) referencing `PLOT_VARIANTS`
1. Uncomment the `workflow.onComplete {...}` section in [main.nf](main.nf)
1. Run [main.nf](main.nf) and look at the output plot
    ```
    nextflow run ~/wehi-nextflow-training/module_4/main.nf -resume
    ```
    <details>
    <summary>Solution</summary>

    ```nextflow
    process PLOT_VARIANTS {
        cpus 1
        memory '2 GB'
        time '1 h'
        container 'library://jemunro/training/tidyverse-pheatmap'
        publishDir "results", mode: 'copy'

        input:
        path(vcf)
        path(metadata)

        output:
        path(plot)

        script:
        plot = 'plot.png'
        """
        plot_variants.R $vcf $metadata $plot
        """
    }
    ```
    </details>