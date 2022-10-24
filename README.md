# WEHI Nextflow Training

## Prerequisites
* To complete the training, you will require an appropriate IDE/editor installed on your computer.
* [Visual Studio Code](https://code.visualstudio.com/) is recommended, but for other options see https://nf-co.re/developers/editor_plugins. 
* You will be editing code stored in your unix home directory through a WEHI NAS connection, see [how to connect to NAS](https://wehieduau.sharepoint.com/sites/rc2/Shared%20Documents/Forms/AllItems.aspx?id=%2Fsites%2Frc2%2FShared%20Documents%2FUsing%20Milton%2Fhow%20to%20connect%20to%20NAS%5FPC%2DMac%5F%20Mar2022%2Epdf&parent=%2Fsites%2Frc2%2FShared%20Documents%2FUsing%20Milton).

If using Visual Studio Code, the following extensions should be installed:
* https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow
* https://marketplace.visualstudio.com/items?itemName=ms-vsliveshare.vsliveshare
* https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh
    * Using this extension is a good alternative to using NAS to edit remote files.
You can connect using vc7-shared as the remote host, for more details see [remote coding with Visual Studio](https://wehieduau.sharepoint.com/sites/RCPNewsletter/SitePages/June-2022-Issue.aspx#remote-coding-with-visual-studio-code)

## Setup (on the day)
1. Open a termnial and SSH into vc7-shared
1. Clone this repostiory to your home directory
    ```
    cd ~
    git clone https://github.com/jemunro/wehi-nextflow-training.git
    ```
2. Create and navigate to a run directory on your vast scractch space
    ```
    mkdir -p /vast/scratch/users/$USER/nf-training-run
    cd /vast/scratch/users/$USER/nf-training-run
    ```
3. (Optional) start an interactive slurm session
    ```
    salloc --partition interactive --mem 8G --cpus-per-task 4 --time 2:00:00
    ```
4. Load the nextflow module
    ```
    module load nextflow/22.04.5
    ```
5. Open the VSCode liveshare link in your browser


## [Module 1: Hello World](module_1/README.md)
A brief introduction to the Nextflow language.

## [Module 2: Processes and Software Dependencies](module_2/README.md)
Configuration of processes and software dependencies using an example RNA-seq pipeline.

## [Module 3: Nextflow Scripting](module_3/README.md)
Scripting with Groovy and Nextflow. Creating channels and using operators.

## [Module 4: NGS Variant Calling Pipeline](module_4/README.md)
Work though an example pipeline calling SARS-Cov-2 variants using bwa, samtools, bcftools and then visualising with R.

## Using Nextflow on Milton
* Use `screen` or `tmux` to leave Nextflow pipelines running on vc7-shared
* Alternatively, submit pipelines to the long slurm partition
* Be aware of the VAST scratch deletion policy, make sure to transfer pipeline outputs to permanent storage

## Additional Resources
### Documentation
* [Nextflow documentaion](https://www.nextflow.io/docs/latest/basic.html)
### Training
* [Seqera Nextflow Training](https://training.seqera.io/)
### Community Nextflow Pipelines
* [nf-core](https://nf-co.re/) 
### Getting Help
* [Nextflow GitHub Discussions](https://github.com/nextflow-io/nextflow/discussions)
* [Nextflow Slack](https://www.nextflow.io/slack-invite.html)
* [nf-core Slack](https://nf-co.re/join/slack)