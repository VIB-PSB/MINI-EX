# FAQ

## Q: MINI-EX failed, what should I do now?
A: 
* Check the [config file](/docs/configuration.md): did you specify the correct [executor](https://www.nextflow.io/docs/latest/executor.html) (e.g. SGE, SLURM, ...)?
* Check the Nextflow output to see which specific process failed. Try to increase the memory allocated to that process in the config file. For example, if the process "makeBorda" failed, change `-l h_vmem=10G` to `-l h_vmem=20G` for the makeBorda process in the config file.
