# FAQ

## Q: MINI-EX failed, what should I do now?
A: 
* Check the [config file](/docs/configuration.md): did you specify the correct [executor](https://www.nextflow.io/docs/latest/executor.html) (e.g. SGE, SLURM, ...)?
* Check the Nextflow output to see which specific process failed. Try to increase the memory allocated to that process in the config file. For example, if the process "makeBorda" failed, change `-l h_vmem=10G` to `-l h_vmem=20G` for the makeBorda process in the config file.

## Q: The global Borda rank and the cluster-specific Borda rank are not following each other exactly. Is this a bug?
A:
No, this is not a bug. There is not necessarily a monotonic relationship between the global Borda rank and the Borda rank of a subset. The example below shows this in more detail for a simple example of Borda ranking.
![global-vs-subset_Borda](/docs/Difference_global-vs-subset_Borda.PNG)
