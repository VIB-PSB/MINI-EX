# FAQ

## Q: MINI-EX failed, what should I do now?
A: 
* Check the [config file](/docs/configuration.md): did you specify the correct [executor](https://www.nextflow.io/docs/latest/executor.html) (e.g. SGE, SLURM, ...)?
* Check the Nextflow output to see which specific process failed. Try to increase the memory allocated to that process in the config file. For example, if the process `make_borda` failed, change `memory = '10 GB'` to `memory = '20 GB'` for that process in the config file.
* If the process failed after approximately 1 hour, then you might be running MINI-EX on a cluster that applies a default wall time of 1 hour. This means that processes are automatically killed if they did not finish after an hour, unless you explicitly define a longer wall time in the config file. This can be done by adding the following line to the process in the config file: `time = '48 h'`. See [this issue](https://github.com/VIB-PSB/MINI-EX/issues/24) for more details.

## Q: The global Borda rank and the cluster-specific Borda rank are not following each other exactly. Is this a bug?
A:
No, this is not a bug. There is not necessarily a monotonic relationship between the global Borda rank and the Borda rank of a subset. The example below shows this in more detail for a simple example of Borda ranking.
![global-vs-subset_Borda](/docs/Difference_global-vs-subset_Borda.PNG)
