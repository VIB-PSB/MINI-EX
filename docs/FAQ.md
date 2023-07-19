# FAQ

## Q: MINI-EX failed, what should I do now?
A: 
* Check the [config file](docs/configuration.md). Did you already substitute `<pe_name>` with the parallel environment name of your system (e.g. `serial`)? If needed, `-pe <pe_name> 5` can also be omitted. This will run the process "run_grnboost" without parallellization. Note that this will significantly increase the overall runtime.
* Check the Nextflow output to see which specific process failed. Try to increase the memory allocated to that process in the config file. For example, if the process "makeBorda" failed, change `-l h_vmem=10G` to `-l h_vmem=20G` for the makeBorda process in the config file.
