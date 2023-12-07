# Configuring MINI-EX

The config file of MINI-EX is composed by **three main sections**:

## **Exectutor scope**
MINI-EX is configured to run on a SLURM computer cluster as specified in the executor configuration scope at the top of the config file (to use a different executor, please check [the Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html)).
MINI-EX uses a Singularity container built on Docker image with the necessary Python 3.6 modules.

```
// SPECIFY THE SYSTEM EXECUTOR //
executor {
    name = 'slurm'
    //  --> common examples: local, slurm, sge, azurebatch
    queueSize = 5
}

// SINGULARITY CONTAINERS SETTINGS // (not to be adapted)
process.container = "vibpsb/mini-ex:latest"
singularity {
    enabled = true
    cacheDir = "singularity_cache"
    autoMounts = true
}
```

**name**, **queueSize** and [optional settings](https://www.nextflow.io/docs/latest/config.html) can be changed according to resource availability. 

## **Params scope**
The params section is meant to define parameters and paths that will be used by the pipeline scripts.  
The **first block** defines the paths to the single-cell [input files](https://github.com/VIB-PSB/MINI-EX/tree/main/example/INPUTS) required and provided by the user. 

```
    // INPUT FILES DERIVED FROM SINGLE-CELL DATASETS //
    expressionMatrix = "$baseDir/example/INPUTS/*_matrix.txt"
    markersOut = "$baseDir/example/INPUTS/*_allMarkers.txt"
    cellsToClusters = "$baseDir/example/INPUTS/*_cells2clusters.txt"
    clustersToIdentities = "$baseDir/example/INPUTS/*_identities.txt"
```
The **second block** allows to specify the path to a precomputed GRNBoost2 network. By default, this is set to the GRNBoost2 network of the example run. However, when running MINI-EX for the first time on your own data, `grnboostOut` should be set to `null` in order to run GRNBoost2 on your dataset.

```
    // CACHED FILES // (allows to reuse the initial coexpression network)
    grnboostOut = "/$baseDir/example/OUTPUTS/GRNBoost2_output/*_grnboost2.txt" 
    // --> to infer coexpression network de novo, replace the line above by: grnboostOut = null
```

The **third block** is composed of paths to species-specific files provided by MINI-EX for supported species.  
There are currently four directories, one for each species supported by MINI-EX: *Arabidopsis thaliana* (data/ath), *Oryza sativa* (data/osa), *Solanum lycopersicum* (data/sly) and *Zea mays* (data/zma).

This consists in:
* Ensemble motif-mapping file obtained by selecting the top motif matches for [Cluster Buster](https://github.com/weng-lab/cluster-buster) and [FIMO](https://meme-suite.org/meme/doc/fimo.html) matches for a collection of motifs mapped to the regulatory regions of the species.
* TF-motif Family file linking each TF to its TF family and reporting wheather direct motif information is available or not, and the motifs directly associated to the TF.
* Gene-GO links for BP (biological process).
For *A. thaliana* only those supported by experimental evidence ("EXP", "IMP", "IDA", "IPI", "IGI", "IEP") and manually curated ("TAS", "NAS", "IC") are included, while for other species all the evidences are kept.
Note: all ancestral terms are included and terms associated with more than 30% of all the genes of the species are filtered out.
* Gene-alias file downloaded from [TAIR](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FPublic_Data_Releases%2FTAIR_Data_20140331), [funRiceGenes](https://funricegenes.github.io/), [Sol Genomics Network](https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_descriptions.txt) and [MaizeGDB](https://www.maizegdb.org/associated_genes?type=all&style=tab).

```	
    // SPECIES SPECIFIC INFORMATION //
    tfList = "$baseDir/data/ath/ath_TF_list.txt"
    geneAliases = "$baseDir/data/ath/ath_gene_aliases.txt"
    infoTf = "$baseDir/data/ath/ath_TF2fam2mot.txt"
    featureFileMotifs = "$baseDir/data/ath/ath_2021.1_motifMapping.out.gz"
    goFile = "$baseDir/data/ath/ath_full_BP_expcur_ext_names.txt" 
    // --> if GO data is not available, replace the line above by: goFile = null 
    //     (when doing so, termsOfInterest should also be set to <null>)
```
The **fourth block** consists of MINI-EX parameters that can be adapted to change the behavior of MINI-EX.

When the first parameter `doMotifAnalysis` is set to `false`, the motif mapping step in the MINI-EX workflow (step 2) is omitted. This yields less precise GRNs, but allows to run MINI-EX without species-specific information which are only shipped with MINI-EX for supported species.

The `termsOfInterest` parameter is the path to a user-defined file with GO terms of interest. This file should contain terms (e.g. root, phloem, xylem, etc.) related to the process or tissue of interest and will affect the ranking of the final regulons. When set to `null`, no GO information is taken into account for ranking the regulons and MINI-EX will use its standard ranking procedure.

The `topMarkers` parameter refers to the number of upregulated genes per cluster (sorted by q-value) to use during the cell cluster enrichment. `expressionFilter` refers to the percentage of cells that need to express the TF to retain the regulon for the specific cell cluster. Default values for both parameters have been chosen by benchmarking different filters against a root gold standard of known protein-DNA interactions.

* `motifFilter` can be set to either **TF-F_motifs** (default) or **TF_motifs**. The default option keeps regulons in the motif filtering step (step 2) if the regulon is enriched for any motif of that TF family. When setting the parameter to **TF_motifs**, then only regulons are retained if they are enriched for direct TF-motifs.

* `topRegs` defines the top regulons to show in the two output heatmaps. It can be changed according to the user needs.

```	
    // PARAMETERS //
    doMotifAnalysis = true 
    // --> set to <false> if no motif mapping data is available 
    //     [CAUTION: without motif data MINI-EX is less reliable]
    termsOfInterest = "$baseDir/example/INPUTS/GOsIwant.txt"
    // --> to use the standard ranking procedure, replace the line above by: termsOfInterest = null
    topMarkers = "700"
    expressionFilter = "10"
    motifFilter = "TF-F_motifs" 
    // --> to use the motifs of the TF family: motifFilter = "TF-F_motifs"
    // --> to only use the motifs known for a TF: motifFilter =  "TF_motifs"
    topRegulons = "150"
```

The **last block** specifies the output directory.  

```	
    // SPECIFY OUTPUT DIRECTORY //
    outputDir = "$baseDir/output"
}
```

## **Process scope**
The process scope allows the user to define the configuration (memory, parallel environment) for each of the processes run by MINI-EX.  
These can be changed according to the resources availabe/needed. 
   
```	
process {
    ...
    withName: run_grnboost {
        memory = '20 GB'
        cpus = 5
    }
    ...
}
```
