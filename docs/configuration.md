# Configuring MINI-EX

The config file of MINI-EX is composed by **three main sections**:

## **Exectutor scope**
MINI-EX is configured to run on a SLURM computer cluster as specified in the executor configuration scope at the top of the config file (to use a different executor, please check [the Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html)).
MINI-EX uses a Singularity container built on Docker image with the necessary Python 3.6 modules.

```
executor {
    name = 'slurm'
    queueSize = 5
}

process.container = "vibpsb/mini-ex:latest"
singularity.enabled = true
```

**name**, **queueSize** and [optional settings](https://www.nextflow.io/docs/latest/config.html) can be changed according to resource availability. 

## **Params scope**
The params section is meant to define parameters and paths that will be used by the pipeline scripts.  
The **first block** of parameters defines the paths to the [input files](https://github.com/VIB-PSB/MINI-EX/tree/main/example/INPUTS) required and provided by the user.

```
params {
    expressionMatrix = "$baseDir/example/INPUTS/*_matrix.txt"
    markersOut = "$baseDir/example/INPUTS/*_allMarkers.txt"
    cell2clusters = "$baseDir/example/INPUTS/*_cells2clusters.txt"
    cluster2ident = "$baseDir/example/INPUTS/*_identities.txt"
    TF_list = "$baseDir/example/INPUTS/TF_list.txt"
    termsOfInterest = "$baseDir/example/INPUTS/GOsIwant.txt"
//  termsOfInterest = null

//  grnboostOut = "/$baseDir/example/OUTPUTS/GRNBoost2_output/*_grnboost2.txt"
    grnboostOut = null
```

The **second block** of parameters is composed by files provided and necessary for the pipeline to run.  
There are currently four directories, one for each species supported by MINI-EX: *Arabidopsis thaliana* (data/ath), *Oryza sativa* (data/osa), *Solanum lycopersicum* (data/sly) and *Zea mays* (data/zma).

This consists in:
* Ensemble motif-mapping file obtained by selecting the top motif matches for [Cluster Buster](https://github.com/weng-lab/cluster-buster) and [FIMO](https://meme-suite.org/meme/doc/fimo.html) matches for a collection of motifs mapped to the regulatory regions of the species.
* TF-motif Family file linking each TF to its TF family and reporting wheather direct motif information is available or not, and the motifs directly associated to the TF.
* Gene-GO links for BP (biological process).
For *A. thaliana* only those supported by experimental evidence ("EXP", "IMP", "IDA", "IPI", "IGI", "IEP") and manually curated ("TAS", "NAS", "IC") are included, while for other species all the evidences are kept.
Note: all ancestral terms are included and terms associated with more than 30% of all the genes of the species are filtered out.
* Gene-alias file downloaded from [TAIR](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FPublic_Data_Releases%2FTAIR_Data_20140331), [funRiceGenes](https://funricegenes.github.io/), [Sol Genomics Network](https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/ITAG4.0_descriptions.txt) and [MaizeGDB](https://www.maizegdb.org/associated_genes?type=all&style=tab).

```	
    doMotifAnalysis = true // set to <false> if no motif mapping data is available [CAUTION: without motif data MINI-EX is less reliable]
    featureFile_motifs = "$baseDir/data/ath/ath_2021.1_motifMapping.out.gz"
    infoTF = "$baseDir/data/ath/ath_TF2fam2mot.txt"
    GOfile = "$baseDir/data/ath/ath_full_BP_expcur_ext_names.txt"
    alias = "$baseDir/data/ath/ath_gene_aliases.txt"
```
The **third block** of parameters consists in all the scripts used in the pipeline:  
 
```	
    script_enricher = "$baseDir/bin/enricherv2.4"
    script_checkInput = "$baseDir/bin/MINIEX_checkInput.py"
    script_grnboost = "$baseDir/bin/MINIEX_grnboostMultiprocess.py"
    script_motifs = "$baseDir/bin/MINIEX_filterForMotifs.py"
    script_topDEGs = "$baseDir/bin/MINIEX_selectTopDEGs.py"
    script_expTFs = "$baseDir/bin/MINIEX_filterForTFExp.py"
    script_info = "$baseDir/bin/MINIEX_makeInfoFile.py"
    script_clustermap = "$baseDir/bin/MINIEX_clustermap.py"
    script_networkCentrality = "$baseDir/bin/MINIEX_network_analysis.py"
    script_checkReference = "$baseDir/bin/MINIEX_checkRef.py"
    script_filesEnrichment = "$baseDir/bin/MINIEX_makeFilesEnrichment.py"	
    script_makedfRef = "$baseDir/bin/MINIEX_makeRankingDf_ref.py"
    script_makedfStd = "$baseDir/bin/MINIEX_makeRankingDf_std.py"
    script_makeborda = "$baseDir/bin/MINIEX_makeBorda.py"
    script_scoreEdges = "$baseDir/bin/MINIEX_scoreEdges.py"
    script_heatmapTops = "$baseDir/bin/MINIEX_visual_heatmap_top150.py"
    script_regmaps = "$baseDir/bin/MINIEX_regmap.py"
```

The **last block** of parameters defines the filters used along the GRN inference.  
  
The first two filters (tops and expressionFilter) have been chosen by benchmarking different filters against a root gold standard of known protein-DNA interactions.  
The first refers to the number of upregulated genes per cluster (sorted by q-value) to use during the cell cluster enrichment, while the second refers to the percentage of cells that need to express the TF to retain the regulon for the specific cell cluster:    
  
* motifFilter can be set to **TF_motifs** if the user wishes not to extend the retention of regulons enriched for family motifs, but only to direct TF-motifs  
* topRegs defines the top regulons to show in the two output heatmaps. It can be changed according to the user needs

```	
    tops = "700"
    expressionFilter = "10"
    motifFilter = "TF-F_motifs"	
    topRegs = "150"
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
