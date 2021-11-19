# Configuring MINI-EX

The config file of MINI-EX is composed by **three main sections**:

## **Exectutor scope**
MINI-EX is built to run on a SGE computer cluster as specified in the executor configuration scope at the top of the config file.  
It uses a Singularity container built on Docker image with the necessary Python 3.6 moduels.

```
executor {
    name = 'sge'
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
//	termsOfInterest = null
```

The **second block** of parameters is composed by files provided and necessary for the pipeline to run.  
There are currently three directory, one for each species supported by MINI-EX: *Arabidopsis thaliana* (data_ath), *Oryza sativa* (data_osa), and *Zea mays* (data_zma).

This consists in:
* Ensemble motif-mapping file obtained by selecting the top motif matches for [Cluster Buster](https://github.com/weng-lab/cluster-buster) and [FIMO](https://meme-suite.org/meme/doc/fimo.html) matches for a collection of motifs mapped to the regulatory regions of the species   
* TF-motif Family file linking each TF to its TF family and reporting wheather direct motif information is available or not, and the motifs directly associated to the TF
* gene-GO links for BP (biological process) supported by experimental evidence ("EXP", "IMP", "IDA", "IPI", "IGI", "IEP") and manually curated ("TAS", "NAS", "IC")  
Note: all ancestral terms are included
* gene-alias file downloaded from [TAIR](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FPublic_Data_Releases%2FTAIR_Data_20140331) and [funRiceGenes](https://funricegenes.github.io/)

```	
	featureFile_motifs = "$baseDir/data_ath/ath_featureFile_ensemble_top4000top7000.out.gz"
	infoTF = "$baseDir/data_ath/ath_TF2fam2mot.txt"
	GOfile = "$baseDir/data_ath/ath_full_BP_expcur_ext_names.txt"
	alias = "$baseDir/data_ath/ath_gene_aliases.txt"
```
The **third block** of parameters consists in all the scripts used in the pipeline  
 
```	
	script_enricher = "$baseDir/bin/enricherv2.4"
	script_grnboost = "$baseDir/bin/MINIEX_grnboostMultiprocess.py"
	script_motifs = "$baseDir/bin/MINIEX_filterForMotifs.py"
	script_topDEGs = "$baseDir/bin/MINIEX_selectTopDEGs.py"
	script_expTFs = "$baseDir/bin/MINIEX_filterForTFExp.py"
	script_clustermap = "$baseDir/bin/MINIEX_clustermap.py"
	script_networkCentrality = "$baseDir/bin/MINIEX_network_analysis.py"
	script_filesEnrichment = "$baseDir/bin/MINIEX_makeFilesEnrichment.py"	
	script_makedfRef = "$baseDir/bin/MINIEX_makeRankingDf_ref.py"
	script_makedfStd = "$baseDir/bin/MINIEX_makeRankingDf_std.py"
	script_makebordaRef = "$baseDir/bin/MINIEX_makeBorda_ref.py"
	script_makebordaStd = "$baseDir/bin/MINIEX_makeBorda_std.py"
	script_heatmapTops = "$baseDir/bin/MINIEX_visual_heatmap_top150.py"
```

The **last block** of parameters defines the filters used along the GRN inferece.  
  
The first two filters (tops and expressionFilter) have been chosen by benchmarking different filters against a root gold standard of known protein-DNA interactions. The first refers to number of upregulated genes (sorted by q-value) per cluster to use during the cell cluster enrichment, while the second refers to the percentage of cells that need to express the TF to retain the regulon for the specific cell cluster    
  
* motifFilter can be set to **TF_motifs** if the user wishes not to extend the retention of regulons enriched for family motifs, but only to direct TF-motifs.  
* topRegs defines the top regulons to show in the two output heatmaps. It can be changed according to the user needs.

```	
	tops = "700"
	expressionFilter = "10"
	motifFilter = "TF-F_motifs"	
	topRegs = "150"
}
```

## **Process scope**
The process scope allows the user to define the configuration (memory, parallel environment) for each of the processes run by MINI-EX.  
  
These can be changed according to the resources availabe/needed (i.e. below change <pe_name> with the parallel environment name of your system). 
   
```	
process {
    withName: run_grnboost {
        clusterOptions = '-l h_vmem=15G -pe <pe_name> 5'
    }
	...
}
```