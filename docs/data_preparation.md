# Prepare your files

Following constrains must be satisfied for the input files:
1. MINI-EX has been designed to run on multiple datasets at once. Different files from the same dataset **must have the same prefix followed by an underscore** (i.e. dataset1_matrix.txt, dataset1_allMarkers.txt, dataset1_cells2clusters.txt, dataset1_identities.txt).
2. **Underscores are not allowed** in cluster names.
3. **Row names must be unique** in all the input files.
4. **Quotes are not allowed** in the input files.
5. **Cluster identifiers must correspond** between the "*_cells2clusters" and "*_identities" files.
6. **Only up-regulated markers** can be provided in the "*_allMarkers" file.

A detailed example of the necessary input files can be found [here](/example/).  
  
The paths to the folders containing the different input files need to be stated in the [config file](/docs/configuration.md).  

**Note:** if GRNBoost2 was previously run separately, its output can be specified in the configuration file and the corresponding MINI-EX's step will be skipped.   
```
params {

    // INPUT FILES DERIVED FROM SINGLE-CELL DATASETS //
    expressionMatrix = "$baseDir/example/INPUTS/*_matrix.tsv"
    markersOut = "$baseDir/example/INPUTS/*_allMarkers.tsv"
    cellsToClusters = "$baseDir/example/INPUTS/*_cells2clusters.tsv"
    clustersToIdentities = "$baseDir/example/INPUTS/*_identities.tsv"
    
    // CACHED FILES // (allows to reuse the initial coexpression network)
    grnboostOut = "/$baseDir/example/OUTPUTS/grnboost2/*_grnboost2.tsv" 
    // --> to infer coexpression network de novo, replace the line above by: grnboostOut = null
```

The **expressionMatrix** points to the gene-to-cell count matrix and can be extracted from the [Seurat](https://satijalab.org/seurat/) object using the command below:

```
expression.matrix <- as.data.frame(as.matrix(GetAssayData(object = SEURAT_OBJECT, assay = "RNA", slot = "counts")))
write.table(expression.matrix, "./dataset1_matrix.tsv", sep='\t', quote = FALSE)
```

```
geneid	CTACGTCTCGGATGGA-1	TCGGGACGTCGAATCT-1	CGGACTGGTGATGATA-1	GTTCGGGTCACCATAG-1	GACTAACTCGCTGATA-1	GACTAACGTCTCCATC-1	GCCTCTATCTCGCATC-1	TACTCATAGCAATATG-1	ACTGAACAGTTAAGTG-1
AT1G67265	10	0	0	3	0	0	7	0	1
AT3G12100	0	0	0	1	0	3	0	1	18
AT4G03340	14	0	0	0	0	0	0	0	0
AT2G28360	0	2	0	0	25	0	0	0	0
```

The **markersOut** points to the output obtained by Seurat [FindAllMarkers](https://www.rdocumentation.org/packages/Seurat/versions/3.1.2/topics/FindAllMarkers) using the command below:  

```
ath_markers_all <- FindAllMarkers(SEURAT_OBJECT, only.pos=TRUE)
write.table(ath_markers_all, "./dataset1_allMarkers.tsv", sep = "\t", quote=FALSE)
```

```
geneID	p_val	avg_logFC	pct.1	pct.2	p_val_adj	cluster	gene
AT4G26620	1.7566143109614e-70	0.273351093463235	0.532	0.214	3.9999864474902e-66	13	AT4G26620
AT3G11400	1.77798086156283e-70	0.529710021315979	0.774	0.482	4.04864021986471e-66	13	AT3G11400
AT3G48380	3.94993156215768e-70	0.442472302557895	0.658	0.343	8.99438916018925e-66	13	AT3G48380
AT1G20330	4.26798715389231e-70	0.5128778779104	0.79	0.501	9.71863354812817e-66	13	AT1G20330
AT1G67030	4.34924991762384e-70	0.462602175019394	0.706	0.363	9.90367698742125e-66	13	AT1G67030
AT3G02360	5.20422679108741e-12	0.439837431384755	0.226	0.54	1.18505448259851e-07	14	AT3G02360
AT4G03110	6.13871513228955e-12	0.380868416173313	0.084	0.274	1.39784682277365e-07	14	AT4G03110
AT1G08190	6.19518803968411e-12	0.323064627623938	0.072	0.255	1.41070626851647e-07	14	AT1G08190
AT5G15470	6.73638999455266e-12	0.323992721599046	0.054	0.223	1.53394336565959e-07	14	AT5G15470
```

The **cellsToClusters** points to a tab-separated file containing the cluster annotation for each cell in the expression matrix. It can be obtained using the [Seurat](https://satijalab.org/seurat/) command [FetchData](https://www.rdocumentation.org/packages/Seurat/versions/3.1.2/topics/FetchData) as shown below:  

```
object <-FetchData(SEURAT_OBJECT, vars = 'ident')
write.table(object,"./dataset1_cell2cluster.tsv", sep='\t',quote=FALSE,col.names = FALSE) 
```

```
CTACGTCTCGGATGGA-1	1
TCGGGACGTCGAATCT-1	4
CGGACTGGTGATGATA-1	31
GTTCGGGTCACCATAG-1	5
GACTAACTCGCTGATA-1	10
GACTAACGTCTCCATC-1	9
GCCTCTATCTCGCATC-1	20
TACTCATAGCAATATG-1	10
ACTGAACAGTTAAGTG-1	20
```

The **clustersToIdentities** points to a tab-separated file containing the cell type annotation for each cluster. Optionally, this file can contain a third column which specifies a cluster index. This index is used to indicate the position of a cluster along a known developmental trajectory (see [miniexExample_identities_with_idx.tsv](/example/INPUTS/miniexExample_identities_with_idx.tsv) in the [INPUTS](/example/INPUTS) folder), and translates to the column index in the regulator heatmaps ([example](/example/OUTPUTS/figures/miniexExample_regmap_8.svg)). If this column is not specified, an automatic index is created by sorting the cluster identities alphabetically.

```
13	cortex
14	endodermis
15	epidermis
16.0	pericycle
```

The **TF_list** is a list of TFs which is used in the GRNBoost2 run.  
  
The [TF_list.txt](/example/INPUTS/TF_list.tsv) contained in the [INPUTS](/example/INPUTS) folder contains 1879 TFs collected from [PlantRegMap/PlantTFDB v5.0](http://planttfdb.gao-lab.org/), [PlnTFDB v3.0](http://plntfdb.bio.uni-potsdam.de/v3.0/) and [TF2Network](http://bioinformatics.psb.ugent.be/webtools/TF2Network/) for which either direct TF-motif information was available or motif information related to the TF family.

```
AT1G18790
AT1G18835
AT1G18860
AT1G19040
AT1G19350
AT1G19490
AT1G19790
AT1G20640
AT1G20980
AT1G21200
AT1G21340
```

The **termsOfInterest** parameter can be set to null - if we do not expect specific functions popping up in our regulons - or alternatively a list of expected terms can be provided.  
  
In the example dataset we analyzed Arabidopsis's root samples and therefore we looked for terms related to root and vasculature:

```
xylem
phloem
trichome
root 
vasculature
trichoblast
vascular
stele
tracheary
procambium
sieve
```
