# MINI-EX

Motif-Informed Network Inference based on single-cell EXpression data  

The pipeline is built using Nextflow DSL2 and has the purpose of infer cell-type specific gene regulatory network using scRNA-Seq data in plants.  
  
MINI-EX uses a GNU GENERAL PUBLIC LICENSE version 3 within a dual license to offer the distribution of the software under a proprietary model as well as an open source model.  

## **Pipeline summary**
1. Run expression-based GRN inference ([GRNBoost2](https://arboreto.readthedocs.io/en/latest/algorithms.html#grnboost2)) given a list of TFs and a gene-to-cell count matrix
2. Run TFBS enrichment on the expression-based regulons  
3. Filter the TFBS-enriched regulons for TF or TF-Family motifs (default TF-Family)
4. Filter the previously identified regulons by target genes' expression among the defined cell clusters (cell cluster enrichment)
5. Filter the cell cluster specific regulons by TF expression 
6. Calculate network centrality measures (out-degree, betweenness, closeness)
7. Calculate functional enrichment of the target genes of each regulon (if a list of expected GO terms is provided)
8. Generate a list of ranked regulons based on Borda ranking

If a list of expected GO terms is provided:
- First all the combinations of weighted metrics (4 network centrality measures, q-value from the cell cluster enrichment, q-value from the functional enrichment) are evaluated
- The combination which returns half of the expected regulons earlier in the ranks (R50) is chosen for the weighted Borda ranking

else:
- The 4 network centrality measures and q-value from the cell cluster enrichment are used to calculate the Borda ranking (caluclated on the geometric mean of the single metrics)

## **Inputs**
* Gene-to-cell count matrix (genes as rows and cells as columns)
* List of TFs
* [Seurat](https://satijalab.org/seurat/) output from [FindAllMarkers](https://www.rdocumentation.org/packages/Seurat/versions/3.1.2/topics/FindAllMarkers)
* Tab-separated file containing the cluster identity of each cell (cell_barcode \t cluster_id)
* Tab-separated file containing the cluster annotation (cluster_id \t cluster_annotation)
* (Optional) List of GO terms of interest

As the pipeline can be run in parallel for multiple datasets all the inputs can be provided as a path to the dedicated directories.  
All input files should have specific extensions and names as shown in [here](docs/data_preparation.md).  

## **Outputs**
* **regulons_output folder** containing a tab-separated file with the inferred regulons and an excel file with the ranked regulons and relative metadata
* **figures folder** containing a clustermap reporting the distribution of the different regulons across the cell clusters, and two heatmaps showing the cell cluster specificity and DE calls of the top 150 regulons, respectively. 
* **GOenrichment_output folder** containing a tab-separated file with GO enrichment for the different regulons with relative statistics
* **GRNBoost2_output folder** containing a TF-TG tab-separated file resulted from the GRNBoost2 run   

##   
A detailed overview on necessary input files and expected output files can be found [here](example/).
## 
Requirements:

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html)
## 
How to run it:

* Define paths in the [config file](docs/configuration.md) to all the required imputs

```
nextflow -C miniex.config run miniex.nf
```

##  
MINI-EX uses a GNU GENERAL PUBLIC LICENSE version 3 within a dual license  