#!/usr/bin/env nextflow

nextflow.enable.dsl=2
/*
=====================================================
        MINI-EX
=====================================================
Motif-Informed Network Inference from gene EXpression

-----------------------------------------------------
*/

if (params.doMotifAnalysis){
    motifLog = "Running TF motif enrichment filtering on ${params.motifFilter}"
} else {
    motifLog = "Skipping motif enrichment filtering"
}

log.info"""\
         Motif-Informed Network Inference from gene EXpression v.2.2
         ===========================================================
         ${motifLog}
         Running single-cell cluster enrichment using the top ${params.topMarkers} upregulated genes per cluster
         Filtering out regulons of single-cell clusters where the TF is expressed in less than ${params.expressionFilter} % of the cells
         Plotting expression specificity and DE calls for the top ${params.topRegulons} regulons
         """
         .stripIndent()


scriptEnricher = "$baseDir/bin/enricherv3.2.4"

grnboostOutput = "$baseDir/GRNBoost2_output"
regOutput = "$baseDir/regulons_output"
goEnrichmentOutput = "$baseDir/GOenrichment_output"
figs = "$baseDir/figures"


process check_input_files {
    input:
    path expressionMatrix
    path markersOut
    path cellsToClusters
    path clustersToIdentities
    path tfList
    path termsOfInterest
    path grnboostOut
    path featureFileMotifs
    path infoTf
    path goFile
    path geneAliases

    output:
    stdout

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_checkInput.py" "$expressionMatrix" "$markersOut" "$cellsToClusters" "$clustersToIdentities" "$tfList" "$termsOfInterest" "$grnboostOut" "$featureFileMotifs" "$infoTf" "$goFile" "$geneAliases"
    """
}


process get_expressed_genes {

    input:
    tuple val(datasetId), path(matrix)

    output:
    tuple val(datasetId), path("${datasetId}_expressedGenes.txt")

    """
    tail -n +2 "$matrix" | cut -f 1 > "${datasetId}_expressedGenes.txt"
    """
}


process run_grnboost {
    publishDir grnboostOutput, mode: 'copy'

    input:
    path tfList
    tuple val(datasetId), path(matrix)
    
    output:
    tuple val("${datasetId}"), path("${datasetId}_grnboost2.txt")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_grnboostMultiprocess.py" $tfList "$matrix" "${task.cpus}" "${datasetId}_grnboost2.txt"
    """
}

process unzip_motifMappings {
    
    input:
    path featureFileMotifs

    output:
    path "MotifMappings.csv"

    """
    gunzip -c $featureFileMotifs > "MotifMappings.csv"
    """
}

process run_enricher_motifs {
    
    input:
    path scriptEnricher
    path featureFileMotifsUnzipped
    tuple val(datasetId), path(modules), path(expressedGenes)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enricherRegulons.txt")

    """
    cat "$featureFileMotifsUnzipped" | grep -f "$expressedGenes" > featureFileMotifsfiltered.txt
    $scriptEnricher featureFileMotifsfiltered.txt "$modules" -f 0.001 -n \$(cat "$expressedGenes" | wc -l) --min-hits 2 --print-hits -o "${datasetId}_enricherRegulons.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enricherRegulons.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process filter_motifs {

    input:
    path infoTf
    val motifFilter
    tuple val(datasetId), path(enrichedModules)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enrichedRegulons.txt")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_filterForMotifs.py" $infoTf "$enrichedModules" "${datasetId}_enrichedRegulons.txt" "$motifFilter"
    """
}

process filter_motifs_dummy {

    input:
    tuple val(datasetId), path(modules)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enrichedRegulons.txt")
    
    """
    cp "$modules" "${datasetId}_enrichedRegulons.txt"
    """
}

process get_topDEGs {
    
    input:
    val topMarkers
    tuple val(datasetId), path(allMarkers)

    output:
    tuple val("${datasetId}"), path("${datasetId}_top${topMarkers}cellClusters.out")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_selectTopDEGs.py" $allMarkers "$topMarkers" "${datasetId}_top${topMarkers}cellClusters.out"
    """
}

process run_enricher_cluster {
    
    input:
    path scriptEnricher
    tuple val(datasetId), path(featureFileCellClusters), path(filteredRegulons), path(expressedGenes)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enricherCelltypes.txt") 

    """
    cat "$featureFileCellClusters" | grep -f "$expressedGenes" > feature_file_cell_clusters_filtered.txt
    $scriptEnricher feature_file_cell_clusters_filtered.txt "$filteredRegulons" -f 0.001 -n \$(cat "$expressedGenes" | wc -l) --min-hits 2 --print-hits -o "${datasetId}_enricherCelltypes.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enricherCelltypes.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process filter_expression {
    publishDir regOutput, mode: 'copy'

    input:
    path infoTf
    val expressionFilter
    tuple val(datasetId), path(expressionMatrix), path(cellClusters), path(regulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_regulons.txt")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_filterForTFExp.py" "$expressionMatrix" $infoTf $cellClusters "$expressionFilter" "$regulons" "${datasetId}_regulons.txt"
    """
}

process make_info_file {
    publishDir regOutput, mode: 'copy'

    input:
    tuple val(datasetId), path(expressionMatrix), path(grnboostRegulons), path(motifEnrichedRegulons), path(finalRegulons), path(cellClusters), path(clusterIdentities)
    path infoTf

    output:
    tuple val("${datasetId}"), path("${datasetId}_TF_info_file.txt")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeInfoFile.py" "$expressionMatrix" "$grnboostRegulons" "$motifEnrichedRegulons" "$finalRegulons" $infoTf $cellClusters $clusterIdentities "${datasetId}_TF_info_file.txt"
    """
}

process clustermap_regs {
    publishDir figs, mode: 'copy'
    
    input:
    tuple val(datasetId), path(clusterIdentities), path(finalRegulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_clustermap.svg")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_visualClustermap.py" "$clusterIdentities" "$finalRegulons" "${datasetId}_clustermap.svg"
    """
}


process network_centrality {
    
    input:
    tuple val(datasetId), path(finalRegulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_networkCentrality.txt")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_network_analysis.py" "$finalRegulons" "${datasetId}_networkCentrality.txt"
    """
}

process getFiles_enrichment {
    
    input:
    tuple val(datasetId), path(finalRegulons), path(goFile)

    output:
    tuple val("${datasetId}"), path("${datasetId}_setFileRegulons.out"), path("${datasetId}_featureFileGO.out")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeFilesEnrichment.py" "$finalRegulons" "$goFile" "${datasetId}_setFileRegulons.out" "${datasetId}_featureFileGO.out"
    """
}

process GO_enricher {
    publishDir goEnrichmentOutput, mode: 'copy'
    
    input:
    path scriptEnricher
    tuple val(datasetId), path(set), path(featureFile), path(expressedGenes)
    
    output:
    tuple val("${datasetId}"), path("${datasetId}_enricherGO.txt")
     
    """
    cat "$featureFile" | grep -f "$expressedGenes" > feature_file_filtered.txt
    $scriptEnricher feature_file_filtered.txt "$set" -f 0.05 -n \$(cat "$expressedGenes" | wc -l) --min-hits 2 -p -o "${datasetId}_enricherGO.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enricherGO.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process check_reference {
    
    input:
    tuple val(datasetId), path(finalRegulons), path(goFile)
    path termsOfInterest

    output:
    tuple val("${datasetId}"), stdout
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_checkRef.py" "$finalRegulons" "$goFile" "$termsOfInterest" 
    """
}

process ranking_df_ref {
    
    input:
    path geneAliases
    path termsOfInterest
    tuple val(datasetId), path(clusterIdentities), path(finalRegulons), path(qValueCluster), path(networkCentrality), path(goEnrichment), path(allMarkers), path(goFile)

    output:
    tuple val("${datasetId}"), path("${datasetId}_dfForRanking.txt")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeRankingDf_ref.py" "$clusterIdentities" "$finalRegulons" "$geneAliases" "$qValueCluster" "$goFile" "$termsOfInterest" "$networkCentrality" "$goEnrichment" "$allMarkers" "${datasetId}_dfForRanking.txt"
    """  
}

process ranking_df_std {
    
    input:  
    path geneAliases
    tuple val(datasetId), path(clusterIdentities), path(finalRegulons), path(qValueCluster), path(networkCentrality), path(allMarkers), val(goFile)

    output:
    tuple val("${datasetId}"), path("${datasetId}_dfForRanking.txt")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeRankingDf_std.py" "$clusterIdentities" "$finalRegulons" "$geneAliases" "$qValueCluster" "$goFile" "$networkCentrality" "$allMarkers" "${datasetId}_dfForRanking.txt"
    """  
}


process makeBorda {
    publishDir regOutput, mode: 'copy'
    echo true
    
    input:
    tuple val(datasetId), path(regulonsDataframe), val(metrics), val(ref)

    output:
    tuple val("${datasetId}"), path("${datasetId}_rankedRegulons.xlsx")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeBorda.py" "$regulonsDataframe" "$metrics" "${datasetId}_rankedRegulons.xlsx" "$ref"
    """
}


process scoreEdges {
    publishDir regOutput, mode: 'copy'
    
    input:
    tuple val(datasetId), path(regulons), path(rankedRegulons), path(modules)
    
    output:
    tuple val(datasetId), path("${datasetId}_edgeTable.tsv")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_scoreEdges.py" $regulons $rankedRegulons $modules ${datasetId}_edgeTable.tsv
    """
}


process heatmap_tops {
    publishDir figs, mode: 'copy'
    
    input:
    tuple val(datasetId), path(rankedRegulons)
    val topRegulons

    output:
    tuple val(datasetId), path("${datasetId}_heatmapSpecificity.svg"), path("${datasetId}_heatmapDEcalls.svg")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_visualHeatmapTop150.py" "$rankedRegulons" "${datasetId}_heatmapSpecificity.svg" "${datasetId}_heatmapDEcalls.svg" "$topRegulons"  
    """
}


process regmaps {
    publishDir figs, mode: 'copy'
    
    input:
    tuple val(datasetId), path(expressionMatrix), path(cellClusters), path(clusterIdentities), path(rankedRegulons)
    val topRegulons

    output:
    tuple val(datasetId), path("${datasetId}_regmap_*.svg")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_visualRegmap.py" -c $cellClusters \
                                                                    -i $clusterIdentities \
                                                                    -r $rankedRegulons \
                                                                    -m $expressionMatrix \
                                                                    -t 10,25,50,100,$topRegulons \
                                                                    -d $datasetId
    """
}


workflow {
    check_input_files(
        Channel.fromPath(params.expressionMatrix, checkIfExists:true).collect(),
        Channel.fromPath(params.markersOut, checkIfExists:true).collect(),
        Channel.fromPath(params.cellsToClusters, checkIfExists:true).collect(),
        Channel.fromPath(params.clustersToIdentities, checkIfExists:true).collect(),
        Channel.fromPath(params.tfList, checkIfExists:true).collect(),
        params.termsOfInterest != null ? Channel.fromPath(params.termsOfInterest, checkIfExists:true).collect() : "/dummy_path_tof",
        params.grnboostOut != null ? Channel.fromPath(params.grnboostOut, checkIfExists:true).collect() : "/dummy_path_grnb",
        params.doMotifAnalysis == true? Channel.fromPath(params.featureFileMotifs, checkIfExists:true) : "/dummy_path_motif",
        Channel.fromPath(params.infoTf, checkIfExists:true).collect(),
        params.goFile != null ? Channel.fromPath(params.goFile, checkIfExists:true).collect() : "/dummy_path_go",
        Channel.fromPath(params.geneAliases, checkIfExists:true).collect())

    check_input_files.out.view()

    matrix_ch = Channel.fromPath(params.expressionMatrix).map { n -> [ n.baseName.split("_")[0], n ] }
    matrix_ch.view()
    
    if (params.grnboostOut == null){
        run_grnboost(params.tfList,matrix_ch)
        grnboost_ch = run_grnboost.out
        grnboost_ch.view()
    }
    if (params.grnboostOut != null){
        grnboost_ch = Channel.fromPath(params.grnboostOut).map { n -> [ n.baseName.split("_")[0], n ] }
        grnboost_ch.view()
    }

    expressed_genes_ch = get_expressed_genes(matrix_ch)
    grnboost_combined_ch = grnboost_ch.join(expressed_genes_ch)
    
    if (params.doMotifAnalysis){
        unzip_motifMappings(params.featureFileMotifs)
        run_enricher_motifs(scriptEnricher,unzip_motifMappings.out,grnboost_combined_ch)
        
        filter_motifs(params.infoTf,params.motifFilter,run_enricher_motifs.out)
        filter_motifs_ch = filter_motifs.out
    } else {
        filter_motifs_dummy(grnboost_ch)
        filter_motifs_ch = filter_motifs_dummy.out
    }

    deg_ch = Channel.fromPath(params.markersOut).map { n -> [ n.baseName.split("_")[0], n ] }
    cluster_enrich_ch = get_topDEGs(params.topMarkers,deg_ch)

    cluster_enrich_combined_ch = cluster_enrich_ch.join(filter_motifs_ch).join(expressed_genes_ch)
    run_enricher_cluster(scriptEnricher,cluster_enrich_combined_ch)  

    cluster_ch = Channel.fromPath(params.cellsToClusters).map { n -> [ n.baseName.split("_")[0], n ] }
    filter_combined_ch = matrix_ch.join(cluster_ch).join(run_enricher_cluster.out)
    filter_expression(params.infoTf,params.expressionFilter,filter_combined_ch)
    
    cluster_ids_ch = Channel.fromPath(params.clustersToIdentities).map { n -> [ n.baseName.split("_")[0], n ] }
    info_ch = matrix_ch.join(grnboost_ch).join(filter_motifs_ch).join(filter_expression.out).join(cluster_ch).join(cluster_ids_ch)
    make_info_file(info_ch,params.infoTf)    
    
    regulons_ident_ch = cluster_ids_ch.join(filter_expression.out)
    clustermap_regs(regulons_ident_ch)
    
    network_centrality(filter_expression.out)

    if ( params.goFile != null ){

        getFiles_enrichment_input_ch = filter_expression.out.combine(Channel.fromPath(params.goFile))
        getFiles_enrichment_input_ch.view()
        getFiles_enrichment_ch = getFiles_enrichment(getFiles_enrichment_input_ch)

        getFiles_enrichment_combined_ch = getFiles_enrichment_ch.join(expressed_genes_ch)
        GO_enricher(scriptEnricher,getFiles_enrichment_combined_ch)   
    }

    if ( params.goFile != null && params.termsOfInterest != null ){        

        check_reference(getFiles_enrichment_input_ch,params.termsOfInterest)
        check_reference_trimmed = check_reference.out.map { n,it -> [ n, it.trim() ] }
    
        rankingRef_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(network_centrality.out).join(GO_enricher.out).join(deg_ch).combine(Channel.fromPath(params.goFile))
        ranking_df_ref(params.geneAliases,params.termsOfInterest,rankingRef_combined_ch)
    
        metrics_ch = Channel.value('qval_cluster,out-degree,betweenness,closeness,GO_enrich_qval')
        metrics_combi_ch = ranking_df_ref.out.combine(metrics_ch).join(check_reference_trimmed)
    
    } else {
        rankingStd_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(network_centrality.out).join(deg_ch).combine(Channel.value(false))
        ranking_df_std(params.geneAliases,rankingStd_combined_ch)    
        
        metrics_ch = Channel.value('qval_cluster,out-degree,betweenness,closeness')
        std_ch = Channel.value('std')
        metrics_combi_ch = ranking_df_std.out.combine(metrics_ch).combine(std_ch)        
    }

    makeBorda(metrics_combi_ch)
    
    scoreEdges_ch = filter_expression.out.join(makeBorda.out).join(grnboost_ch)
    scoreEdges(scoreEdges_ch)

    heatmap_tops(makeBorda.out,params.topRegulons)

    regmaps_ch = matrix_ch.join(cluster_ch).join(cluster_ids_ch).join(makeBorda.out)    
    regmaps(regmaps_ch, params.topRegulons)
}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
