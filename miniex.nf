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

miniexVersion = "v2.2"
log.info"""\
         Motif-Informed Network Inference from gene EXpression ${miniexVersion}
         ===========================================================
         ${motifLog}
         Running single-cell cluster enrichment using the top ${params.topMarkers} upregulated genes per cluster
         Filtering out regulons of single-cell clusters where the TF is expressed in less than ${params.expressionFilter} % of the cells
         Plotting expression specificity and DE calls for the top ${params.topRegulons} regulons
         """
         .stripIndent()


scriptEnricher = "$baseDir/bin/enricher_v.3.3.1"

grnboostDir = "$params.outputDir/grnboost2"
regulonsDir = "$params.outputDir/regulons"
goEnrichmentDir = "$params.outputDir/go_enrichment"
figuresDir = "$params.outputDir/figures"
logDir = "$params.outputDir"


process check_user_input {

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
    path enrichmentBackground
    val doMotifAnalysis
    val topMarkers
    val expressionFilter
    val motifFilter
    val topRegulons

    output:
    stdout emit: stdoutLog
    path("processLog.log"), emit: processLog

    """
    echo -n "MINI-EX" "$miniexVersion" "\nPipeline started on: " > "processLog.log"
    date >> "processLog.log"
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_checkUserInput.py" "$expressionMatrix" "$markersOut" "$cellsToClusters" "$clustersToIdentities" \
                                                                  "$tfList" "$termsOfInterest" "$grnboostOut" "$featureFileMotifs" "$infoTf" \
                                                                  "$goFile" "$geneAliases" "$enrichmentBackground" "$doMotifAnalysis" "$topMarkers" \
                                                                   "$expressionFilter" "$motifFilter" "$topRegulons" >> "processLog.log"
    # print input validation statistics on stdout
    awk '/== INPUT VALIDATION/,/== INPUT FILES/ {if (!/== INPUT (VALIDATION|FILES)/) print}' processLog.log
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
    publishDir grnboostDir, mode: 'copy'

    input:
    path tfList
    tuple val(datasetId), path(matrix)
    
    output:
    tuple val("${datasetId}"), path("${datasetId}_grnboost2.tsv")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_runGrnboost.py" $tfList "$matrix" "${task.cpus}" "${datasetId}_grnboost2.tsv"
    """
}


process unzip_motif_mappings {
    
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
    tuple val("${datasetId}"), path("${datasetId}_enrichedRegulons.txt")

    """
    $scriptEnricher "$featureFileMotifsUnzipped" "$modules" -f 0.001 -b "$expressedGenes" --min-hits 2 --print-hits -o "${datasetId}_enrichedRegulons.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enrichedRegulons.txt | grep -P '^[^#]' | wc -l)" -eq 1 ]; then 
        echo 'ERROR: Enricher output is empty!'
        exit 1
    fi
    """
}


process filter_motifs {

    input:
    path infoTf
    val motifFilter
    tuple val(datasetId), path(enrichedRegulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enrichedRegulonsFiltered.txt")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_filterMotifs.py" $infoTf "$enrichedRegulons" "${datasetId}_enrichedRegulonsFiltered.txt" "$motifFilter"
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


process get_top_degs {
    
    input:
    val topMarkers
    tuple val(datasetId), path(allMarkers)

    output:
    tuple val("${datasetId}"), path("${datasetId}_top${topMarkers}cellClusters.out")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_getTopDegs.py" $allMarkers "$topMarkers" "${datasetId}_top${topMarkers}cellClusters.out"
    """
}


process run_enricher_cluster {
    
    input:
    path scriptEnricher
    tuple val(datasetId), path(featureFileCellClusters), path(filteredRegulons), path(expressedGenes)

    output:
    tuple val("${datasetId}"), path("${datasetId}_enrichedCelltypes.txt") 

    """
    $scriptEnricher "$featureFileCellClusters" "$filteredRegulons" -f 0.001 -b "$expressedGenes" --min-hits 2 --print-hits -o "${datasetId}_enrichedCelltypes.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enrichedCelltypes.txt | grep -P '^[^#]' | wc -l)" -eq 1 ]; then 
        echo 'ERROR: Enricher output is empty!'
        exit 1
    fi
    """
}


process filter_expression {
    publishDir regulonsDir, mode: 'copy'

    input:
    path infoTf
    val expressionFilter
    tuple val(datasetId), path(expressionMatrix), path(cellClusters), path(regulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_regulons.tsv")
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_filterExpression.py" "$expressionMatrix" $infoTf $cellClusters "$expressionFilter" "$regulons" "${datasetId}_regulons.tsv"
    """
}


process make_info_file {
    publishDir regulonsDir, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(datasetId), path(expressionMatrix), path(grnboostRegulons), path(motifEnrichedRegulons), path(finalRegulons), path(cellClusters), path(clusterIdentities)
    path tfList

    output:
    tuple val("${datasetId}"), path("${datasetId}_TF_info_file.tsv")
    tuple val("${datasetId}"), path("${datasetId}_regulonInfoLog.log"), emit: processLog
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeInfoFile.py" "$expressionMatrix" "$grnboostRegulons" "$motifEnrichedRegulons" "$finalRegulons" $tfList $cellClusters $clusterIdentities "${datasetId}_TF_info_file.tsv" > "${datasetId}_regulonInfoLog.log"
    """
}


process make_regulon_clustermap {
    publishDir figuresDir, mode: 'copy'
    
    input:
    tuple val(datasetId), path(clusterIdentities), path(finalRegulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_clustermap.svg"), path("${datasetId}_clustermap.png")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeRegulonClustermap.py" "$clusterIdentities" "$finalRegulons" "${datasetId}_clustermap"
    """
}


process get_network_centrality {
    
    input:
    tuple val(datasetId), path(finalRegulons)

    output:
    tuple val("${datasetId}"), path("${datasetId}_networkCentrality.txt")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_getNetworkCentrality.py" "$finalRegulons" "${datasetId}_networkCentrality.txt"
    """
}


process make_go_enrichment_files {
    
    input:
    tuple val(datasetId), path(finalRegulons), path(goFile)

    output:
    tuple val("${datasetId}"), path("${datasetId}_setFileRegulons.out"), path("${datasetId}_featureFileGO.out")

    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeGoEnrichmentFile.py" "$finalRegulons" "$goFile" "${datasetId}_setFileRegulons.out" "${datasetId}_featureFileGO.out"
    """
}


process run_enricher_go {
    publishDir goEnrichmentDir, mode: 'copy'
    
    input:
    path scriptEnricher
    tuple val(datasetId), path(set), path(featureFile), path(expressedGenes)
    
    output:
    tuple val("${datasetId}"), path("${datasetId}_enricherGO.txt")
     
    """
    $scriptEnricher "$featureFile" "$set" -f 0.05 -b "$expressedGenes" --min-hits 2 -p -o "${datasetId}_enricherGO.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${datasetId}_enricherGO.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}


process select_borda_procedure {
    
    input:
    tuple val(datasetId), path(finalRegulons), path(goFile)
    path termsOfInterest

    output:
    tuple val("${datasetId}"), stdout
    
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_selectBordaProcedure.py" "$finalRegulons" "$goFile" "$termsOfInterest" 
    """
}


process make_ranking_dataframe {
    
    input:  
    path geneAliases
    tuple val(datasetId), path(clusterIdentities), path(finalRegulons), path(regulonEnrichment), path(networkCentrality), path(allMarkers), path(goEnrichment)
    path goAnnotations
    path termsOfInterest

    output:
    tuple val("${datasetId}"), path("${datasetId}_dfForRanking.txt")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeRankingDataframe.py" "$geneAliases" "$clusterIdentities" "$finalRegulons" "$regulonEnrichment" "$networkCentrality" "$allMarkers" "$goEnrichment" "$goAnnotations" "$termsOfInterest" "${datasetId}_dfForRanking.txt"
    """  
}


process make_borda {
    publishDir regulonsDir, mode: 'copy', pattern: '*.{xlsx,tsv}'
    
    input:
    tuple val(datasetId), path(regulonsDataframe), val(ref)

    output:
    tuple val("${datasetId}"), path("${datasetId}_rankedRegulons.xlsx"), emit: processOut
    tuple val("${datasetId}"), path("${datasetId}_rankedRegulons.tsv")
    tuple val("${datasetId}"), path("${datasetId}_bordaProcessLog.log"), emit: processLog
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeBorda.py" "$regulonsDataframe" "${datasetId}_rankedRegulons" "$ref" > "${datasetId}_bordaProcessLog.log"
    """
}


process score_edges {
    publishDir regulonsDir, mode: 'copy'
    
    input:
    tuple val(datasetId), path(regulons), path(rankedRegulons), path(modules)
    
    output:
    tuple val(datasetId), path("${datasetId}_edgeTable.tsv")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_scoreEdges.py" $regulons $rankedRegulons $modules ${datasetId}_edgeTable.tsv
    """
}


process make_top_regulons_heatmaps {
    publishDir figuresDir, mode: 'copy'
    
    input:
    tuple val(datasetId), path(rankedRegulons)
    val topRegulons

    output:
    tuple val(datasetId), path("${datasetId}_heatmapSpecificity.svg"), path("${datasetId}_heatmapSpecificity.png"), path("${datasetId}_heatmapDEcalls.svg"), path("${datasetId}_heatmapDEcalls.png")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeTopRegulonsHeatmap.py" "$rankedRegulons" "${datasetId}_heatmapSpecificity" "${datasetId}_heatmapDEcalls" "$topRegulons"
    """
}


process make_regmaps {
    publishDir figuresDir, mode: 'copy'
    
    input:
    tuple val(datasetId), path(expressionMatrix), path(cellClusters), path(clusterIdentities), path(rankedRegulons)
    val topRegulons

    output:
    tuple val(datasetId), path("${datasetId}_regmap_*.svg"), path("${datasetId}_regmap_*.png")
     
    """
    OMP_NUM_THREADS=1 python3 "$baseDir/bin/MINIEX_makeRegmaps.py" -c $cellClusters \
                                                                   -i $clusterIdentities \
                                                                   -r $rankedRegulons \
                                                                   -m $expressionMatrix \
                                                                   -t 10,25,50,100,$topRegulons \
                                                                   -d $datasetId
    """
}


process make_log_file {
    publishDir logDir, mode: 'copy'

    input:
    tuple path(checkInputLog), val(datasetId), path(regulonInfoLog), path(bordaLog)

    output:
    path("${datasetId}_log.txt")

    """
    cat "$checkInputLog" "$regulonInfoLog" "$bordaLog" > "${datasetId}_log.txt"
    echo -n "Pipeline ended on: " >> "${datasetId}_log.txt"
    date >> "${datasetId}_log.txt"
    """
}


workflow {
    // handle files with null values (cannot be provided as parameters to a Nextflow process -> replaced with dummy names)
    terms_of_interest_file     = params.termsOfInterest      != null ? Channel.fromPath(params.termsOfInterest, checkIfExists:true).collect()      : "/.dummy_path_terms_of_interest"
    grnboost_file              = params.grnboostOut          != null ? Channel.fromPath(params.grnboostOut, checkIfExists:true).collect()          : "/.dummy_path_grnboost"
    motif_mapping_file         = params.doMotifAnalysis      == true ? Channel.fromPath(params.featureFileMotifs, checkIfExists:true).collect()    : "/.dummy_path_motif_mapping"
    go_file                    = params.goFile               != null ? Channel.fromPath(params.goFile, checkIfExists:true).collect()               : "/.dummy_path_go_annotations"
    gene_aliases_file          = params.geneAliases          != null ? Channel.fromPath(params.geneAliases, checkIfExists:true).collect()          : "/.dummy_path_gene_aliases"
    enrichment_background_file = params.enrichmentBackground != null ? Channel.fromPath(params.enrichmentBackground, checkIfExists:true).collect() : "/.dummy_path_enrichment_background"

    check_user_input(
        Channel.fromPath(params.expressionMatrix, checkIfExists:true).collect(),
        Channel.fromPath(params.markersOut, checkIfExists:true).collect(),
        Channel.fromPath(params.cellsToClusters, checkIfExists:true).collect(),
        Channel.fromPath(params.clustersToIdentities, checkIfExists:true).collect(),
        Channel.fromPath(params.tfList, checkIfExists:true).collect(),
        terms_of_interest_file,
        grnboost_file,
        motif_mapping_file,
        Channel.fromPath(params.infoTf, checkIfExists:true).collect(),
        go_file,
        gene_aliases_file,
        enrichment_background_file,
        params.doMotifAnalysis,
        params.topMarkers,
        params.expressionFilter,
        params.motifFilter,
        params.topRegulons)

    check_user_input.out.stdoutLog.view()  // print the output of check_user_input to the terminal

    matrix_ch = Channel.fromPath(params.expressionMatrix).map { n -> [ n.baseName.split("_")[0], n ] }
    
    if (params.grnboostOut == null){
        run_grnboost(params.tfList,matrix_ch)
        grnboost_ch = run_grnboost.out
    }
    if (params.grnboostOut != null){
        grnboost_ch = Channel.fromPath(params.grnboostOut).map { n -> [ n.baseName.split("_")[0], n ] }
    }

    if (params.enrichmentBackground == null){
        // use expressed genes as the enrichment background, if it is not specified by the user
        enrichment_background_ch = get_expressed_genes(matrix_ch)
        grnboost_combined_ch = grnboost_ch.join(enrichment_background_ch)
    } else {
        enrichment_background_ch = Channel.fromPath(params.enrichmentBackground)
        grnboost_combined_ch = grnboost_ch.combine(enrichment_background_ch)
    }
    
    if (params.doMotifAnalysis){
        unzip_motif_mappings(params.featureFileMotifs)
        run_enricher_motifs(scriptEnricher,unzip_motif_mappings.out,grnboost_combined_ch)
        
        filter_motifs(params.infoTf,params.motifFilter,run_enricher_motifs.out)
        filter_motifs_ch = filter_motifs.out
    } else {
        filter_motifs_dummy(grnboost_ch)
        filter_motifs_ch = filter_motifs_dummy.out
    }

    deg_ch = Channel.fromPath(params.markersOut).map { n -> [ n.baseName.split("_")[0], n ] }
    cluster_enrich_ch = get_top_degs(params.topMarkers,deg_ch)

    // add enrichment background (use join in case of expressed genes (as dataset dependent) or combine in case of a user specified enrichment background)
    if (params.enrichmentBackground == null) { cluster_enrich_combined_ch = cluster_enrich_ch.join(filter_motifs_ch).join(enrichment_background_ch) }
    else { cluster_enrich_combined_ch = cluster_enrich_ch.join(filter_motifs_ch).combine(enrichment_background_ch) }

    run_enricher_cluster(scriptEnricher,cluster_enrich_combined_ch)  

    cluster_ch = Channel.fromPath(params.cellsToClusters).map { n -> [ n.baseName.split("_")[0], n ] }
    filter_combined_ch = matrix_ch.join(cluster_ch).join(run_enricher_cluster.out)
    filter_expression(params.tfList,params.expressionFilter,filter_combined_ch)
    
    cluster_ids_ch = Channel.fromPath(params.clustersToIdentities).map { n -> [ n.baseName.split("_")[0], n ] }
    info_ch = matrix_ch.join(grnboost_ch).join(filter_motifs_ch).join(filter_expression.out).join(cluster_ch).join(cluster_ids_ch)
    make_info_file(info_ch,params.tfList)    
    
    regulons_ident_ch = cluster_ids_ch.join(filter_expression.out)
    make_regulon_clustermap(regulons_ident_ch)
    
    get_network_centrality(filter_expression.out)

    if ( params.goFile != null ){
        enrichment_files_input_ch = filter_expression.out.combine(Channel.fromPath(params.goFile))
        make_go_enrichment_files_ch = make_go_enrichment_files(enrichment_files_input_ch)

        // add enrichment background (use join in case of expressed genes (as dataset dependent) or combine in case of a user specified enrichment background)
        if (params.enrichmentBackground == null) { make_go_enrichment_files_combined_ch = make_go_enrichment_files_ch.join(enrichment_background_ch)}
        else { make_go_enrichment_files_combined_ch = make_go_enrichment_files_ch.combine(enrichment_background_ch) }

        run_enricher_go(scriptEnricher,make_go_enrichment_files_combined_ch)
        ranking_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(get_network_centrality.out).join(deg_ch).join(run_enricher_go.out)
    }
    else {
        // no GO enrichment was performed: add dummy paths for GO-related files
        enrichment_files_input_ch = filter_expression.out.combine(Channel.fromPath("/.dummy_path_go_annotations"))
        ranking_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(get_network_centrality.out).join(deg_ch).combine(Channel.fromPath("/.dummy_path_go_enrichment"))
    }
    
    make_ranking_dataframe(gene_aliases_file, ranking_combined_ch, go_file, terms_of_interest_file)

    select_borda_procedure(enrichment_files_input_ch, terms_of_interest_file)
    select_borda_procedure_trimmed = select_borda_procedure.out.map { n,it -> [ n, it.trim() ] }
    borda_input_ch = make_ranking_dataframe.out.join(select_borda_procedure_trimmed)

    make_borda(borda_input_ch)
    
    score_edges_ch = filter_expression.out.join(make_borda.out.processOut).join(grnboost_ch)
    score_edges(score_edges_ch)

    make_top_regulons_heatmaps(make_borda.out.processOut,params.topRegulons)

    make_regmaps_input_ch = matrix_ch.join(cluster_ch).join(cluster_ids_ch).join(make_borda.out.processOut)    
    make_regmaps(make_regmaps_input_ch, params.topRegulons)

    make_log_file(check_user_input.out.processLog.combine(make_info_file.out.processLog.join(make_borda.out.processLog)))
}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
