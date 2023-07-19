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
    motif_log = "Running TF motif enrichment filtering on ${params.motifFilter}"
} else {
    motif_log = "Skipping motif enrichment filtering"
}

log.info"""\
         Motif-Informed Network Inference from gene EXpression v.2.0
         ===========================================================
         ${motif_log}
         Running single-cell cluster enrichment using the top ${params.tops} upregulated genes per cluster
         Filtering out regulons of single-cell clusters where the TF is expressed in less than ${params.expressionFilter} % of the cells
         Plotting expression specificity and DE calls for the top ${params.topRegs} regulons
         """
         .stripIndent()
         
grnboostOutput = "$baseDir/GRNBoost2_output"
regOutput = "$baseDir/regulons_output"
GOenrichmentOutput = "$baseDir/GOenrichment_output"
figs = "$baseDir/figures"


process check_input_files {
    echo true

    input:
    path script_checkInput
    path expressionMatrix
    path markersOut
    path cells2clusters
    path cluster2ident
    path tfList
    path termsOfInterest
    path grnboostOut
    path featureFileMotifs
    path infoTf
    path goFile
    path alias

    """
    OMP_NUM_THREADS=1 python3 $script_checkInput "$expressionMatrix" "$markersOut" "$cells2clusters" "$cluster2ident" "$tfList" "$termsOfInterest" "$grnboostOut" "$featureFileMotifs" "$infoTf" "$goFile" "$alias"
    """
}


process get_expressed_genes {

    input:
    tuple val(dataset_id), path(matrix)

    output:
    tuple val(dataset_id), path("${dataset_id}_expressedGenes.txt")

    """
    tail -n +2 "$matrix" | cut -f 1 > "${dataset_id}_expressedGenes.txt"
    """
}


process run_grnboost {
    publishDir grnboostOutput, mode: 'copy'

    input:
    path script_grnboost
    path TF_list
    tuple val(dataset_id), path(matrix)
    
    output:
    tuple val("${dataset_id}"), path("${dataset_id}_grnboost2.txt")

    """
    OMP_NUM_THREADS=1 python3 $script_grnboost $TF_list "$matrix" 5 "${dataset_id}_grnboost2.txt"
    """
}

process unzip_motifMappings {
    
    input:
    path ff_motifs

    output:
    path "MotifMappings.csv"

    """
    gunzip -c $ff_motifs > "MotifMappings.csv"
    """
}

process run_enricher_motifs {
    
    input:
    path script_enricher
    path ff_motifs_unzipped
    tuple val(dataset_id), path(modules), path(expressed_genes)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_enricherRegulons.txt")

    """
    cat "$ff_motifs_unzipped" | grep -f "$expressed_genes" > ff_motifs_filtered.txt
    $script_enricher ff_motifs_filtered.txt "$modules" -f 0.001 -n \$(cat "$expressed_genes" | wc -l) --min-hits 2 --print-hits -o "${dataset_id}_enricherRegulons.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${dataset_id}_enricherRegulons.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process filter_motifs {

    input:
    path script_motifs
    path infoTF
    val motifFilter
    tuple val(dataset_id), path(enrichedModules)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_enrichedRegulons.txt")
    
    """
    OMP_NUM_THREADS=1 python3 $script_motifs $infoTF "$enrichedModules" "${dataset_id}_enrichedRegulons.txt" "$motifFilter"
    """
}

process filter_motifs_dummy {

    input:
    tuple val(dataset_id), path(modules)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_enrichedRegulons.txt")
    
    """
    cp "$modules" "${dataset_id}_enrichedRegulons.txt"
    """
}

process get_topDEGs {
    
    input:
    path script_topDEGs
    val tops
    tuple val(dataset_id), path(allMarkers)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_top${tops}cellClusters.out")
    
    """
    OMP_NUM_THREADS=1 python3 $script_topDEGs $allMarkers "$tops" "${dataset_id}_top${tops}cellClusters.out"
    """
}

process run_enricher_cluster {
    
    input:
    path script_enricher
    tuple val(dataset_id), path(ff_celltypes), path(filteredRegulons), path(expressed_genes)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_enricherCelltypes.txt") 

    """
    cat "$ff_celltypes" | grep -f "$expressed_genes" > ff_celltypes_filtered.txt
    $script_enricher ff_celltypes_filtered.txt "$filteredRegulons" -f 0.001 -n \$(cat "$expressed_genes" | wc -l) --min-hits 2 --print-hits -o "${dataset_id}_enricherCelltypes.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${dataset_id}_enricherCelltypes.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process filter_expression {
    publishDir regOutput, mode: 'copy'

    input:
    path script_expTFs
    path infoTF
    val expressionFilter
    tuple val(dataset_id), path(expMatrix), path(cellClusters), path(regulons)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_regulons.txt")
    
    """
    OMP_NUM_THREADS=1 python3 $script_expTFs "$expMatrix" $infoTF $cellClusters "$expressionFilter" "$regulons" "${dataset_id}_regulons.txt"
    """
}

process make_info_file {
    publishDir regOutput, mode: 'copy'

    input:
    path script_info
    tuple val(dataset_id), path(expMatrix), path(grnboostRegulons), path(motenrichRegulons), path(finalRegulons), path(cellClusters), path(identClust)
    path infoTF

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_TF_info_file.txt")
    
    """
    OMP_NUM_THREADS=1 python3 $script_info "$expMatrix" "$grnboostRegulons" "$motenrichRegulons" "$finalRegulons" $infoTF $cellClusters $identClust "${dataset_id}_TF_info_file.txt"
    """
}

process clustermap_regs {
    publishDir figs, mode: 'copy'
    
    input:
    path script_clustermap
    tuple val(dataset_id), path(identClust), path(finalRegulons)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_clustermap.svg")
     
    """
    OMP_NUM_THREADS=1 python3 $script_clustermap "$identClust" "$finalRegulons" "${dataset_id}_clustermap.svg"
    """
}


process network_centrality {
    
    input:
    path script_networkCentrality
    tuple val(dataset_id), path(finalRegulons)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_networkCentrality.txt")

    """
    OMP_NUM_THREADS=1 python3 $script_networkCentrality "$finalRegulons" "${dataset_id}_networkCentrality.txt"
    """
}

process getFiles_enrichment {
    
    input:
    path script_filesEnrichment
    tuple val(dataset_id), path(finalRegulons), path(GOfile)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_setFileRegulons.out"), path("${dataset_id}_featureFileGO.out")

    """
    OMP_NUM_THREADS=1 python3 $script_filesEnrichment "$finalRegulons" "$GOfile" "${dataset_id}_setFileRegulons.out" "${dataset_id}_featureFileGO.out"
    """
}

process GO_enricher {
    publishDir GOenrichmentOutput, mode: 'copy'
    
    input:
    path script_enricher
    tuple val(dataset_id), path(s), path(f), path(expressed_genes)
    
    output:
    tuple val("${dataset_id}"), path("${dataset_id}_enricherGO.txt")
     
    """
    cat "$f" | grep -f "$expressed_genes" > f_filtered.txt
    $script_enricher f_filtered.txt "$s" -f 0.05 -n \$(cat "$expressed_genes" | wc -l) --min-hits 2 -p -o "${dataset_id}_enricherGO.txt"

    # Throw an error if no enrichment is found
    if [ "\$(head ${dataset_id}_enricherGO.txt | grep -P '^[^#]' | wc -l)" -eq 0 ]; then 
        echo 'ERROR: Enricher output empty!'
        exit 1
    fi
    """
}

process check_reference {
    
    input:
    path script_checkReference
    tuple val(dataset_id), path(finalRegulons), path(GOfile)
    path termsOfInterest

    output:
    tuple val("${dataset_id}"), stdout
    
    """
    OMP_NUM_THREADS=1 python3 $script_checkReference "$finalRegulons" "$GOfile" "$termsOfInterest" 
    """
}

process ranking_df_ref {
    
    input:
    path script_makedfRef
    path alias
    path termsOfInterest
    tuple val(dataset_id), path(identClust), path(finalRegulons), path(qvalct), path(netCent), path(enrichGO), path(allMarkers), path(GOfile)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_dfForRanking.txt")
     
    """
    OMP_NUM_THREADS=1 python3 $script_makedfRef "$identClust" "$finalRegulons" "$alias" "$qvalct" "$GOfile" "$termsOfInterest" "$netCent" "$enrichGO" "$allMarkers" "${dataset_id}_dfForRanking.txt"
    """  
}

process ranking_df_std {
    
    input:
    path script_makedfStd   
    path alias
    tuple val(dataset_id), path(identClust), path(finalRegulons), path(qvalct), path(netCent), path(allMarkers), val(GOfile)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_dfForRanking.txt")
     
    """
    OMP_NUM_THREADS=1 python3 $script_makedfStd "$identClust" "$finalRegulons" "$alias" "$qvalct" "$GOfile" "$netCent" "$allMarkers" "${dataset_id}_dfForRanking.txt"
    """  
}


process makeBorda {
    publishDir regOutput, mode: 'copy'
    echo true
    
    input:
    path script_makeborda
    tuple val(dataset_id), path(regulonsDF), val(metrics), val(ref)

    output:
    tuple val("${dataset_id}"), path("${dataset_id}_rankedRegulons.xlsx")
     
    """
    OMP_NUM_THREADS=1 python3 $script_makeborda "$regulonsDF" "$metrics" "${dataset_id}_rankedRegulons.xlsx" "$ref"
    """
}


process scoreEdges {
    publishDir regOutput, mode: 'copy'
    
    input:
    path script_scoreEdges
    tuple val(dataset_id), path(regulons), path(rankedRegulons), path(modules)
    
    output:
    tuple val(dataset_id), path("${dataset_id}_edgeTable.tsv")
     
    """
    OMP_NUM_THREADS=1 python3 $script_scoreEdges $regulons $rankedRegulons $modules ${dataset_id}_edgeTable.tsv
    """
}


process heatmap_tops {
    publishDir figs, mode: 'copy'
    
    input:
    path script_heatmapTops
    tuple val(dataset_id), path(rankedRegulons)
    val topRegs

    output:
    tuple val(dataset_id), path("${dataset_id}_heatmapSpecificity.svg"), path("${dataset_id}_heatmapDEcalls.svg")
     
    """
    OMP_NUM_THREADS=1 python3 $script_heatmapTops "$rankedRegulons" "${dataset_id}_heatmapSpecificity.svg" "${dataset_id}_heatmapDEcalls.svg" "$topRegs"  
    """
}


process regmaps {
    publishDir figs, mode: 'copy'
    
    input:
    path script_regmaps
    tuple val(dataset_id), path(expMatrix), path(cellClusters), path(identClust), path(rankedRegulons)
    val topRegs

    output:
    tuple val(dataset_id), path("${dataset_id}_regmap_*.svg")
     
    """
    OMP_NUM_THREADS=1 python3 $script_regmaps -c $cellClusters \
                                              -i $identClust \
                                              -r $rankedRegulons \
                                              -m $expMatrix \
                                              -t 10,25,50,100,$topRegs \
                                              -d $dataset_id
    """
}


workflow {
    check_input_files(params.script_checkInput,
        Channel.fromPath(params.expressionMatrix).collect(),
        Channel.fromPath(params.markersOut).collect(),
        Channel.fromPath(params.cell2clusters).collect(),
        Channel.fromPath(params.cluster2ident).collect(),
        Channel.fromPath(params.TF_list).collect(),
        params.termsOfInterest != null ? Channel.fromPath(params.termsOfInterest).collect() : "/dummy_path_tof",
        params.grnboostOut != null ? Channel.fromPath(params.grnboostOut).collect() : "/dummy_path_grnb",
        params.doMotifAnalysis == true? Channel.fromPath(params.featureFile_motifs) : "/dummy_path_motif",
        Channel.fromPath(params.infoTF).collect(),
        params.GOfile != null ? Channel.fromPath(params.GOfile).collect() : "/dummy_path_go",
        Channel.fromPath(params.alias).collect())

    matrix_ch = Channel.fromPath(params.expressionMatrix).map { n -> [ n.baseName.split("_")[0], n ] }
    matrix_ch.view()
    
    if (params.grnboostOut == null){
        run_grnboost(params.script_grnboost,params.TF_list,matrix_ch)
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
        unzip_motifMappings(params.featureFile_motifs)
        run_enricher_motifs(params.script_enricher,unzip_motifMappings.out,grnboost_combined_ch)
        
        filter_motifs(params.script_motifs,params.infoTF,params.motifFilter,run_enricher_motifs.out)
        filter_motifs_ch = filter_motifs.out
    } else {
        filter_motifs_dummy(grnboost_ch)
        filter_motifs_ch = filter_motifs_dummy.out
    }

    deg_ch = Channel.fromPath(params.markersOut).map { n -> [ n.baseName.split("_")[0], n ] }
    cluster_enrich_ch = get_topDEGs(params.script_topDEGs,params.tops,deg_ch)

    cluster_enrich_combined_ch = cluster_enrich_ch.join(filter_motifs_ch).join(expressed_genes_ch)
    run_enricher_cluster(params.script_enricher,cluster_enrich_combined_ch)  

    cluster_ch = Channel.fromPath(params.cell2clusters).map { n -> [ n.baseName.split("_")[0], n ] }
    filter_combined_ch = matrix_ch.join(cluster_ch).join(run_enricher_cluster.out)
    filter_expression(params.script_expTFs,params.infoTF,params.expressionFilter,filter_combined_ch)
    
    cluster_ids_ch = Channel.fromPath(params.cluster2ident).map { n -> [ n.baseName.split("_")[0], n ] }
    info_ch = matrix_ch.join(grnboost_ch).join(filter_motifs_ch).join(filter_expression.out).join(cluster_ch).join(cluster_ids_ch)
    make_info_file(params.script_info,info_ch,params.infoTF)    
    
    regulons_ident_ch = cluster_ids_ch.join(filter_expression.out)
    clustermap_regs(params.script_clustermap,regulons_ident_ch)
    
    network_centrality(params.script_networkCentrality,filter_expression.out)

    if ( params.GOfile != null ){

        getFiles_enrichment_input_ch = filter_expression.out.combine(Channel.fromPath(params.GOfile))
        getFiles_enrichment_input_ch.view()
        getFiles_enrichment_ch = getFiles_enrichment(params.script_filesEnrichment,getFiles_enrichment_input_ch)

        getFiles_enrichment_combined_ch = getFiles_enrichment_ch.join(expressed_genes_ch)
        GO_enricher(params.script_enricher,getFiles_enrichment_combined_ch)   
    }

    if ( params.GOfile != null && params.termsOfInterest != null ){        

        check_reference(params.script_checkReference,getFiles_enrichment_input_ch,params.termsOfInterest)
        check_reference_trimmed = check_reference.out.map { n,it -> [ n, it.trim() ] }
    
        rankingRef_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(network_centrality.out).join(GO_enricher.out).join(deg_ch).combine(Channel.fromPath(params.GOfile))
        ranking_df_ref(params.script_makedfRef,params.alias,params.termsOfInterest,rankingRef_combined_ch)
    
        metrics_ch = Channel.value('qval_cluster,out-degree,betweenness,closeness,GO_enrich_qval')
        metrics_combi_ch = ranking_df_ref.out.combine(metrics_ch).join(check_reference_trimmed)
    
    } else {
        rankingStd_combined_ch = cluster_ids_ch.join(filter_expression.out).join(run_enricher_cluster.out).join(network_centrality.out).join(deg_ch).combine(Channel.value(false))
        ranking_df_std(params.script_makedfStd,params.alias,rankingStd_combined_ch)    
        
        metrics_ch = Channel.value('qval_cluster,out-degree,betweenness,closeness')
        std_ch = Channel.value('std')
        metrics_combi_ch = ranking_df_std.out.combine(metrics_ch).combine(std_ch)        
    }

    makeBorda(params.script_makeborda,metrics_combi_ch)
    
    scoreEdges_ch = filter_expression.out.join(makeBorda.out).join(grnboost_ch)
    scoreEdges(params.script_scoreEdges, scoreEdges_ch)

    heatmap_tops(params.script_heatmapTops,makeBorda.out,params.topRegs)

    regmaps_ch = matrix_ch.join(cluster_ch).join(cluster_ids_ch).join(makeBorda.out)    
    regmaps(params.script_regmaps, regmaps_ch, params.topRegs)
}

workflow.onComplete {
    log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
