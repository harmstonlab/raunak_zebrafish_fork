#################

# DE-functions

##################

# Description: Functions used in 02_DE for a typical differential expression analysis. 
# Author: Qian Hui TAN
# Date last modified: 18-Apr-2023


## This function gets a list of lfcShrink fold changes from a wald dds
## It takes in an nbinomWaldTest object, a character string contrast, and an ensembl.genes object containing gene biotype for annotations. 

get_dds_res <- function(wald_dds, contrast, ensembl.genes, 
                        lfcshrinktype = "apeglm", 
                        shrink = TRUE, parallel = TRUE) {
  
  # Run: get_dds_res(wald_dds, 
  #                  contrast = c("condition", "mdd", "ctrl"),
  #                  ensembl.genes = ensembl.genes,
  #                  shrink = TRUE)
  
  # Get the results
  results_c1_control = results(wald_dds, 
                               contrast = contrast,  
                               filter = rowMeans(counts(wald_dds, 
                                                        normalized = TRUE)), 
                               test = "Wald", alpha = 0.1, 
                               independentFiltering = TRUE)
  
  # Make the condition_mdd_vs_ctrl string
  deseq_coef <- paste(contrast[1], contrast[2], "vs", 
                      contrast[3], sep = "_")
  print(deseq_coef)
  
  # Shrink
  if (shrink == TRUE){
    results_c1_control <- lfcShrink(wald_dds, 
                                    coef = deseq_coef, 
                                    res = results_c1_control, 
                                    type = lfcshrinktype, parallel = TRUE)
  }
  
  # Add gene annotations
  results_c1_control$gene_biotype = ensembl.genes$gene_biotype[match(row.names(results_c1_control), ensembl.genes$gene_id)]
  results_c1_control$external_gene_name = ensembl.genes$external_gene_name[match(row.names(results_c1_control), ensembl.genes$gene_id)]
  
  return(results_c1_control)
}




## Takes a gsea_collection, a results_table from nbinomWaldTest, and creates a GSEA plot. 
# For gsea_collection options: see http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
# Use warning: FALSE in the code chunk label to avoid repeated warning messages for tied terms. 

# Run: 
#res_pc <- results_table[rownames(results_table) %in% 
#           ensembl.genes$gene_id[ensembl.genes$gene_biotype == "protein_coding"], ]
# Run GSEA
# gsea_res_hallmark <- plotGSEA_hs(gsea_collection = "H", 
#                   results_table = res_pc, 
#            ensembl.genes = ensembl.genes, collapse = TRUE)

plotGSEA_hs <- function(gsea_collection, results_table, 
                        ensembl.genes, 
                        top_n = 10, collapse = FALSE) {
  
  collection = gsea_collection
  subcat = NULL
  
  message("Running GSEA for organism = Homo sapiens")
  
  if (gsea_collection == "reactome") {
    
    # Get entrez IDs
    entrez_ids <- as.character(ensembl.genes$entrezgene_id[match(rownames(results_table), ensembl.genes$gene_id)])
    
    
    # Get reactome pathways for given entrez IDs
    fgsea_sets <- reactomePathways(entrez_ids)

    # Create vector of fold changes, rank by lfc
    df_ranked <- results_table$log2FoldChange

    names(df_ranked) <- entrez_ids
    df_ranked <- df_ranked[order(df_ranked, decreasing = T)]
    
  } else {
    message("Obtaining gene sets from msigdbr")
    # Get the gene sets
    genesets = msigdbr(species = "Homo sapiens", 
                       category = collection, 
                       subcategory = subcat)
    
    # Split into lists
    fgsea_sets = genesets %>% 
      split(x = .$ensembl_gene, f = .$gs_name) 
    
    # Create vector of fold changes, rank by lfc
    df_ranked <- results_table$log2FoldChange
    names(df_ranked) <- rownames(results_table)
    df_ranked <- df_ranked[order(df_ranked, decreasing = T)]
  }
  
  ### -- Run GSEA -- ###
  
  message("Running GSEA")
  # Run gsea in parallel - use all but 1 core
  fgseaRes <- fgseaMultilevel(pathways = fgsea_sets, 
                              stats = df_ranked,
                              minSize = 1,
                              maxSize = 500,
                              eps = 0, 
                              BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores() - 1))
  
  #df_fgseaRes <- fgseaRes
  #df_fgseaRes <- fgseaRes %>% 
  #  dplyr::select(-leadingEdge, -ES) %>% 
  #  arrange(desc(NES), padj)

  if (collapse == TRUE) {
    # Collapse pathways
    # something wrong here- to fix tmr. 
    message("Collapsing pathways")
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.1], 
                                          fgsea_sets, df_ranked)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways, ]
    
    # Pick top n pathways
    topPathwaysUp <- mainPathways[ES > 0][head(order(pval), n = top_n), pathway]
    topPathwaysDown <- mainPathways[ES < 0][head(order(pval), n = top_n), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    # Plot 
    print(paste0(nrow(mainPathways), " collapsed pathways found.", 
                 "Plotting top ", top_n, " and bottom ", top_n, " pathways"))
    
    print(plotGseaTable(pathways = fgsea_sets[topPathways],
                        stats = df_ranked, 
                        fgseaRes = fgseaRes, 
                        pathwayLabelStyle = list(size = 10),
                        valueStyle = list(size = 10),
                        headerLabelStyle = list(fontface = "bold"),
                        gseaParam = 0.4))
    
    #df_mainpathways <- mainPathways
    #df_mainpathways <- mainPathways %>% 
    #  dplyr::select(-leadingEdge, -ES) %>% 
    #  arrange(desc(NES), padj)
    
    #return(df_mainpathways)
    return(mainPathways)
    
  } else {
    
    # Pick top n pathways
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = top_n), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = top_n), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    
    print(paste0(nrow(fgseaRes), " GSEA terms found. ", 
                 "Plotting top ", top_n, " and bottom ", top_n, " pathways"))
    
    print(plotGseaTable(pathways = fgsea_sets[topPathways],
                        stats = df_ranked, 
                        fgseaRes = fgseaRes, 
                        pathwayLabelStyle = list(size = 10),
                        valueStyle = list(size = 10),
                        headerLabelStyle = list(fontface = "bold"),
                        gseaParam = 0.4))
    
    
    return(fgseaRes)
    #return(df_fgseaRes)
  }
  
  
}

make_gsea_readable_hs <- function(fgseaRes, id_type){
  
  if(id_type == "ENTREZ"){
    # Map IDs
    fgseaRes <- fgseaRes[, leadingEdge := mapIdsList(
      x = org.Hs.eg.db, 
      keys = leadingEdge,
      keytype = "ENTREZID", 
      column = "SYMBOL")]
  } else if (id_type == "ENSEMBL") {
    # Map IDs
    fgseaRes <- fgseaRes[, leadingEdge := mapIdsList(
      x = org.Hs.eg.db, 
      keys = leadingEdge,
      keytype = "ENSEMBL", 
      column = "SYMBOL")]
  } else {
    stop("id_type must be either ENTREZ or ENSEMBL")
  }

  
  # Convert list to string
  fgsea_export <- as.data.frame(fgseaRes) %>% 
    mutate(leadingEdge = sapply(leadingEdge, toString)) 
  
  
  return(fgsea_export)
}



custom_gsea_msigdb <- function(interesting_pathways, 
                               msigdb_collection, 
                               fgseaRes, results_table,
                                 font_size = 10,
                                 gsea_param = 0.4){
  
  collection = msigdb_collection
  subcat = NULL
  
  message("Obtaining gene sets from msigdbr")
  # Get the gene sets
  genesets = msigdbr(species = "Homo sapiens", 
                     category = collection, 
                     subcategory = subcat)
  
  # Split into lists
  fgsea_sets = genesets %>% 
    split(x = .$ensembl_gene, f = .$gs_name) 
  
  # Create vector of fold changes, rank by lfc
  df_ranked <- results_table$log2FoldChange
  names(df_ranked) <- rownames(results_table)
  df_ranked <- df_ranked[order(df_ranked, decreasing = T)]
  
  # Arrange pathways
  
  fgsea_interesting <- fgseaRes[pathway %in% interesting_pathways, ]
  
  pathways_up <- fgsea_interesting[ES > 0][order(pval), pathway]
  pathways_down<- fgsea_interesting[ES < 0][order(pval), pathway]
  pathways_combined <- c(pathways_up, rev(pathways_down))
  
  
  #pathways_up <- fgseaRes[ES > 0][order(pval), pathway]
  #pathways_down<- fgseaRes[ES < 0][order(pval), pathway]
  #pathways_combined <- c(pathways_up, pathways_down)
  
  
  
  print(plotGseaTable(pathways = fgsea_sets[pathways_combined],
                      stats = df_ranked, 
                      fgseaRes = fgseaRes, 
                      pathwayLabelStyle = list(size = font_size),
                      valueStyle = list(size = font_size),
                      headerLabelStyle = list(fontface = "bold",
                                              size = font_size),
                      gseaParam = gsea_param))
}



custom_gsea_reactome <- function(interesting_pathways, fgseaRes, results_table,
                                 font_size = 10,
                                 gsea_param = 0.4){
  

  
  message("Obtaining gene sets from msigdbr for reactome")
  # Get entrez IDs
  entrez_ids <- as.character(ensembl.genes$entrezgene_id[match(rownames(results_table), ensembl.genes$gene_id)])
  # Get reactome pathways for given entrez IDs
  fgsea_sets <- reactomePathways(entrez_ids)
  
  
  # Create vector of fold changes, rank by lfc
  df_ranked <- results_table$log2FoldChange
  #names(df_ranked) <- rownames(results_table)
  names(df_ranked) <- as.character(ensembl.genes$entrezgene_id[match(rownames(results_table), ensembl.genes$gene_id)])
  df_ranked <- df_ranked[order(df_ranked, decreasing = T)]
  
  # Arrange pathways
  
  fgsea_interesting <- fgseaRes[pathway %in% interesting_pathways, ]
  
  pathways_up <- fgsea_interesting[ES > 0][order(pval), pathway]
  pathways_down<- fgsea_interesting[ES < 0][order(pval), pathway]
  pathways_combined <- c(pathways_up, rev(pathways_down))
  
  
  #pathways_up <- fgseaRes[ES > 0][order(pval), pathway]
  #pathways_down<- fgseaRes[ES < 0][order(pval), pathway]
  #pathways_combined <- c(pathways_up, pathways_down)
  

  
  print(plotGseaTable(pathways = fgsea_sets[pathways_combined],
                      stats = df_ranked, 
                      fgseaRes = fgseaRes, 
                      pathwayLabelStyle = list(size = font_size),
                      valueStyle = list(size = font_size),
                      headerLabelStyle = list(fontface = "bold",
                                              size = font_size),
                      gseaParam = gsea_param))
}














