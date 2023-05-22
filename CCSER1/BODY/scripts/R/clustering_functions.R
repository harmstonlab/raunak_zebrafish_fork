
## make_svaplot ## 
## Run: make_svaplot(rld_sva, var = "SV1")

make_svaplot <- function(rld_sva, var){
  boxplot(rld_sva[[var]] ~ rld_sva$condition, 
          xlab="sample", ylab = var)
  
  stripchart(rld_sva[[var]] ~ rld_sva$condition,  
             method = "jitter",       
             pch = 19,          
             col = 4,            
             vertical = TRUE,      
             add = TRUE)
  abline(h = 0, col = "red", lty = 2)
}

make_pca <- function(rld, intgroup,  
                     title = "title",         
                     xlimits = c(-70, 70),
                     ylimits = c(-25, 25),
                     label = FALSE){
  # Calculations
  ntop = 500
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  data <- plotPCA(rld, intgroup = intgroup, 
                  returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"), digits = 2)
  
  # Plot 
  pca_plot <- ggplot(as.data.frame(data), aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = as.data.frame(data)[[intgroup]]), 
               size = 2, alpha = 0.8) +
    labs(title = title,
         subtitle = paste0("By ", intgroup), 
         colour = intgroup) +
    
    # Add scale annotations
    scale_x_continuous(paste0("PC1: ",percentVar[1],"% variance"), 
                       limits = xlimits) +
    scale_y_continuous(paste0("PC2: ",percentVar[2],"% variance"), 
                       limits = ylimits) +
    scale_color_brewer(palette = "Paired") + 
    
    # Make 1 unit on x-axis equal to 1 unit on y-axis
    coord_fixed(ratio = 1) +
    theme_classic()
  
  if(label == TRUE){
    pca_plot <- pca_plot + 
      geom_text_repel(data = data, aes(PC1,PC2, label = name), 
                      hjust = 0.5, box.padding = 0.5, size = 3,
                      max.overlaps = Inf)
    return(pca_plot)
  } else {
    return(pca_plot)
  }
  
}



plot_gene <- function(gene_of_interest, title = "title"){
  
  # Subset
  df_nc <- counts(dds, normalized = TRUE)[gene_of_interest, ] %>% melt()
  df_nc$sample <- rownames(df_nc)
  
  df_nc$condition <- colData(dds)$condition[match(df_nc$sample, colData(dds)$sample_id)]
  
  # Plot
  df_plot <- ggplot(df_nc, aes(x = condition, y = value)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = title,
         subtitle = gene_of_interest, 
         x = "", y = "Normalized counts")
  
  return(df_plot)
  
}

plot_gene_bysample <- function(gene_of_interest, title = "title",
                               sample_order){
  
  # Subset
  df_nc <- counts(dds, normalized = TRUE)[gene_of_interest, ] %>% melt()
  df_nc$sample <- rownames(df_nc)
  
  df_nc$sample <- factor(df_nc$sample, levels = sample_order)
  df_nc$condition <- colData(dds)$condition[match(df_nc$sample, colData(dds)$original_label)]
  
  # Plot
  df_plot <- ggplot(df_nc, aes(x = sample, y = value,
                               fill = condition)) +
    geom_col() +
    scale_fill_brewer(palette = "Paired") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = title,
         subtitle = gene_of_interest, 
         x = "", y = "Normalized counts") +
    theme(legend.position = "none")
  
  return(df_plot)
  
}

# Run: zscore_boxplot(kmeans_cl_k6, clust_num = 1, sample_order)

zscore_boxplot <- function(kmeans_cl, clust_num, 
                           sample_order = sample_order){
  # Get the genes
  c3_genes <- names((kmeans_cl)[[clust_num]])
  print(paste0(length(c3_genes), " genes in cluster ", clust_num))
  
  # Plot the normalized counts
  c3_all <- rld_z[c3_genes, ] %>% 
    melt()
  colnames(c3_all) <- c("gene_id", "sample", "vst_z")
  # Reorder samples
  sample_order
  c3_all$sample <- factor(c3_all$sample, levels = sample_order)
  
  # Plot the zscores
  zscore_plot <- ggplot(c3_all, aes(x = sample, y = vst_z)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_boxplot() +
    scale_y_continuous(limits = c(-4, 4)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste0("Zscore boxplot, cluster ", clust_num),
         x = "", 
         y = "z-score") 
  
  return(zscore_plot)
  
}


custom_ego <- function(interesting_pathways, 
                       ego_tibble,
                       title = "title", 
                       font_size = 14) {
  
  # Get relevant columns
  c2_interesting_egos <- ego_tibble%>% 
    dplyr::select("Description", "GeneRatio", "p.adjust", "Count") %>% 
    filter(Description %in% interesting_pathways) %>% 
    # Convert GeneRatio from fraction to decimal
    mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>% 
    arrange(GeneRatio) %>% 
    # Preserve sorted order
    mutate(Description = factor(Description, levels = unique(Description)))
  
  
  # Plot
  df_plot <- ggplot(c2_interesting_egos, aes(x = GeneRatio, y = Description,
                                             color = p.adjust)) +
    geom_point(aes(size = Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = title,
      y = "") +
    theme_light() +
    theme(axis.text.y = element_text(size = font_size))
  
  return(df_plot)
  
}



cluster_kmeans <- function(rld_z, nclust, plot_sil = FALSE){
  
  ## Perform k means clustering
  set.seed(1)
  nclust = nclust
  kmeans_coef = kmeans(rld_z, nclust, nstart = 1000, iter.max = 50)
  
  ## Add silhouette plots for diagnosis 
  if(plot_sil == TRUE){
    # Set colors
    c8 <- c("#EB6424","#FA9500", "#FFBB00",
                     "#FFEE00", "#CCBA00",
                     "#B99900", "#AA9900", "#AAAA99") 
                     
    # Calculate distance
    d <- dist(rld_z)
    
    # Plot
    plot(cluster::silhouette(kmeans_coef$cluster, d),
         col= c8[1:nclust],
         border = NA, 
         main = paste("Silhouette plot, k = ", nclust)
    )
    
    return(kmeans_coef)
    
  } else {
    return(kmeans_coef)
  }
}



plot_kmeans_heatmap <- function(rld_z, k_coef,
                                sample_order){
  
  # Arrange samples in correct order
  #order_samples <- as.factor(sample_order)
  
  #rld_z <- rld_z[ ,order_samples]
  
  # Set up heatmap
  # Thresholding it to 3 standard deviations away from the median. 
  thr = 3
  rld_z[rld_z > thr] = thr
  rld_z[rld_z < -thr] = -thr
  
  # Sort out color scheme
  paletteLength = 20 
  breaksList <- c(seq(-thr, 0, length.out = ceiling(paletteLength/2) + 1), 
                  seq(thr/paletteLength, thr, length.out = floor(paletteLength/2)))
  
  color = c(colorRampPalette(c("mediumblue", "white"))(14), 
            colorRampPalette(c("white", "firebrick2"))(14))
  breaksList = seq(-3, 3, length.out = 29)
  
  cs = k_coef$cluster
  #cs <- factor(cs, levels = c(2, 1))
  
  # order things in correct order
  z.toplot = rld_z[order(cs), sample_order]
  
  
  # Heatmap according to clusters
  heatmap_cl = pheatmap::pheatmap(z.toplot, 
                                  breaks = breaksList, 
                                  cluster_col = FALSE,
                                  cluster_rows = FALSE, 
                                  show_rownames = FALSE,
                                  show_colnames = TRUE, 
                                  #color = color,
                                  #annotation_col = annotation,
                                  #annotation_colours = anno_colours,
                                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(29),
                                  fontsize_col = 12, 
                                  legend = TRUE,
                                  border_color = NA, 
                                  main = paste("Heatmap"),
                                  angle_col = 45
  )
  
  heatmap_cl
  
}


get_cluster_genes <- function(k_coef, nclust){
  
  # Loop to obtain genes for each cluster
  kmeans_cl = c()
  
  for (i in 1:nclust){
    kmeans_cl[[i]] = k_coef$cluster[
      which(k_coef$cluster == i)
    ]
  }
  
  return(kmeans_cl)
}



plotEGO_dr = function(clust_target_genes,
                   universe = universe, 
                   ont = "BP", 
                   title = "title", 
                   subtitle = NULL, 
                   font_size = 14){
  
  # Run KEGG enrichment
  message("Running GO for organism = zebrafish")
  
  cl_target_ego = enrichGO(gene = clust_target_genes, 
                           universe = universe,
                           OrgDb = org.Dr.eg.db,
                           keyType = 'ENSEMBL', 
                           ont = ont, 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.1,
                           qvalueCutoff = 0.1,
                           readable = TRUE)
  
  # If no GO terms found, return warning message and a tibble with NA
  if(nrow(data.frame(cl_target_ego)) == 0) {
    warning(paste0("No GO enrichment found. Returning a NA tibble."
    ))
    return(tibble(`Description` = "NA"))
  } else {
    print(clusterProfiler::dotplot(cl_target_ego,
                  title = title, 
                  font.size = font_size))
    # Print number of enrichments found
    print(paste0(nrow(cl_target_ego), " enrichments found"))
    
    return(as_tibble(cl_target_ego))
  }
  
  
  
}

# Searches a given tibble for a term. Basically grep but with fewer characters to type. 

findEGO <- function(ego_tibble, term, print_top_matches = FALSE){
  # Grep search for given term 
  matches <- ego_tibble[grep(term, ego_tibble$Description, ignore.case = TRUE), ]
  
  # Output number of matches
  print(paste0(nrow(matches), " matches found."))
  
  print(matches$Description)
  
  # If TRUE, print top 3 matches along with the genes
  
  if(print_top_matches == TRUE) {
    # If < 4 matches found, print all of them
    if (nrow(matches) < 4) {
      print(paste0("Printing all ", nrow(matches), " matches."))
    
        as.data.frame(matches[ , colnames(matches) %in% c("Description", "GeneRatio", "p.adjust", "geneID")])
      } else {
    
       print(paste0("Matches are: "))
        print(matches$Description)
    
    print(paste0("Printing top 3 matches."))
    
    df_match <- as.data.frame(matches[ , colnames(matches) %in% c("Description", "GeneRatio", "p.adjust", "geneID")])
    
    df_match <- df_match[order(df_match$p.adjust, decreasing = FALSE), ]
    
    head(df_match, n = 3)
    }
  }
 
}


get_entrez <-function(ensembl_ids, ensembl_genes){
  entrez <- na.omit(ensembl.genes[ensembl.genes$gene_id %in% ensembl_ids, ]$entrezgene_id)
  return(entrez)
}

# Takes in an ego tibble, subsets interesting pathways. 
custom_ego_table <- function(ego_tibble, interesting_pathways) {
  ego_tbl <- ego_tibble[ego_tibble$Description %in% interesting_pathways, 
                        colnames(ego_tibble) %in% c("Description", "GeneRatio", "p.adjust", "geneID")]
  return(ego_tbl)
}



### This function ONLY WORKS FOR HUMAN GENES

# Run: 
## Get the entrez IDs for cluster 2
# c1_entrez <- na.omit(ensembl.genes[ensembl.genes$gene_id %in% names((kmeans_cl_k2)[[1]]), ]$entrezgene_id)
# k2_c1_kegg <- plotKEGG(c1_entrez, title = "KEGG, cluster 1")

plotKEGG_dr <- function(target_genes_entrez, title = "title", font_size = 14){
  
  # Run KEGG enrichment
  message("Running KEGG for organism = zebrafish")
  ekegg = clusterProfiler::enrichKEGG(target_genes_entrez,
                                      organism = "dre",
                                      keyType = "ncbi-geneid")
  
  # If no GO terms found, return warning message and a tibble with NA
  if(nrow(data.frame(ekegg)) == 0) {
    warning(paste0("No KEGG enrichment found. Returning a NA tibble."
    ))
    return(tibble(`Description` = "NA"))
  } else {
    
    # Print the kegg dotplot
    print(dotplot(ekegg, title = title, 
                  font.size = font_size))
    
    # Convert entrez to gene symbols
    ekegg = DOSE::setReadable(ekegg, OrgDb = "org.Dr.eg.db", keyType = "ENTREZID")
    
    # Print number of enrichments found
    print(paste0(nrow(ekegg), " enrichments found"))
    
    # Return the result as a data frame
    return(as.data.frame(ekegg))
  }
}



custom_ekegg <- function(interesting_pathways, 
                       ekegg_tibble,
                       title = "title", 
                       font_size = 14) {
  
  # Get relevant columns
  interesting_ekeggs  <- ekegg_tibble%>% 
    dplyr::select("Description", "GeneRatio", "p.adjust", "Count") %>% 
    filter(Description %in% interesting_pathways) %>% 
    # Convert GeneRatio from fraction to decimal
    mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>% 
    arrange(GeneRatio) %>% 
    # Preserve sorted order
    mutate(Description = factor(Description, levels = unique(Description)))
  
  
  # Plot
  df_plot <- ggplot(interesting_ekeggs, aes(x = GeneRatio, y = Description,
                                             color = p.adjust)) +
    geom_point(aes(size = Count)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = title,
      y = "") +
    theme_light() +
    theme(axis.text.y = element_text(size = font_size))
  
  return(df_plot)
  
}






# Plots a nice heatmap
# Run: make_heatmap(rld_z_plot, names(kmeans_cl[[1]]), title = "Cluster 1")

make_heatmap <- function(count_matrix, 
                         genes,
                         show_row_names = TRUE,
                         rownames_fontsize = 12,
                         title = "title"){
  
  
  # Heatmap according to clusters
  heatmap_cl = ComplexHeatmap::Heatmap(
    #matrix = rld_z_plot[names(kmeans_cl[[1]]), ], 
    matrix = count_matrix[genes, ],
    col = colorRamp2(c(min(count_matrix), 0, max(count_matrix)), 
                     c("mediumblue", "white", 
                                   "firebrick2")
    ),
    column_title = title,
    column_title_side = "top",
    name = "z-score",
    
    # Graphic parameters: rows
    cluster_rows = TRUE,
    row_names_gp =  gpar(fontsize = rownames_fontsize),
    row_names_side = "right",
    #row_dend_gp = gpar(fontsize = 12),  
    #row_dend_width = unit(20, "mm"),
    
    clustering_distance_rows = "euclidean", # clust by genes   
    
    # Graphic parameters: columns
    column_names_gp = gpar(fontsize = 14), 
    column_names_side = "bottom",
    column_dend_gp = gpar(fontsize = 12),
    column_dend_side = "top",
    #column_dend_height = unit(10, "mm"),
    
    show_column_names = TRUE, 
    column_names_rot = 45,
    show_row_names = show_row_names,
    
    # column_order = c("control_1", "control_2", 
    #                 "damaged_1", "damaged_2")
  )
  
  heatmap_cl
  
}

### ------ clusterHeatmap ------ ###

clusterHeatmap = function(rld_z,
                          kmeans_cl, 
                          clust_num, 
                          sample_order,
                          cluster_rows = TRUE,
                          cluster_columns = FALSE, 
                          show_row_names = FALSE,
                          show_column_names = TRUE,
                          font_size = 14
) {
  
  # Set up heatmap
  # Thresholding it to 3 standard deviations away from the median. 
  thr = 3
  rld_z[rld_z > thr] = thr
  rld_z[rld_z < -thr] = -thr
  
  # Subset rld   
  rld_z_cl = rld_z[which(rownames(rld_z) %in% 
                           names(kmeans_cl[[ clust_num ]])), ]

  # Create a heatmap
  
  heatmap_cl = ComplexHeatmap::Heatmap(
    matrix = rld_z_cl, 
    col = colorRamp2(c(-3, 0, 3), 
                     c("mediumblue", "white", 
                                   "firebrick2")
    ),
    column_title = paste0("Heatmap for cluster ", clust_num),
    column_title_side = "top",
    name = "z-score",
    ## Top annotation (couldn't get this to display properly)
    #top_annotation = ha,
    
    ## -- Cluster rows --- ## 
    cluster_rows = cluster_rows,
    clustering_distance_rows = "euclidean", # clust by genes 
    
    ## Rows: graphic parameters  
    
    show_row_names = show_row_names,
    row_names_gp =  gpar(fontsize = 11),
    row_names_side = "right",
    row_dend_gp = gpar(fontsize = 12),  
    row_dend_width = unit(20, "mm"),
    
    ## -- Cluster columns -- ##
    cluster_columns = cluster_columns,
    
    ## Columns: graphic parameters
    show_column_names = show_column_names, 
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 14), 
    column_names_side = "bottom",
    column_dend_gp = gpar(fontsize = 12),
    column_dend_side = "top",
    column_dend_height = unit(10, "mm"),
    
    column_order = factor(sample_order, levels = sample_order)
  )
  heatmap_cl
  
}