#######################################
## Differential Expression Functions ##
#######################################

## This function is to write out the files of differentially expressed genes
## It takes in a DESeqResults object, and two strings -- the numerator and denominator used in the analysis -- and writes out csv files

write_files <- function(results, numerator, denominator, output_directory = output_dir){
  # these are all the genes that are differentially expressed between the two conditions, not just the significant ones
  write.csv(results, paste0(output_directory,numerator,"_",denominator,"_all.csv"), row.names = TRUE, col.names = TRUE, sep = " ")
  
  # these are the genes that are significantly differentially expressed by FDR 10% and abs(log2fc) > log2(1.5)
  sig_padj_genes <- results[!is.na(results$padj),]
  sig_padj_genes <- sig_padj_genes[sig_padj_genes$padj < 0.1,]
  sig_padj_fc_genes <- sig_padj_genes[abs(sig_padj_genes$log2FoldChange) > lfc.threshold,]
  write.csv(sig_padj_fc_genes, paste0(output_dir,numerator,"_",denominator,"_significant.csv"), row.names = TRUE, col.names = TRUE)
}

## This function plots the volcano plot
## It takes in a data frame and two strings which are used for the title of the plot

generate_volcano <- function(data_frame, numerator, denominator){
  lfc.threshold = log2(1.5)
  tmp = as.data.frame(data_frame)
  tmp$signif = ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.01, "U1", 
                      ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.05, "U2",
                             ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.1, "U3",
                                    ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.01, "D1", 
                                           ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.05, "D2",
                                                  ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.1, "D3",                                                  "N"))))))
  tmp$signif = factor(tmp$signif, c("N", "U1", "U2", "U3", "D3", "D2", "D1"))
  
  x = ggplot(data=tmp, aes(x=log2FoldChange, y=-log10(padj), colour= signif)) + geom_point(alpha = 1,  size = 1) +
    ggtitle(paste("Volcano Plot:", numerator, "vs.", denominator)) + scale_x_continuous("log2(fold change)", limits=c(-10, 10)) +    
    scale_y_continuous("-log10(FDR)") + geom_vline(xintercept = lfc.threshold, linetype="dotdash") +
    geom_vline(xintercept = -1*(lfc.threshold), linetype="dotdash") +
    geom_hline(yintercept = -log10(0.1), colour="gray40", linetype="dotdash") +   
    geom_hline(yintercept = -log10(0.05), colour="gray40", linetype="dotdash") + 
    geom_hline(yintercept = -log10(0.01), colour="gray40", linetype="dotdash") + 
    scale_colour_manual("", values=c("#666666", "#d73027", "#f46d43", "#fdae61", "#abd9e9", "#74add1", "#4575b4" ), labels = c("N", "U1", "U2", "U3", "D3", "D2", "D1")) + theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 20), axis.title=element_text(size=16,face="bold"))
  return(x)
  print(table(tmp$signif))
}


## This function is similar to the above, but changes the alpha values to highlight genes of interest. 
## Genes of interest must be rownames(data_frame) - usually the ensembl ID (ENSG.....)

highlight_volcano <- function(data_frame, numerator, denominator, 
                              genes_of_interest,
                              subtitle = ""){
  
  lfc.threshold = log2(1.5)
  tmp = as.data.frame(data_frame)
  
  tmp$signif = ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.01, "U1", 
                      ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.05, "U2",
                             ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.1, "U3",
                                    ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.01, "D1", 
                                           ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.05, "D2",
                                                  ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.1, "D3",                                                  "N"))))))
  tmp$signif = factor(tmp$signif, c("N", "U1", "U2", "U3", "D3", "D2", "D1"))
  
  # Change alpha for genes of interest
  tmp$label = ifelse(rownames(tmp) %in% genes_of_interest, TRUE, FALSE)
  alpha <- ifelse(tmp$label, 0.9, 0.01)
  
  labels <- tmp %>% filter(label == TRUE)
  
  # Plot
  x = ggplot(data = tmp, 
             aes(x=log2FoldChange, y=-log10(padj), colour = signif)) + 
    geom_point(alpha = alpha, size=2.00) +
    ggtitle(paste("Volcano Plot:", numerator, "vs.", denominator)) +
    labs(subtitle = subtitle) + 
    scale_x_continuous("log2(fold change)", limits=c(-5, 5)) +    
    scale_y_continuous("-log10(FDR)", limits = c (0, 100)) + 
    geom_vline(xintercept = lfc.threshold, linetype="dotdash") +
    geom_vline(xintercept = -1*(lfc.threshold), linetype="dotdash") +
    geom_hline(yintercept = -log10(0.1), colour="gray40", linetype="dotdash") +   
    geom_hline(yintercept = -log10(0.05), colour="gray40", linetype="dotdash") + 
    geom_hline(yintercept = -log10(0.01), colour="gray40", linetype="dotdash") + 
    scale_colour_manual("", 
                        values=c("#666666", "#d73027", "#f46d43", "#fdae61",
                                          "#abd9e9", "#74add1", "#4575b4" ), 
                                          labels = c("N", "U1", "U2", "U3", 
                                                     "D3", "D2", "D1")) + 
    theme_classic() + 
    theme(legend.position = "none", plot.title = element_text(size = 20),
          axis.title=element_text(size=16,face="bold")) +
    geom_label_repel(data = labels, aes(label = external_gene_name), 
                     max.overlaps = Inf, 
                     force = 20, box.padding = 2)
  print(x)
  
  return(labels)
  

}

## This function generates the MA plots with significant changes above the threshold coloured in red and significant changes below the threshold coloured in blue
## It takes in a DESeqResults object, uses the plotMA function from DESeq2 to obtain the necessary data frame to plot

generate_ma <- function(results){
  df <- DESeq2::plotMA(results, ylim = c(-10,10), colSig = "red", returnData = TRUE)
  plot <- df %>%
    mutate(signif = ifelse(lfc > lfc.threshold & isDE == TRUE, "U", 
                           ifelse(lfc < -lfc.threshold & isDE == TRUE, "D", "N"))) %>%
    ggplot(aes(x=mean, y=lfc, colour = signif)) + 
    geom_point(size = 1.5, alpha = 0.8) + 
    theme_classic() + 
    geom_hline(yintercept=0, colour="grey40", lwd = 1) + 
    #stat_smooth(se = FALSE, method = "loess", color = "red3") + 
    theme_classic() + 
    scale_colour_manual(values=c("#4575b4","#a3a3a3","#d73027"), labels = c("D", "N", "U")) +
    ylim(c(-10,10)) +
    theme(legend.position = "none") +
    ylab("Log fold change") +
    xlab("Mean of normalized counts") +
    scale_x_log10()
  return(plot)
}

############################################
## Visualisation and Enrichment Functions ##
############################################

## This function generates the zscore table ordered by cluster number
## It takes in the full zscore matrix and the number of clusters desired

generate_data <- function(zscores, n_clust, fun = c("kmeans","pam")){
  set.seed(2)
  if(fun == "kmeans") {
    km = kmeans(zscores, n_clust) ## Running kmeans will return 9 components and var[1] will give back the first component which is "cluster", var[2] will return "centers" etc.
    clust <- cbind(zscores, km$cluster)
  } else {
    dd = as.dist((1 - cor(t(zscores)))/2)
    pam_clust = pam(dd, n_clust, diss = TRUE)
    clust <- cbind(zscores, pam_clust$clustering)
  }
  n_col = dim(clust)[2]
  colnames(clust)[n_col] = "Cluster"
  ordered <- order(clust[,n_col])
  clust <- clust[ordered, ]
  return(clust)
}

## This function generates the heatmap using the output from the above function, annotation data, and annotation colours

generate_heatmap <- function(clust_data, annotation, annotation_col){
  return(pheatmap(clust_data[,1:(ncol(clust_data)-1)],
                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                  fontsize_row = 5.5,
                  annotation_col = annotation,
                  annotation_colors = anno_colours,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE))
}

## This function plots the barplots using results from EnrichR 
## It takes in an enrichR result, name of the plot, and the number of terms to show (default being 20)

plot_enrichr <- function(data.frame, name, showCategory = 20){
  plot = data.frame %>%
    arrange(Adjusted.P.value) %>%
    head(showCategory) %>%
    mutate(Term = gsub("\\([^()]*\\)", "", Term),
           Term = unique(Term)) %>%
    mutate(Annotated = as.numeric(str_extract(as.character(Overlap), "\\d+$")),
           Significant = as.numeric(str_extract(as.character(Overlap), "^\\d+")),
           Ratio = Significant/Annotated) %>%
    ggplot(aes(x = reorder(Term, dplyr::desc(Adjusted.P.value)), y = Ratio, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    xlab(NULL) + 
    ylab("Gene Ratio") +
    scale_fill_continuous(low = "red", high = "blue") + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(title = name,
         fill = "Adjusted p-value") +
    guides(fill = guide_colorbar(title = "Adjusted p-value", reverse = TRUE)) +
    theme(axis.text.x = element_text(size=8))
  return(plot)
}

## This function is purely for extracting the legend of a plot to be used in plot_grid
extracting_legend <- function(df, name){
  fig = df[name,] %>%
    gather() %>%
    `colnames<-`(c("Sample", "Value")) %>%
    mutate(Condition = ifelse(grepl("WT", Sample), "WT", 
                              ifelse(grepl("Wg_Ab", Sample), "Wg_Ab",
                                     ifelse(grepl("Wg", Sample, fixed = TRUE), "Wg", "Ab")))) %>%
    mutate(Condition = factor(Condition, levels = c("WT", "Ab", "Wg", "Wg_Ab"))) %>%
    ggplot(aes(x = Condition, y = Value, fill = Condition)) + 
    geom_boxplot() + 
    theme_classic() + 
    theme(axis.text.x = element_blank(), legend.position="top") +
    scale_fill_manual(values = c("#C05746", "#E9B44C", "#7EB09B", "#086788")) +
    theme(axis.text.x = element_text(angle = 90, colour="black", hjust = 1)) + 
    xlab(NULL) +
    ylab(NULL) +
    labs(title = name)
  
  legend <- get_legend(fig)
  return(legend)
}

## This function plots the boxplots for gene expression profiles
## It takes in a data frame from which to get the normalised counts, the name of the gene to plot, the lower y-axis limit and upper y-axis limit

gene_boxplot <- function(df, name, lower_y_limit = 0, upper_y_limit = NA){
  plot = df[name,] %>%
    gather() %>%
    `colnames<-`(c("Sample", "Value")) %>%
    mutate(condition = ifelse(grepl("WT", Sample), "WT", 
                              ifelse(grepl("Wg_Ab", Sample), "Wg_Ab",
                                     ifelse(grepl("Wg", Sample, fixed = TRUE), "Wg", "Ab")))) %>%
    mutate(condition = factor(condition, levels = c("WT", "Ab", "Wg", "Wg_Ab"))) %>%
    ggplot(aes(x = condition, y = Value, fill = condition)) + 
    geom_boxplot() + 
    theme_classic() + 
    scale_fill_manual(values = c("#C05746", "#E9B44C", "#7EB09B", "#086788")) +
    theme(axis.text.x = element_text(angle = 90, colour="black", hjust = 1)) + 
    theme(axis.text.x = element_blank(), legend.position="none") +
    xlab(NULL) +
    ylab(NULL) +
    labs(title = name) +
    scale_y_continuous(limits = c(lower_y_limit, upper_y_limit), breaks = seq(lower_y_limit,upper_y_limit,by = (upper_y_limit-lower_y_limit)/5))
  return(plot)
}

cluster_boxplot <- function(cluster_df, name){
  return(cluster_df %>%
           gather() %>%
           mutate(condition = ifelse(grepl("WT", variable), "WT", 
                              ifelse(grepl("Wg_Ab", variable), "Wg_Ab",
                                     ifelse(grepl("Wg", variable, fixed = TRUE), "Wg", "Ab")))) %>%
           mutate(condition = factor(condition, levels = c("WT", "Ab", "Wg", "Wg_Ab"))) %>%
           ggplot(aes(x = condition, y = value, fill = condition)) + 
           geom_boxplot() + 
           theme_classic() + 
           scale_fill_manual(values = c("#C05746", "#E9B44C", "#7EB09B", "#086788")) +
           theme(axis.text.x = element_text(angle = 90, colour="black", hjust = 1)) + 
           theme(axis.text.x = element_blank(), legend.position="none") +
           xlab(NULL) +
           ylab(NULL) +
           labs(title = name) #+
           #scale_y_continuous(limits = c(lower_y_limit, upper_y_limit))
  )
}


## This function takes the cluster of interest and returns the fasta file in the stated output directory
generate_fasta <- function(data_frame, cluster_name){
  granges <- sig_de_granges[rownames(data_frame),]
  fasta <- getSeq(Dmelanogaster, granges)
  writeXStringSet(fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}


## This function takes the cluster of interest and returns the fasta file of EVERYTHING BUT THE CLUSTER OF INTEREST in the stated output directory (differentially expressed genes)
extract_control_deg <- function(cluster_of_interest, cluster_name){
  deg_control <- sig_de_granges[-which(sig_de_granges$gene_id %in% rownames(cluster_of_interest)),] #note the order
  deg_fasta <- getSeq(Dmelanogaster, deg_control)
  writeXStringSet(deg_fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}


## This function takes the cluster of interest and returns the fasta file of EVERYTHING BUT THE CLUSTER OF INTEREST in the stated output directory (all expressed genes)
extract_control_all <- function(cluster_of_interest, cluster_name){
  all_control <- all_granges[-which(all_granges$gene_id %in% rownames(cluster_of_interest)),] #note the order
  all_fasta <- getSeq(Dmelanogaster, all_control)
  writeXStringSet(all_fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}
