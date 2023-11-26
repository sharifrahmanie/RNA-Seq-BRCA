# Package installation
# install.packages(c("BiocManager",
#                    "readr", "ggplot2", "magick"))
# BiocManager::install(c("DESeq2", 
#                        "genefilter",
#                        "org.Hs.eg.db", 
#                        "clusterProfiler",
#                        "ComplexHeatmap"))
# library loading
suppressMessages({
  require(readr)
  require(ggplot2)
  require(DESeq2)
  require(genefilter)
  require(org.Hs.eg.db)
  require(ComplexHeatmap)
  require(clusterProfiler)
})
set.seed(1234)
# Reading and preparing files
counts <- data.frame(read_tsv("TCGA-BRCA.htseq_counts_gene_name.tsv"))
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- round(2^counts -1)
length(which(is.na(counts)))

pheno <- read_tsv("TCGA-BRCA.pheno.tsv")
colnames(pheno)[1] <- "sample"
pheno <- pheno[order(pheno$sample), ]
pheno$type <- factor(pheno$type)
counts <- counts[, which(colnames(counts) %in% pheno$sample)]
counts <- counts[, order(colnames(counts))]
identical(colnames(counts), pheno$sample)

# Filtering low count genes
counts <- varFilter(as.matrix(counts))

# DESeq2 analysis
deseq <- DESeqDataSetFromMatrix(countData = counts,
                                colData = pheno,
                                design = ~ type)
deseq <- DESeq(deseq)
deseq_result <- function(deseq, 
                         FC, 
                         p_adj,
                         ntop) {
  res <- results(deseq)
  res <- data.frame(res)
  res <- cbind(rownames(res), res)
  colnames(res)[1] <- "GeneName" 
  up_down <- res[which(abs(res$log2FoldChange) > FC & res$padj < p_adj), ]
  up <- up_down[which(up_down$log2FoldChange > FC), ]
  down <- up_down[which(!up_down$log2FoldChange > FC), ]
  up <- up[order(up$log2FoldChange, decreasing = TRUE),]
  down <- down[order(down$log2FoldChange, decreasing = FALSE),]
  up_down <- rbind(up, down)
  up_down <- up_down[, c(1,3,6:7)]
  message(cat("Up-regulate: ", nrow(up), "|", "Down-rgulated: ", nrow(down)))
  message(cat("Top 10 Up-regulated", "Top 10 Down-regulated"))
  write_csv(up_down, "diff_result_EdrisSharifRahmani.csv")
  write_csv(res, "ALL_genes.csv")
  return(list(UP = head(up, 10), DOWN = head(down, 10)))
}


# Drawing volcano plot 
volcano_plot <- function(padjlevel,
                         Up_FC,
                         Down_FC,
                         ntop) {
  volcano <- data.frame(read_csv(paste0("ALL_genes.csv")))
  volcano <- na.omit(volcano)
  colnames(volcano)
  volcano <- volcano %>%
    mutate(log10padj = -log10(padj)) %>%
    mutate(DEG = "NotSignificance") %>%
    mutate(DEG = ifelse(log10padj > -log10(padjlevel) & log2FoldChange > Up_FC, "Upregulated", DEG)) %>%
    mutate(DEG = ifelse(log10padj > -log10(padjlevel) & log2FoldChange < Down_FC, "Downregulated", DEG))
  my_pal <- c("#943126",
              "#839192",
              "#A6761D",
              "#1E8449" 
              )
  volcano_up <- volcano[which(volcano$DEG == "Upregulated"),]
  volcano_up <- volcano_up[order(volcano_up$log2FoldChange, decreasing = TRUE), ]
  volcano_down <- volcano[which(volcano$DEG == "Downregulated"),]
  volcano_down <- volcano_down[order(volcano_down$log2FoldChange, decreasing = F), ]
  volcano_NA <- volcano[which(volcano$DEG == "NotSignificance"),]
  volcano_all <- rbind(volcano_up, volcano_down, volcano_NA)
  ntop_rep <- rep("Top", ntop)
  volcano_all$DEG <- c(ntop_rep,
                       volcano_up$DEG[-c(1:ntop)], 
                       ntop_rep, volcano_down$DEG[-c(1:ntop)],
                       volcano_NA$DEG)
  volcano_all$Top <- ifelse(volcano_all$DEG == "Top",
                            volcano_all$GeneName, "")
  g <- ggplot(data = volcano_all, aes(x = log2FoldChange, y = log10padj, color = DEG, fill = DEG, label = Top)) +
    geom_point(size = 2, shape = 21) +
    geom_text(check_overlap = T, vjust = -0.1, nudge_y = 0.1) +
    scale_color_manual(values = my_pal) +
    scale_fill_manual(values = my_pal) +
    theme_classic() +
    labs(x = "Log Fold Change", y = "-log10Padjusted", title = "Volcano PLot")
  g <- g + theme(axis.line = element_line(linetype = "solid"),
                 axis.title = element_text(family = "Times",
                 size = 14, face = "bold"), axis.text = element_text(family = "Times",
                 size = 14, face = "bold", colour = "black"),
                 plot.title = element_text(family = "Times",
                 size = 20, face = "bold", hjust = 0.5)) +
    labs(title = "Volcano Plot")
  pdf("Volcano.pdf", width = 10, height = 8)
  print(g)
  dev.off() 
  return(g)
}



# Heatmap using all DEG 
heatmap_plot <- function(count, 
                         groups){
  DEG_genes <- data.frame(read_csv("diff_result_EdrisSharifRahmani.csv"))[,1]
  hm_df <- counts[which(rownames(counts) %in% DEG_genes), ]
  hm_df <- data.frame(t(hm_df))
  hm_df$class <- groups
  hm_df <- hm_df[order(hm_df$class), ]
  class <- hm_df$class
  hm_df <- as.matrix(scale(t(hm_df[, -ncol(hm_df)])))
  colnames(hm_df) <- class
  ha <- HeatmapAnnotation(
    df = data.frame(Groups = colnames(hm_df)),
    col = list(Groups = c("Tumor" = "#D4AC0D",
                          "Normal" = "#0E6251")),
    annotation_height = unit(4, "mm")
  )
  h <- Heatmap(hm_df,
               bottom_annotation = ha,
               use_raster = TRUE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_column_names = FALSE,
               heatmap_legend_param = list(legend_direction = "vertical",
                                           legend_height = unit(50, "mm"),
                                           grid_width = unit(6, "mm"),
                                           grid_height = unit(50, "cm"),
                                           title = "Eexpression"))
  pdf("Heatmap.pdf", width = 10, height = 8)
  print(h)
  dev.off() 
  return(h)
}

# KEGG pathway
pathway_result <- function(p_adj, top = 10){
  up_down <- data.frame(read_csv("diff_result_EdrisSharifRahmani.csv"))
  up_down$entrez <- mapIds(org.Hs.eg.db, 
                           keys = up_down$GeneName, 
                           column = "ENTREZID", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  kegg_enrich <- enrichKEGG(gene = up_down$entrez)
  kegg_enrich <- kegg_enrich@result
  kegg_enrich <- kegg_enrich[which(kegg_enrich$p.adjust < p_adj), ]
  write_csv(kegg_enrich, "pathway_result_EdrisSharifRahmani.csv")
  kegg_enrich$total <- as.numeric(gsub("[0-9]+/", "", kegg_enrich$GeneRatio))
  kegg_enrich$GeneRatio <- kegg_enrich$Count / kegg_enrich$total
  kegg_enrich <- kegg_enrich[1:top, ]
  g <- ggplot(data = kegg_enrich, aes(x = reorder(Description, GeneRatio, sum),
                             y = GeneRatio, color = p.adjust, size = Count)) +
    geom_point() +
    coord_flip() +
    scale_color_gradient(low = "red", high = "blue") +
    xlab("Pathway") +
    ylab("Gene Ratio") +
    theme_classic() +
    theme(axis.line = element_line(linetype = "solid"),
          axis.title = element_text(family = "Times",
          size = 14, face = "bold"), axis.text = element_text(family = "Times",
          size = 14, face = "bold", colour = "black"),
          plot.title = element_text(family = "Times",
          size = 20, face = "bold", hjust = 0.5))
  pdf("Pathway.pdf", width = 14, height = 8)
  print(g)
  dev.off()
  return(g)
}

deseq_result(deseq, 
             FC = 1, 
             p_adj = 0.05,
             ntop = 10)

volcano_plot(padjlevel = 0.05,
             Up_FC = 1,
             Down_FC = -1,
             ntop = 10) 
heatmap_plot(count = counts, 
             groups = pheno$type)

pathway_result(p_adj = 0.05,
               top = 10)
