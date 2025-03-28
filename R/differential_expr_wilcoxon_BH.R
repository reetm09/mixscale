# Wilcoxon Differential Expression Test using presto (similar to Seurat 5.0)

#install.packages('presto')
#library('presto')
#' @export
#' Wilcoxon Differential Expression Test using presto (similar to Seurat 5.0)
#' 
#' @import presto
#' 
#' @param object An object of class Seurat.
#' @param ident.1 Assay to use for mixscape classification.
#' @param ident.2 Assay data slot to use.
#' @param labels metadata column with target gene labels.
#' @param de.assay Assay to use when performing differential expression analysis.
#' Usually RNA.
#' @param padj.cutoff specify the DE test p-value cutoff (after Benjamini-Hochberg correction) 
#' to select top large-effect DE genes. Default is 0.05. 
#' @param logfc.threshold the log-fold-change threshold to select the large-effect
#' DE genes. Only DE genes with log-fold-change larger than this value will be 
#' selected. Default is 0.25.
#' @param verbose Display messages
#' 
#' @return Returns a list of top DEGs found based on the wilcoxon test and logfc & pvalue cutoff thresholds. 
#' @concept differential_expre_wilcoxon_BH
TopDEGenes_wilcoxBH <- function(object, 
                          ident.1,
                          ident.2 = "NT",
                          labels = 'gene', #in my case perturbed_gene, ident.1 & ident.2 should be part of labels
                          de.assay = "rna",
                          padj.cutoff = 5e-2,
                          logfc.threshold = 0.25,
                          verbose = TRUE
                          
) {
    #print(c(ident.1))
    if (verbose) {
        message("Finding differentially expressed genes")
        #message("Finding differentially expressed genes for " + chr(ident.1) + " vs " + chr(ident.2))
    }
    
    de.genes.df <- data.frame()
    #tryCatch(
       # expr = {
    #print(c(ident.1, "HELLO"))
    #print(de.assay)
    de.genes.df = wilcoxauc(object, labels, groups_use = c(ident.1, ident.2), seurat_assay = de.assay)
    #print(de.genes.df)
    #print(dim(de.genes.df))
    de.genes.df.woNTC = de.genes.df[de.genes.df$group == ident.1, ]
    
    de.genes.df.woNTC  <-  de.genes.df.woNTC [ de.genes.df.woNTC$padj < padj.cutoff, ]
    de.genes.df.woNTC  <- de.genes.df.woNTC[ abs(de.genes.df.woNTC$logFC) > logfc.threshold, ]
       # },
    #    error = function(e) {}
    #)
    object@tools$TopDEGenes_wilcoxBH[[de.assay]][[ident.1]] <- de.genes.df.woNTC
    
    #print(de.genes.df.woNTC)
    #print(c(ident.1, de.genes.df.woNTC$feature))
    return(list(object=object, de.genes=de.genes.df.woNTC$feature))
    
}
