##Differential gene detecting
Select_DEgene <- function(gene_count,cond1, cond2,num_cond1, num_cond2,out_dir){
  coldata <- data.frame(condition = c(rep(cond2,num_cond2),rep(cond1,num_cond1)))
  dds <- DESeqDataSetFromMatrix(countData = gene_count,
                                colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- DESeq(dds)
  dds <- dds[keep,]
  size_factor <- sizeFactors(dds)
  resLFC <- lfcShrink(dds, coef=2)
  resSig <- subset(resLFC, padj < 0.05)
  DE_gene <- as.data.frame(resSig)
  write.table(DE_gene,file=paste(out_dir,"DE_gene.xls",sep="/"), sep="\t",row.names =TRUE,quote = FALSE)
  DE_result <- list(DE_gene, size_factor)
  return(DE_result)
  
  
}