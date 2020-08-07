##Differential gene detecting
Select_DEgene <- function(gene_count_infor,cond1, cond2,num_cond1, num_cond2,
                          DIFF_GENE_CUTOFF_PVALUE,DIFF_GENE_cutoff_FDR,DE_CUTOFF_TYPE){
  gene_count_data <- gene_count_infor[[1]]
  gene_count <- gene_count_data[,-1]
  countMatrix <- sapply(as.matrix(gene_count), as.numeric)
  gene_countmatrix <- matrix(countMatrix, nrow=nrow(gene_count), ncol = ncol(gene_count))
  colnames(gene_countmatrix) <- paste0("sample_",1:length(gene_count))
  rownames(gene_countmatrix) <- as.character(gene_count_data$gene_name)
  coldata <- data.frame(condition = c(rep(cond1,num_cond1),rep(cond2,num_cond2)))
  dds <- DESeqDataSetFromMatrix(countData = gene_countmatrix,
                                colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- DESeq(dds)
  dds <- dds[keep,]
  size_factor <- sizeFactors(dds)
  resLFC <- lfcShrink(dds, coef=2)
  DE_gene <- as.data.frame(resLFC)

  if (DE_CUTOFF_TYPE =="padj") {sig_DEgene <- DE_gene[(DE_gene$padj<DIFF_GENE_cutoff_FDR),]}
  if (DE_CUTOFF_TYPE =="pvalue") {sig_DEgene <- DE_gene[(DE_gene$pvalue<DIFF_GENE_CUTOFF_PVALUE),]}

  DE_gene_count <- gene_count_data[which(!is.na(match(gene_count_data$gene_name, rownames(sig_DEgene)))),]
  DE_result <- list(DE_gene_count,size_factor,sig_DEgene)
  return(DE_result)


}
