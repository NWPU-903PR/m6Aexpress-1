##Differential gene detecting
Select_DEgene <- function(gene_count_infor,cond1, cond2,num_cond1, num_cond2,pvalue=0.05,FDR=0.05,out_dir){
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
  sig_DEgene <- DE_gene[(DE_gene$pvalue<pvalue)|(DE_gene$padj<FDR),]
  DE_gene_count <- gene_count_data[which(!is.na(match(gene_count_data$gene_name, rownames(sig_DEgene)))),]
  write.table(sig_DEgene,file=paste(out_dir,"signif_DE_gene.xls",sep="/"), sep="\t",row.names =TRUE,quote = FALSE)
  DE_result <- list(DE_gene_count,size_factor)
  return(DE_result)


}
