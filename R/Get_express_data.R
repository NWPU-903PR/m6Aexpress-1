Get_express_data <- function(INPUT_BAM, TREATED_INPUT_BAM=character(0),annot_file="hg19", 
                            isPairedEnd=FALSE, nthreads=2){
  Input_data <- c(INPUT_BAM, TREATED_INPUT_BAM)
  gene_count <- featureCounts(Input_data,useMetaFeatures=TRUE, annot.inbuilt=annot_file, isPairedEnd=isPairedEnd, nthreads=nthreads)
  counts_data <- gene_count$counts
  countMatrix <- sapply(as.matrix(counts_data), as.numeric)
  gene_countmatrix <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
  gene_ID <- as.character(rownames(countmatrixgene)) 
  tans_name <- select(org.Hs.eg.db, keys=gene_ID, columns = c("SYMBOL"),keytype="ENTREZID")
  gene_name <- as.character(tans_name[,2]) 
  gene_countsdata <- cbind(gene_name, countmatrixgene)
  gene_countdata <- as.data.frame(gene_countsdata)
  gene_countdata <- na.omit(gene_countdata)
  rownames(gene_countdata)<-NULL
  return(gene_countdata)
}
