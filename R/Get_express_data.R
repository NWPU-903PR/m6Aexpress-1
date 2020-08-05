Get_express_data <- function(INPUT_BAM, TREATED_INPUT_BAM=character(0),annot_file="hg19", 
                            isPairedEnd=FALSE, nthreads=2){
  Input_data <- c(INPUT_BAM, TREATED_INPUT_BAM)
  gene_count <- featureCounts(Input_data,useMetaFeatures=TRUE, annot.inbuilt=annot_file, isPairedEnd=isPairedEnd, nthreads=nthreads)
  counts_data <- gene_count$counts
  countMatrix <- sapply(as.matrix(counts_data), as.numeric)
  gene_countmatrix <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
  return(gene_countmatrix)
}