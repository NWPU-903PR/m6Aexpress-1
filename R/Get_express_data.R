Get_express_data <- function(INPUT_BAM, TREATED_INPUT_BAM=character(0),annot_file="hg19",
                             species="human",isPairedEnd=FALSE, nthreads=2){
  Input_data <- c(INPUT_BAM, TREATED_INPUT_BAM)
  gene_count <- featureCounts(Input_data,useMetaFeatures=TRUE, annot.inbuilt=annot_file, isPairedEnd=isPairedEnd, nthreads=nthreads)
  counts_data <- gene_count$counts
  countMatrix <- sapply(as.matrix(counts_data), as.numeric)
  gene_countmatrix <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
  colnames(gene_countmatrix) <- paste0("sample_",1:length(Input_data))
  rownames(gene_countmatrix) <- rownames(counts_data)
  conds <- factor(colnames(counts_data))
  cds <- newCountDataSet(countmatrixgene, conds )
  size_factor <-  sizeFactors(estimateSizeFactors( cds ))
  gene_ID <- rownames(gene_countmatrix)
  if(SPECIES="human"){
    org_db <- org.Hs.eg.db
    tans_name <- select(org.Hs.eg.db, keys=gene_ID, columns = c("SYMBOL"),keytype= "ENTREZID")
  }
  if(SPECIES="mouse"){
    org_db <- org.Mm.eg.db
    tans_name <- select(org_db, keys=gene_ID, columns = c("SYMBOL"),keytype= "ENTREZID")
  }
  if(SPECIES="yeast"){
    org_db <- org.Sc.sgd.db
    KeyType="ORF"
    tans_name <- select(org_db, keys=gene_ID, columns = c("GENENAME"),keytype= "ORF")
  }
  gene_name <- as.character(tans_name[,2])
  gene_countsdata <- cbind(gene_name, gene_countmatrix)
  gene_countdata <- as.data.frame(gene_countsdata)
  gene_countdata <- na.omit(gene_countdata)
  rownames(gene_countdata)<-NULL
  Gene_count_infor <- list(gene_countdata, size_factor)
  return(Gene_count_infor)
}
