##match gene count and methylation intensity
match_expr_methy <- function(gene_count_infor, decay_methy){
  ##select DE gene data
  gene_count <- gene_count_infor[[1]]
  size_factor <- gene_count_infor[[2]]
  gene_name <- rownames(as.character(rownames(gene_count)))
  gene_count <- as.data.frame(cbind(gene_name,gene_count))
  rownames(gene_count) <- NULL
  intersect_gene <- intersect(gene_count$gene_name, decay_methy$gene_name)
  match_data <- data.frame()
  for (i in 1:length(intersect_gene)) {

    methy_name <- intersect_gene[i]
    match_gene <- as.data.frame(gene_count[which(!is.na(match(gene_count$gene_name,  methy_name))),])
    select_methy <- as.data.frame(decay_methy[which(!is.na(match(decay_methy$gene_name, methy_name))),])
    match_methy_expr <- cbind(match_gene, select_methy[,-1])
    match_data <- rbind(match_data, match_methy_expr)

  }
  m6A_express_input <- list(match_data,size_factor)
  return(m6A_express_input)
}
