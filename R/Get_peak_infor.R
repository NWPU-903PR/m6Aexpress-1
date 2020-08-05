##Get peak site and peak information for each gene
Get_peakinfor <- function(IP_BAM, Input_BAM,Treated_IP_BAM, Treated_Input_BAM,GENOME = NA, UCSC_TABLE_NAME = "knownGene", GENE_ANNO_GTF=NA, TXDB=NA, OUTPUT_DIR= NA){
  IP_bam <- c(IP_BAM, Treated_IP_BAM)
  INPUT_bam <- c(Input_BAM, Treated_Input_BAM)
  # Get the annotation file
  if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) & 
                       is.na(TXDB) & is.na(GENE_ANNO_GTF))) {
    op <- options(warn = (-1))
    txdb = makeTxDbFromUCSC(genome = GENOME, tablename = UCSC_TABLE_NAME)
    KeyType="ENTREZID"
    options(op)
  }
  if (suppressWarnings(!is.na(GENE_ANNO_GTF) & is.na(TXDB))) {
    op <- options(warn = (-1))
    txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
    KeyType="ENSEMBL"
    options(op)
  }
  
  # use provided annotation data file
  if (suppressWarnings(!is.na(TXDB))) {
    txdb <- loadDb(TXDB)
    KeyType="ENTREZID"
  }

  ##Get peak site from bam file
  result <- exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, IP_BAM=IP_bam, INPUT_BAM=INPUT_bam,OUTPUT_DIR=OUTPUT_DIR)
  peak_data <- load(paste0(OUTPUT_DIR,"/exomePeak.Rdata"))
  consisten_peak <- paste0(OUTPUT_DIR,"/con_peak.bed") 
  read_peak <- import(consisten_peak)
  read_peak <- as.data.frame(read_peak)
  peak_name <- as.character(read_peak$name)
  tans_name <- select(org.Hs.eg.db, keys=peak_name, columns = c("SYMBOL"),keytype= Keytype)
  select_peak <-cbind( as.character(read_peak$seqnames) , read_peak$start, read_peak$end,read_peak$width, as.character(read_peak$strand), tans_name[,2])
  colnames(select_peak) <- c("seqnames", "start", "end", "width", "strand", "gene_name")
  select_peak <- as.data.frame(select_peak)
  peak <- peak_file[["PEAK"]]
  READS_COUNT <- peak_file[["READS_COUNT"]]
  ## get size factor
  reads_count <- READS_COUNT[,-((ncol(READS_COUNT)-1):ncol(READS_COUNT))]
  totalreads <- colSums(reads_count)
  # get number of consistent peaks
  no_peak=length(peak$loci2peak_consistent[,1])
  # peak_reads_count
  peak_reads_count = READS_COUNT[1:no_peak,]
  no_sample=length(peak_reads_count[1,])-2
  # cut the unnecessary information
  peak_reads_count = READS_COUNT[1:no_peak,1:no_sample]
  # count
  for (ipeak in 1:no_peak) {
    temp=peak$loci2peak_consistent[ipeak,]
    temp2=colSums(READS_COUNT[temp[1]:temp[2],1:no_sample])
    peak_reads_count[ipeak,1:no_sample]=temp2
  }
  # remove the overlapping window effects
  peak_reads_count = round (peak_reads_count * 30 / 200);
  colnames(peak_reads_count) <- c(c(paste0(condition1, "_IP", (1:length(IP_BAM))), paste0(condition2,"_IP",(1:length(Treated_IP_BAM)))), 
                                  c(paste0(condition1, "_Input", (1:length(Input_BAM)))), paste0(condition2, "_Input",(1:length(Treated_Input_BAM))))
  ## get fdr
  log_fdr <- peak$PW$log.fdr
  log_fc <- peak$PW$log.fc
  consisten_peak <- peak$Consistent
  consit_log_fdr <- log_fdr[consisten_peak]
  
  
  # initializa the peak reporting
  no_peak=length(peak$loci2peak_consistent[,1])
  # get peak
  peak_report <- data.frame()
  for (i in 1:no_peak) {
    peak_row_id=peak$loci2peak_consistent[i,]
    
    # batch id
    batch_id=unique(READS_COUNT$batch_id[peak_row_id])
    lg.p=min(peak$PW$log.p[peak_row_id[1]:peak_row_id[2]])/log(10)
    lg.fdr=min(peak$PW$log.fdr[peak_row_id[1]:peak_row_id[2]])/log(10)
    fold_enrchment=exp(max(peak$PW$log.fc[peak_row_id[1]:peak_row_id[2]]))
    
    # get sig digits
    lg.p=signif(lg.p, digits = 3)
    lg.fdr=signif(lg.fdr, digits = 3)
    fold_enrchment=signif(fold_enrchment, digits = 3)
    
    test_result <- data.frame(lg.p=lg.p,lg.fdr=lg.fdr, fold_enrchment=fold_enrchment)               
    peak_report=rbind(peak_report,test_result)
  }
  ## peak site read and log(fdr) 
  peak_site_reads <- cbind(select_peak, peak_reads_count, peak_report[,2:3])
  peak_site_reads <- peak_site_reads[which(!is.na(peak_site_reads$gene_name)), ]
  peak_site_infor <- list(peak_site_reads, totalreads)
  return(peak_site_infor)
  
}
