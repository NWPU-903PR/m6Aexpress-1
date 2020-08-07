# m6A-express: Predicting context-specific m6A regulation of gene expression
m6A-express is  a tool to identify m6A methylation significantly regulated expression genes in specific context



# Installation Instructions
The Trumpet package can be installed by the following R commands:
> devtools::install_github("tengzhang123/m6A-express")

> library(m6A-express)

# Usage Example
The following command code will show how to use this package and output m6A methylation regulated expression gene in excel files
## In basic mode: pooled samples from one or multiple conditions together and identify m6A regulated expression gene in a specific context.
### Input the RNA-seq data and MeRIP-seq data from BAM files.
> f1 <- system.file("extdata", "IP1.bam", package="m6A-express")

> f2 <- system.file("extdata", "IP2.bam", package="m6A-express")

> f3 <- system.file("extdata", "IP3.bam", package="m6A-express")

> f4 <- system.file("extdata", "IP4.bam", package="m6A-express")

> f5 <- system.file("extdata", "Input1.bam", package="m6A-express")
 
> f6 <- system.file("extdata", "Input2.bam", package="m6A-express")

> f7 <- system.file("extdata", "Input3.bam", package="m6A-express")
 
> f8 <- system.file("extdata", "Input4.bam", package="m6A-express")
