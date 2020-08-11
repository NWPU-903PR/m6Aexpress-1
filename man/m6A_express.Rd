\name{m6A_express}
\alias{m6A_express}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessment the regulation of m6A methylation on gene expression
}
\description{
The package is designed to predict genes whose expressions are regulaed by m6A methylation based on a negative binomial regression model with pooling prior information across genes to explicitly estimate the sample-to-sample variability between biological replicates and improve the precision and robustness of regression coefficients, which is suitable for MeRIP-seq data in small sample scenario. The inputs of m6A\-express are RNA-seq data BAM files or INPUT control BAM files of MeRIP-seq data, which are used to quantify gene expression, and IP BAM files with paired INPUT BAM files that are used to quantify methylation inensity in m6A\-express model. In addition, users need to provide the origianl gene annotation file, such as GTF or TXDB format file. 
}
\usage{
m6A_Express(control_express,
            treated_express,
            paired = FALSE,
            IP_BAM,
            INPUT_BAM,
            TREATED_IP_BAM,
            TREATED_INPUT_BAM,
            condition1,
            condition2,
            GENOME = NA,
            UCSC_TABLE_NAME = "knownGene",
            GENE_ANNO_GTF=NA,
            TXDB=NA,
            OUTPUT_DIR= NA,
            experiment_name="m6A-Express_output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{control_express}{
the path of all aligned RNA-seq gene expression data files in BAM format under control condition. If the user have no the conventional RNA-seq data under control condition, user can use the Input aligned MeRIP-seq data files in BAM format as the gene expression data.
}
  \item{treated_express}{
the path of all aligned RNA-seq gene expression data files in BAM format under treated condition. If the user have no the conventional RNA-seq data under treated condition, user can use the Input aligned MeRIP-seq data files in BAM format as the gene expression data.
}
  \item{paired}{
logical indicating if paired-end reads are used. If TRUE, fragments (templates or read pairs) will be counted instead of individual reads. FALSE by default.
}
 \item{IP_BAM}{
the path of all IP aligned MeRIP-seq data files in BAM format under control condition
}
 \item{INPUT_BAM}{
the path of all Input aligned MeRIP-seq data files in BAM format under control condition
}
 \item{TREATED_IP_BAM}{
the path of all treated IP aligned MeRIP-seq data files in BAM format under treated condition
}
 \item{TREATED_INPUT_BAM}{
the path of all treated Input aligned MeRIP-seq data files in BAM format under treated condition
}
 \item{condition1}{
A character string is used to give the name of one condition, e.g., condition1="untreated".
}
 \item{condition2}{
A character string is used to give the name of another condition. e.g., condition2="treated".
}
 \item{GENOME}{
A string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will use it directly; otherwise, this function will download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
 \item{UCSC_TABLE_NAME}{
A string, which gives the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
 \item{GENE_ANNO_GTF}{
A gene annotation file in GTF format if available
}
 \item{TXDB}{
An optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
 \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, Trumpet will output the assessment report under the directory of user specified.
}
 \item{EXPERIMENT_NAME}{ a string, which specifies folder name generated in the output directory that contains all the results, default: "m6A-Express_output"
}
\details{
The m6A_Express function is an all-in-one command that performs all the core functions of the m6AExpress R-package.

In order to analysis the regulation of m6A methylation on gene expression, we develop the R package m6AExpress. This package requires user to input gene expression data (RNA-seq data) and MeRIP-seq data. All the input samples are the aligned BAM files. Using the RNA-seq data and MeRIP-seq data, this package can detect the differential expression genes and their corresponding differential methylation under two condition, such as control and treated. Then this package can build the relationship between methylation and gene expression by negative binomial regression model in small sample scenario. And test whether the dynamic m6A methylation level can regulate the change of gene expression.
}
\value{
By default, m6AExpress will output results both
1. Differential gene as EXCELS file on disk (default:"m6A-Express_output") under the specified directory (default: current working directory)
2. Peak calling information of two conditions. This peak calling result is same as exomePeak package. For more information, please click "?exomepeak" in R.
3. The differential m6A methylation peak sites information saved as EXCELS file on disk (default:"m6A-Express_output") under the specified directory (default: current working directory), if it has the differential methylation between two conditions.
4. The gene list saved as EXCELS file on disk (default:"m6A-Express_output") under the specified directory (default: current working directory), whose methlation are significant associated to gene expression.
5. m6A-Express MA-plot

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Teng Zhang <tengzhagn156@126.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

