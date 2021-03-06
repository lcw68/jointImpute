% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_function.R
\name{jointimpute}
\alias{jointimpute}
\title{JointImpute: joint imputation method for microbiome data using DNA and RNA iteratively}
\usage{
jointimpute(
  otu_dna,
  otu_rna,
  covariates_DNA,
  covariates_RNA,
  is_identify,
  epochs = 1,
  coef = 10000
)
}
\arguments{
\item{otu_dna}{A m*n data frame/ matrix that represents OTU count table for DNA}

\item{otu_rna}{A m*n data frame/ matrix that represents OTU count table for RNA}

\item{covariates_DNA}{A n*p data frame/ matrix corresponding to covariates matrix/ meta data for DNA}

\item{covariates_RNA}{A n*p data frame/ matrix corresponding to covariates matrix/ meta data for RNA}

\item{is_identify}{indicator of whether we need identify step}

\item{epochs}{Iteration times for imputation}

\item{coef}{proportional constant for normalizing}
}
\value{
A list of imputed DNA and RNA matrix
}
\description{
The function proposed an imputation method which output the imputed data as DNA and RNA level at the same time.
}
\examples{

covariates_DNA  <- data.frame(batch=factor(meta$batch.DNA), age=meta$age_mos, cariesfree=meta$cariesfree, reads=meta$reads.DNA/1e6)
covariates_RNA = data.frame(batch=factor(meta$batch.RNA),age=meta$age_mos,cariesfree=meta$cariesfree,reads = meta$reads.RNA/1e6)
otu_dna= otu[,,1]
otu_rna= otu[,,2]
otu.imp = jointimpute(otu_dna, otu_rna, covariates_DNA, covariates_RNA,is_identify = TRUE,epochs = 1)



}
