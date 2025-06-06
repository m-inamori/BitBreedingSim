% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Population.R
\name{create_info_pop_from_VCF}
\alias{create_info_pop_from_VCF}
\title{Create BaseInfo and Population from a VCF file}
\usage{
create_info_pop_from_VCF(vcf, seed = -1)
}
\arguments{
\item{vcf}{An external pointer to a VCF object.}

\item{seed}{An integer. The seed value for initializing the BaseInfo object's
pseudo-random number generator. Defaults to -1, which automatically
selects a suitable seed.}
}
\value{
A list containing two elements:
\describe{
\item{info}{An external pointer to a BaseInfo object.}
\item{pop}{An external pointer to a Population object.}
}
}
\description{
This function takes a VCF object and a seed value, and returns both a Population object
and its associated BaseInfo object. It reads the input VCF and initializes the data
accordingly. The seed value is used to initialize the pseudo-random number generator
for the BaseInfo object. If the seed is set to -1, an appropriate value will be
automatically chosen.
}
\examples{
# Assuming 'vcf_file' is a valid VCF file
vcf <- read_VCF(vcf_file)
result <- create_info_pop_from_VCF(vcf, seed = 42)
summary(result$info)
summary(result$pop)
}
