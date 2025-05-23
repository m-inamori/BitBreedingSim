% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Population.R
\name{create_HaploArray_from_pop}
\alias{create_HaploArray_from_pop}
\title{Create a 3-dimensional HaploArray from a Population object}
\usage{
create_HaploArray_from_pop(pop)
}
\arguments{
\item{pop}{An external pointer to a Population object. This object must be properly initialized.}
}
\value{
A 3-dimensional array (HaploArray) containing the following dimensions:
\describe{
\item{Individuals}{Individual names as rows.}
\item{Markers}{Marker names as columns.}
\item{Alleles}{Maternal and Paternal alleles as the third dimension.}
}
}
\description{
This function converts a Population object into a 3-dimensional array (HaploArray).
The resulting HaploArray contains dimensions representing individuals, markers,
and alleles (Maternal and Paternal).
The HaploArray has the following dimensions:
\describe{
\item{dim\link{1}}{Number of individuals.}
\item{dim\link{2}}{Number of markers.}
\item{dim\link{3}}{Size is 2, representing Maternal and Paternal alleles.}
}
The array also includes dimension names such as individual names, marker names,
and allele labels (Maternal and Paternal).
}
\details{
Note: This function is the reverse of \code{create_pop_from_HaploArray}, which creates
a Population object from a HaploArray.
}
\examples{
info <- create_base_info(num_chroms=2, num_markers=10, cM=100, bp=1e6, seed=2)
add_trait_A(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
trait <- get_trait(info, 1)

haploArray <- array(
  data = rbinom(n = 120, size = 1, prob = 0.3),
  dim = c(3, 20, 2),
  dimnames = list(
    paste0("Ind_", 1:3),
    paste0("Mrk_", 1:20),
    c("Maternal", "Paternal")
  )
)
pop <- create_pop_from_HaploArray(haploArray, info)
prog <- cross_randomly(pop, pop, num_inds = 5, name_base = "prog_")
haploArrayProg <- create_HaploArray_from_pop(prog)
}
