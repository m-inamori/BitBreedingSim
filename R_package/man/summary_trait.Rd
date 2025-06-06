% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trait.R
\name{summary_trait}
\alias{summary_trait}
\title{Summarize the details of a Trait object}
\usage{
summary_trait(trait)
}
\arguments{
\item{trait}{A Trait object containing the following elements:
\describe{
\item{name}{The name of the trait as a character string.}
\item{type}{The type of the trait (e.g., "Additive Effect Only").}
\item{mean}{The mean value of the trait as a numeric.}
\item{sd}{The standard deviation of the trait as a numeric.}
\item{h2}{The narrow-sense heritability of the trait as a numeric (between 0 and 1).}
\item{H2}{(Optional) The broad-sense heritability of the trait as a numeric (if applicable).}
\item{hasdominants}{A logical value indicating whether the trait includes dominant effects.}
\item{loci}{A list of loci associated with the trait.}
}}
}
\value{
Prints a summary of the Trait object directly to the console.
}
\description{
This function provides a formatted summary of a Trait object. A Trait object represents
a characteristic or feature in genetic simulations, and contains metadata such as its
name, type, mean, standard deviation, heritability (h2), and information about associated loci.
If the trait includes dominant effects, the broad-sense heritability (H2) is also displayed.
}
\examples{
info <- create_base_info(num_chroms=2, num_markers=10, cM=100, bp=1e6, seed=2)
addTraitA(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
trait <- getTrait(info, 1)
summary_trait(trait)

# Using the generic summary() function
summary(trait)
}
