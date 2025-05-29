#' Summarize the details of a Trait object
#'
#' This function provides a formatted summary of a Trait object. A Trait object represents
#' a characteristic or feature in genetic simulations, and contains metadata such as its 
#' name, type, mean, standard deviation, heritability (h2), and information about associated loci.
#' If the trait includes dominant effects, the broad-sense heritability (H2) is also displayed.
#'
#' @param trait A Trait object containing the following elements:
#' \describe{
#'   \item{name}{The name of the trait as a character string.}
#'   \item{type}{The type of the trait (e.g., "Additive Effect Only").}
#'   \item{mean}{The mean value of the trait as a numeric.}
#'   \item{sd}{The standard deviation of the trait as a numeric.}
#'   \item{h2}{The narrow-sense heritability of the trait as a numeric (between 0 and 1).}
#'   \item{H2}{(Optional) The broad-sense heritability of the trait as a numeric (if applicable).}
#'   \item{hasdominants}{A logical value indicating whether the trait includes dominant effects.}
#'   \item{loci}{A list of loci associated with the trait.}
#' }
#' @return Prints a summary of the Trait object directly to the console.
#' @export
#' @examples
#' info <- create_base_info(num_chroms=2, num_markers=10, cM=100, bp=1e6, seed=2)
#' addTraitA(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
#' trait <- getTrait(info, 1)
#' summary_trait(trait)
#' 
#' # Using the generic summary() function
#' summary(trait)
summary_trait <- function(trait) {
	cat("Trait Summary:\n")
	cat("Name: ", trait$name, "\n")
	cat("Type: ", trait$type, "\n")
	cat("Mean: ", trait$mean, "\n")
	cat("SD: ", trait$sd, "\n")
	cat("h2: ", trait$h2, "\n")
	if(trait$hasdominants) {
		cat("H2: ", trait$H2, "\n")
	}
	cat("Number of Loci: ", length(trait$loci), "\n")
}

#' Print detailed information about a Trait object
#'
#' This function prints the detailed characteristics of a Trait object in a human-readable format.
#' A Trait object represents a characteristic or feature in genetic simulations, and this function
#' displays its metadata such as its name, type, mean, standard deviation, heritability (h2),
#' associated loci, additive effects, and dominant effects (if applicable).
#'
#' The loci are formatted as coordinate pairs, representing their positions. Additive and dominant
#' effects are printed as space-separated values for clarity.
#'
#' @param trait A Trait object containing the following elements:
#' \describe{
#'   \item{name}{The name of the trait as a character string.}
#'   \item{type}{The type of the trait (e.g., "quantitative" or "binary").}
#'   \item{mean}{The mean value of the trait as a numeric.}
#'   \item{sd}{The standard deviation of the trait as a numeric.}
#'   \item{h2}{The narrow-sense heritability of the trait as a numeric (between 0 and 1).}
#'   \item{H2}{(Optional) The broad-sense heritability of the trait as a numeric (if applicable).}
#'   \item{hasdominants}{A logical value indicating whether the trait includes dominant effects.}
#'   \item{loci}{A list of loci associated with the trait, formatted as coordinate pairs.}
#'   \item{additives}{A numeric vector representing the additive effects of the trait.}
#'   \item{dominants}{(Optional) A numeric vector representing the dominant effects of the trait (if applicable).}
#' }
#' @param ... Additional arguments (not currently used).
#' @return Prints the detailed information of the Trait object directly to the console.
#' @export
#' @examples
#' info <- create_base_info(num_chroms=2, num_markers=10, cM=100, bp=1e6, seed=2)
#' addTraitA(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
#' trait <- getTrait(info, 1)
#' print_trait(trait)
#' 
#' # Alternatively, use the generic print() function
#' print(trait)
print_trait <- function(trait, ...) {
	cat("Trait Details:\n")
	cat("Name: ", trait$name, "\n")
	cat("Type: ", trait$type, "\n")
	cat("Mean: ", trait$mean, "\n")
	cat("SD: ", trait$sd, "\n")
	cat("h2: ", trait$h2, "\n")
	if(trait$hasdominants) {
		cat("H2: ", trait$H2, "\n")
	}
	
	format_loci <- function(l) {
		paste("(", l$first + 1, ", ", l$second + 1, ")", sep = "")
	}
	
	cat("Loci: ", paste(sapply(trait$loci, format_loci), collapse = " "), "\n")
	cat("Additive Effects: ", paste(trait$additives, collapse = " "), "\n")
	if(trait$hasdominants) {
		cat("Dominant Effects: ", paste(trait$dominants, collapse = " "), "\n")
	}
}
