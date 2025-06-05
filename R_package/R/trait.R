check_arguments_error_A <- function(trait, h2 = NULL,
										a = NULL, additive_multiplier = NULL) {
	if(is.null(h2) && is.null(a) && is.null(additive_multiplier)) {
		message <- "Error: At least one of 'h2', 'a', or 'additive_multiplier'"
		message <- paste(message, "must be provided.")
		stop(message)
	}
	else if(!is.null(a) && !is.null(additive_multiplier)) {
		message <- "Error: Both 'a' and 'additive_multiplier'"
		message <- paste(message, "cannot be provided simultaneously.")
		stop(message)
	}
	else if(!is.null(a) && !(is.vector(a) && is.numeric(a))) {
		stop("Error: a must be a numeric vector.")
	}
	else if(!is.null(a) && is.vector(a)) {
		num_QTL <- length(trait$loci)
		if(length(a) != num_QTL) {
			stop(sprintf("Error: The length of 'a' must be %d, but it is %d.",
															num_QTL, length(a)))
		}
	}
}

is_at_most_one_valid <- function(a, b, c) {
	counter_valid <- 0
	if(!is.null(a)) { counter_valid <- counter_valid + 1 }
	if(!is.null(b)) { counter_valid <- counter_valid + 1 }
	if(!is.null(c)) { counter_valid <- counter_valid + 1 }
	return(counter_valid <= 1)
}

calc_new_additive_effects <- function(trait, a = NULL,
										additive_multiplier = NULL) {
	if(!is.null(a)) {
		return(a)
	}
	else if(!is.null(additive_multiplier)) {
		return(trait$additives * additive_multiplier)
	}
	else {
		return(trait$additives)
	}
}

calc_new_dominant_effects <- function(trait, d = NULL,
										dominant_multiplier = NULL) {
	if(!is.null(d)) {
		return(d)
	}
	else if(!is.null(dominant_multiplier)) {
		return(trait$dominants * dominant_multiplier)
	}
	else {
		return(trait$dominants)
	}
}

check_arguments_error_AD <- function(trait, h2 = NULL, H2 = NULL,
										a = NULL, additive_multiplier = NULL,
										d = NULL, dominant_multiplier = NULL) {
	if(is.null(h2) && is.null(H2) &&
						is.null(a) && is.null(additive_multiplier) &&
						is.null(d) && is.null(dominant_multiplier)) {
		message <- "Error: At least one of 'h2', 'a', 'additive_multiplier',"
		message <- paste(message, "'d' or 'dominant_multiplier'")
		message <- paste(message, "must be provided.")
		stop(message)
	}
	else if(!is.null(a) && !is.null(additive_multiplier)) {
		message <- "Error: Both 'a' and 'additive_multiplier'"
		message <- paste(message, "cannot be provided simultaneously.")
		stop(message)
	}
	else if(!is_at_most_one_valid(d, dominant_multiplier, H2)) {
		message <- "Error: Both 'd' and 'dominant_multiplier' and 'H2'"
		message <- paste(message, "cannot be provided simultaneously.")
		stop(message)
	}
	else if(!is.null(h2) && !(is.numeric(h2) && length(h2) == 1)) {
		stop("Error: h2 must be a scalar numeric.")
	}
	else if(!is.null(h2) && (h2 < 0.0 || 1.0 < h2)) {
		stop("Error: h2 must be between 0 and 1.")
	}
	else if(!is.null(H2) && !(is.numeric(H2) && length(H2) == 1)) {
		stop("Error: H2 must be a scalar numeric.")
	}
	else if(!is.null(H2) && (H2 < 0.0 || 1.0 < H2)) {
		stop("Error: H2 must be between 0 and 1.")
	}
	else if(!is.null(a) && !(is.vector(a) && is.numeric(a))) {
		stop("Error: a must be a numeric vector.")
	}
	else if(!is.null(d) && !(is.vector(d) && is.numeric(d))) {
		stop("Error: d must be a numeric vector.")
	}
	else if(!is.null(a) && is.vector(a)) {
		num_QTL <- length(trait$loci)
		if(length(a) != num_QTL) {
			stop(sprintf("Error: The length of 'a' must be %d, but it is %d.",
															num_QTL, length(a)))
		}
	}
	else if(!is.null(d) && is.vector(d)) {
		num_QTL <- length(trait$loci)
		if(length(d) != num_QTL) {
			stop(sprintf("Error: The length of 'd' must be %d, but it is %d.",
															num_QTL, length(d)))
		}
	}
	
	# dominat effect is too large
	h2_ = ifelse(!is.null(h2), h2, trait$h2)
	H2_ = ifelse(!is.null(H2), H2, trait$H2)
	if(!is.null(h2) && !is.null(H2) && H2_ < h2_) {
		stop("Error: H2 must be less than h2.")
	}
	add <- calc_new_additive_effects(trait, a, additive_multiplier)
	dom <- calc_new_dominant_effects(trait, d, dominant_multiplier)
	ae2 <- sum(add * add)
	de2 <- sum(dom * dom)
	if((1.0-h2_)*ae2 < de2*h2_/2) {
		stop("Error: dominant effect is too large.")
	}
}

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
