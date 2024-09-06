#' @export
summary.Trait <- function(trait, ...) {
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

#' @export
print.Trait <- function(trait, ...) {
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
