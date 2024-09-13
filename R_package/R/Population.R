#' Create origins for a Population object
#'
#' @param num_inds An integer. The number of individuals.
#' @param info An external pointer to a BaseInfo object.
#' @param name_base A string. The base name for individuals.
#' @return An external pointer to a Population object.
#' @export
createOrigins <- function(num_inds, info, name_base) {
	pop <- .Call('_BitBreedingSim_createOrigins', num_inds, info, name_base)
	class(pop) <- "Population"
	return(pop)
}

#' Cross two Population
#'
#' @param num_inds An integer. The number of individuals.
#' @param mothers An external pointer to a Population object.
#' @param fathers An external pointer to a Population object.
#' @param name_base A string. The base name for individuals.
#' @param num_threads optional An integer. The number of threads to be used.
#'                             If not specified, the function will use
#'                             the maximum number of available threads.
#' @return An external pointer to a Population object.
#' @export
#' @examples
#' # Assuming 'mothers' and 'fathers' are valid Population objects
#' new_population <- cross(100, mothers, fathers, "prog_")
#' summary(new_population)
cross <- function(num_inds, mothers, fathers, name_base, num_threads = 0) {
	if(num_threads < 1) {
		num_threads <- parallel::detectCores()
	}
	cat("num_threads :", num_threads, "\n")
	pop <- .Call('_BitBreedingSim_crossPops', num_inds, mothers, fathers,
														name_base, num_threads)
	class(pop) <- "Population"
	return(pop)
}

#' Get information for a Population object
#'
#' @param pop An external pointer to a BaseInfo object.
#' @return A list containing the number of individuals
#'         and the number of chromosomes in the population.
#' @export
getPopInfo <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	num_inds <- .Call('_BitBreedingSim_getNumInds', pop)
	num_chroms <- .Call('_BitBreedingSim_getNumChromsPop', pop)
	
	return(list(num_individuals = num_inds, num_chromosomes = num_chroms))
}

#' Get phenotypes from a Population object
#'
#' @param pop An external pointer to a Population object.
#' @param i An integer index representing the trait for which phenotypes are to
#'        be retrieved. The index should be between 1
#'        and the total number of traits available in the Population object.
#' @return A vector of phenotypes.
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object and trait index 1 is valid
#' phenotypes <- getPhenotypes(pop, 1)
#' print(phenotypes)
getPhenotypes <- function(pop, i) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop is not a BaseInfo object.")
	}
	
	pop_list <- .Call('_BitBreedingSim_getPopulationInfo', pop)
	num_traits <- pop_list$num_traits
	cat("i :", i, "num_traits :", num_traits, "\n")
	if(i < 1 || i > num_traits) {
		stop(paste("Index out of bounds: i should be between 1 and",
												num_traits, "but got", i))
	}
	.Call('_BitBreedingSim_getPhenotypesCpp', pop, i)
}

#' Get genotypes from a Population object
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to samples
#' and columns correspond to markers.
#'
#' Genotype encoding:\cr
#' - 0/0 is encoded as -1\cr
#' - 0/1 is encoded as 0\cr
#' - 1/1 is encoded as 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are samples and columns are markers.
#' @export
getGenotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getGenotypes', pop)
}

#' Get genotypes from a Population object (slow version)
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to samples
#' and columns correspond to markers.
#'
#' Genotype encoding:\cr
#' - 0/0 is encoded as -1\cr
#' - 0/1 is encoded as 0\cr
#' - 1/1 is encoded as 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are samples and columns are markers.
#' @export
getGenotypesNaive <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getGenotypes_naive', pop)
}

#' Get phased genotypes from a Population object
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to markers
#' and columns correspond to samples.
#'
#' Genotype is 0|0, 0|1, 1|0, or 1|1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are samples and columns are markers.
#' @export
getPhasedGenotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getPhasedGenotypes', pop)
}

#' Get phased integer genotypes from a Population object
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to samples
#' and columns correspond to markers. Each sample has two rows: the first row
#' represents the maternal allele and the second row represents the paternal allele.
#'
#' Genotype is represented as integers: 0 or 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are samples and columns are markers.
#'         Each sample has two rows: the first row is the maternal allele and the
#'         second row is the paternal allele.
#' @export
getPhasedIntGenotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getPhasedGenotypes', pop)
}

#' Select individuals from a Population object
#'
#' @param pop An external pointer to a Population object.
#' @param indices A vector of integer indices representing the individuals
#'        to be selected.
#' @return An external pointer to a new Population object
#'        containing the selected individuals.
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object and indices are valid
#' selected_pop <- selectPop(pop, c(1, 2, 3))
#' print(selected_pop)
selectPop <- function(pop, indices) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop is not a Population object.")
	}
	
	num_inds <- .Call('_BitBreedingSim_getNumInds', pop)
	if(any(indices < 1) || any(indices > num_inds)) {
		error_indices <- indices[indices < 1 | indices > num_inds]
		stop(paste("Index out of bounds: indices should be between 1 and",
					num_inds, "but got", paste(error_indices, collapse = ", ")))
    }
    new.pop <- .Call('_BitBreedingSim_selectPop', pop, indices)
    class(new_pop) <- "Population"
    return (new_pop)
}

#' Get name data from a Population object
#'
#' @param pop An external pointer to a Population object.
#' @return A data.frame containing the names, maternal names,
#'         and paternal names from the Population object.
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object
#' name_data <- getPopNames(pop)
#' print(name_data)
getPopNames <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop is not a Population object.")
	}
	
	.Call('_BitBreedingSim_createNameDataFromPop', pop)
}

#' @export
summary.Population <- function(pop, ...) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop is not a Population object.")
	}
	
	pop_list <- .Call('_BitBreedingSim_getPopulationInfo', pop)
	cat("Population Summary:\n")
	cat("Num individuals: ", pop_list$num_inds, "\n")
	cat("Num chromomoses: ", pop_list$num_chroms, "\n")
	cat("Num markers: ", pop_list$num_markers, "\n")
}
