#' Create a BaseInfo object
#'
#' @param seed An integer. A seed for random number generation.
#' @return An external pointer to a BaseInfo object.
#' @export
createBaseInfo <- function(seed) {
	.Call('_BitBreedingSim_createBaseInfo', seed)
}

#' Get the information of a BaseInfo object
#'
#' @param ptr An external pointer to a BaseInfo object
#' @return The list of the values of a BaseInfo object
#' @export
getInfo <- function(ptr) {
	num_chroms <- .Call('_BitBreedingSim_getNumChroms', ptr)
	num_traits <- .Call('_BitBreedingSim_getNumTraits', ptr)
	
	return(list(num_chroms = num_chroms, num_traits = num_traits))
}

#' Set a default trait on a BaseInfo object
#'
#' @param ptr An external pointer to a BaseInfo object
#' @return No return value, called for side effects.
#' @export
setTrait <- function(ptr) {
	.Call('_BitBreedingSim_setTrait', ptr)
}

#' Add a trait with additive and dominance effects to an object
#'
#' @param ptr An external pointer to a BaseInfo object
#' @param num_loci An integer representing the number of loci
#' @param h2 a double representing a narrow-sense heritability.
#' @param H2 a double representing a broad-sense heritability.
#' @return No return value, called for side effects.
#' @export
setTraitADMulti <- function(ptr, num_loci, h2, H2) {
	.Call('_BitBreedingSim_setTraitADMulti', ptr)
}
