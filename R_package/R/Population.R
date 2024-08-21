#' Create origins for a Population object
#'
#' @param num_inds An integer. The number of individuals.
#' @param info An external pointer to a BaseInfo object.
#' @param name_base A string. The base name for individuals.
#' @return An external pointer to a Population object.
#' @export
createOrigins <- function(num_inds, info, name_base) {
	.Call('_BitBreedingSim_createOrigins', num_inds, info, name_base)
}

#' Cross two Population
#'
#' @param num_inds An integer. The number of individuals.
#' @param mothers An external pointer to a Population object.
#' @param fathers An external pointer to a Population object.
#' @return An external pointer to a Population object.
#' @export
cross <- function(num_inds, mothers, fathers) {
	.Call('_BitBreedingSim_cross', num_inds, mothers, fathers)
}

#' Get information for a Population object
#'
#' @param num_inds An integer. The number of individuals.
#' @param info An external pointer to a BaseInfo object.
#' @param name_base A string. The base name for individuals.
#' @return An external pointer to a Population object.
#' @export
getPopInfo <- function(ptr) {
	num_inds <- .Call('_BitBreedingSim_getNumInds', ptr)
	num_chroms <- .Call('_BitBreedingSim_getPopNumChroms', ptr)
	
	return(list(num_individuals = num_inds, num_chromosomes = num_chroms))
}

#' Get phenotypes from a Population object
#'
#' @param population An external pointer to a Population object.
#' @param i An integer index.
#' @return A vector of phenotypes.
#' @export
getPhenotypes <- function(population, i) {
	.Call('_BitBreedingSim_getPhenotypes', population, i)
}

#' Select individuals from a Population object
#'
#' @param population An external pointer to a Population object.
#' @param indices A vector of indices to select.
#' @return An external pointer to a new Population object containing the selected individuals.
#' @export
selectPop <- function(population, indices) {
  .Call('_BitBreedingSim_selectPop', population, indices)
}
