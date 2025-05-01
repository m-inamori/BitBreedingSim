#' Summarize the details of a Map object
#'
#' This function extracts and displays key information from a Map object. 
#' A Map object represents the genetic map used in simulations, containing 
#' details about chromosomes and points (e.g., markers or loci) along the map. 
#' The summary provides an overview of the structure and size of the map.
#'
#' @param map An external pointer to a Map object. This object must be properly 
#'            initialized and created using the appropriate functions.
#' @param ... Additional arguments (not currently used).
#' @return Prints a formatted summary of the Map object directly to the console.
#' @export
#' @examples
#' # Assuming 'map' is a valid Map object
#' summary_map(map)
summary_map <- function(map, ...) {
	if(!inherits(map, "Map")) {
		stop("Error: map must be a Map object.")
	}
	
	map_list <- .Call('_BitBreedingSim_getMapInfo', map)
	cat("Map Summary:\n")
	cat("Num chromomoses: ", map_list$num_chroms, "\n")
	cat("Num points: ", map_list$num_points, "\n")
}

#' Extract chromosome and position data from a Map object
#'
#' This function retrieves detailed data from a Map object, including chromosome names and 
#' marker positions, and returns the data as a data.frame for easier manipulation.
#'
#' @param map An external pointer to a Map object. This object must be properly initialized 
#'            and created using appropriate functions.
#' @return A data.frame containing the following columns:
#' \describe{
#'   \item{chrom}{A character vector containing the names of chromosomes.}
#'   \item{position}{An integer vector specifying the positions of markers along the chromosomes.}
#' }
#' @export
#' @examples
#' # Assuming 'map' is a valid Map object
#' map_data <- extract_map_data(map)
#' head(map_data) # Display the first few rows of the data.frame
extract_map_data <- function(map) {
	if(!inherits(map, "Map")) {
		stop("Error: map must be a Map object.")
	}
	
	return(.Call('_BitBreedingSim_getMapCpp', map))
}
