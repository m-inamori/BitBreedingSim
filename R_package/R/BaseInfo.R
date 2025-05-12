#' Create a BaseInfo object
#'
#' This function creates a BaseInfo object. If the seed is set to -1, a random
#' seed is generated, resulting in different outcomes each time the function
#' is called. If a specific seed is provided, the random number generation will 
#' be based on that seed, ensuring reproducible results.
#'
#' @name create_base_info
#' @import Rcpp
#' @param positions Optional. A list of numeric vectors, where each vector represents 
#' the marker positions for a chromosome in base pairs. The number of elements in 
#' `positions` defines the number of chromosomes and overrides `num_chroms`. If `positions` 
#' is provided, `num_markers` is ignored.
#' @param num_markers Optional. An integer. Number of markers per chromosome. Ignored if 
#' `positions` is provided. Default is 1000.
#' @param bp Optional. An integer. Length of each chromosome in base pairs. Ignored if 
#' `positions` is provided. Default is 100000000.
#' @param chrom_maps Optional. A list of data.frames, each representing a chromosome map.
#' Each data.frame should have two columns: 'cM' for centiMorgans and
#' 'position' for base pair positions. The list should be named, with each name 
#' corresponding to a chromosome identifier (e.g., "chr1", "chr2", etc.).
#' For a given `cM` value, the corresponding `position` can be interpolated or extrapolated,
#' and vice versa. This effectively represents the mapping as a piecewise linear 
#' approximation of the chromosome's relationship between cM and base pair positions.
#' If `chrom_maps` is provided, the parameter `num_chroms` is ignored.
#' @param num_chroms Optional. An integer. Number of chromosomes. Ignored if either 
#' `chrom_maps` or `positions` is provided. Default is 10.
#' @param cM Optional. A numeric. Length of each chromosome in centiMorgans. Ignored if 
#' `chrom_maps` is provided. Default is 100.
#' @param seed Optional. An integer. A seed for random number generation. Default is -1, 
#' which generates a random seed.
#' @return An external pointer to a BaseInfo object.
#' @useDynLib BitBreedingSim, .registration = TRUE
#' @export
#' @examples
#' # Create a BaseInfo object with a random seed
#' base_info_random <- create_base_info()
#' get_info(base_info_random)
#'
#' # Create a BaseInfo object with a specific seed for reproducible results
#' base_info_reproducible <- create_base_info(seed = 123)
#' get_info(base_info_reproducible)
#'
#' # Create marker positions manually
#' num_chr <- 10
#' position <- sapply(1:100, function(i) i * 1000000)
#' positions <- replicate(10, position, simplify = FALSE)
#'
#' # Create a chromosome map with 100 cM and 100 Mbp
#' cM <- c(14.0, 27.9, 40.2, 48.3, 50.0, 51.7, 59.8, 72.1, 86.0, 100.0)
#' chr_position <- sapply(1:10, function(i) i * 10000000)
#' chrom_map <- data.frame(cM, chr_position)
#' chrom_maps <- replicate(10, chrom_map, simplify = FALSE)
#' names(chrom_maps) <- paste0("chr", 1:10)
#' info <- create_base_info(positions = positions, chrom_maps = chrom_maps, seed = 123)
#' get_info(info)
create_base_info <- function(positions=NULL, num_markers=1000, bp=100000000,
							 chrom_maps=NULL, num_chroms=10, cM=100, seed=-1) {
	if(!is.null(positions) && !is.list(positions)) {
		stop("Error: positions must be a list of data.frames.")
	}
	if(!is.null(chrom_maps) && !is.list(chrom_maps)) {
		stop("Error: chrom_maps must be a list of data.frames.")
	}
	if(!is.numeric(num_chroms) || num_chroms <= 0) {
		stop("Error: num_chroms must be a positive integer.")
	}
	if(!is.null(chrom_maps) && length(chrom_maps) != num_chroms) {
		stop("Error: size of chrom_maps and num_chroms must be the same.")
	}
	if(!is.null(positions) && !is.null(chrom_maps) &&
								length(chrom_maps) != length(positions)) {
		stop("Error: sizes of positions and chrom_maps must be the same.")
	}
	if(!is.numeric(num_markers) || num_markers <= 0) {
		stop("Error: num_markers must be a positive integer.")
	}
	if(!is.numeric(cM) || cM <= 0) {
		stop("Error: cM must be a positive number.")
	}
	if(!is.numeric(bp) || bp <= 0) {
		stop("Error: bp must be a positive integer.")
	}
	if(!is.numeric(seed)) {
		stop("Error: seed must be an integer.")
	}
	
	# create default marker positions
	if(is.null(positions)) {
		step <- (bp - 1) %/% (num_markers - 1)	# common difference
		first_term <- bp - (num_markers - 1) * step
		
		# Generate the arithmetic sequence as natural numbers
		position <- seq(from = first_term, to = bp, by = step)
		positions <- replicate(num_chroms, position, simplify = FALSE)
	}
	
	if(is.null(chrom_maps)) {
		info <- .Call('_BitBreedingSim_createBaseInfoCpp', positions, cM, seed)
	}
	else {
		info <- .Call('_BitBreedingSim_createBaseInfoWithMap',
												positions, chrom_maps, seed)
	}
	class(info) <- "BaseInfo"
	return(info)
}

#' Retrieve key values from a BaseInfo object
#'
#' This function extracts essential information from a BaseInfo object. A BaseInfo object 
#' stores metadata related to genetic simulations, such as the number of chromosomes, traits, 
#' and markers. This function is useful for summarizing or verifying the content of a BaseInfo 
#' object during analysis.
#'
#' @param info An external pointer to a BaseInfo object. This object must be properly initialized 
#'             and created using the appropriate BaseInfo constructor or relevant functions.
#' @return A list containing the following elements:
#' \describe{
#'   \item{num_chroms}{An integer representing the total number of chromosomes contained in the BaseInfo object.}
#'   \item{num_traits}{An integer indicating the number of traits managed by the BaseInfo object.}
#'   \item{num_markers}{An integer specifying the total number of markers tracked in the BaseInfo object.}
#' }
#' @export
#' @examples
#' # Assume 'info' is a valid BaseInfo object
#' values <- get_info(info)
#' print(values$num_chroms)  # Print the number of chromosomes
#' print(values$num_traits)  # Print the number of traits
#' print(values$num_markers) # Print the number of markers
get_info <- function(info) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info must be a BaseInfo object.")
	}
	
	num_chroms <- .Call('_BitBreedingSim_getNumChroms', info)
	num_traits <- .Call('_BitBreedingSim_getNumTraits', info)
	num_markers <- .Call('_BitBreedingSim_getNumAllMarkers', info)
	
	return(list(num_chroms = num_chroms,
				num_traits = num_traits,
				num_markers = num_markers))
}

#' Get a trait from a BaseInfo object
#'
#' @param info An external pointer to a BaseInfo object.
#' @param i An integer. The index of the trait to retrieve.
#' @return The trait at the specified index.
#' @export
#' @examples
#' # Assume 'info' is a valid BaseInfo object
#' trait <- get_trait(info, 1)
#' summary(trait)
get_trait <- function(info, i) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	num_traits <- .Call('_BitBreedingSim_getNumTraits', info)
	if(!is.numeric(i) || i %% 1 != 0 || i < 1) {
		stop("Error: i must be a positive integer.")
	}
	if(i > num_traits) {
		stop(paste("Index out of bounds: i should be between 1 and",
												num_traits, "but got", i))
	}
	trait <- .Call('_BitBreedingSim_getTraitCpp', info, i - 1)
	return(trait)
}

#' Retrieve Genetic Map Information
#'
#' This function retrieves the genetic map information from a `BaseInfo` object.
#'
#' @param info An object of class `BaseInfo`. This object contains the genetic map information.
#' @return A list of data frames, each representing a chromosome. Each data frame contains two columns:
#'   - `cM`: The centiMorgan positions of the markers.
#'   - `position`: The base pair positions of the markers.
#' @examples
#' \dontrun{
#'   # Assuming `info` is a valid BaseInfo object
#'   map <- get_map(info)
#'   print(map)
#' }
#' @export
get_map <- function(info) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getMapfromInfo', info)
}

#' Add Trait with Additive Effects to BaseInfo
#'
#' This function adds a trait with additive effects to the BaseInfo object.
#'
#' @param info External pointer to BaseInfo object
#' @param name Name of the trait
#' @param mean Phenotype mean
#' @param h2 Heritability
#' @param sd Optional. Phenotype standard deviation
#' @param a Optional. Numeric vector of additive effects
#' @param loci Optional. List of loci in the form of a data frame with two columns: 'chrom' and 'marker'. For example:
#' \preformatted{
#' chrom <- c(3, 5, 10)
#' marker <- c(1, 2, 1000)
#' loci <- data.frame(chrom, marker)
#' }
#' @param num_loci Optional. Number of loci (default is 1)
#' @export
#' @examples
#' # Assume 'info' is a valid BaseInfo object
#' add_trait_A(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
#' trait <- get_trait(info, 1)
#' summary(trait)
add_trait_A <- function(info, name, mean, h2, sd = NULL, a = NULL,
												loci = NULL, num_loci = 1) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	else if(is.null(sd) && is.null(a)) {
		stop("Error: Both 'sd' and 'a' cannot be NULL at the same time.")
	}
	else if(!is.null(a) && !is.null(loci) && length(a) != dim(loci)[1]) {
		stop("Error: The length of 'a' and 'loci' must be the same.")
	}
	else if(!is.null(a) && length(a) != num_loci) {
		if(num_loci == 1) {
			num_loci = length(a)
		}
		else {
			stop("Error: The length of 'a' and 'num_loci' must be the same.")
		}
	}
	else if(!is.null(loci) && dim(loci)[1] != num_loci) {
		if(num_loci == 1) {
			num_loci = dim(loci)[1]
		}
		else {
			stop("Error: The length of 'loci' and 'num_loci' must be the same.")
		}
	}
	
	add_Trait_A_wrapper(info, name, mean, h2, sd, a, loci, num_loci)
	num = getNumTraits(info)
	trait <- get_trait(info, num)
	summary_trait(trait)
}

#' Add Trait with Additive and Dominance Effects to BaseInfo
#'
#' This function adds a trait with additive and dominance effects to the BaseInfo object.
#'
#' @param info External pointer to BaseInfo object
#' @param name Name of the trait
#' @param mean Phenotype mean
#' @param sd Optional. Phenotype standard deviation
#' @param h2 Optional. Narrow-sense heritability (proportion of variance due to additive genetic effects)
#' @param H2 Optional. Broad-sense heritability (proportion of variance due to all genetic effects, can be NULL)
#' @param a Optional. Numeric vector of additive effects
#' @param d Optional. Numeric vector of dominance effects
#' @param loci Optional. List of loci in the form of a data frame with two columns: 'chrom' and 'marker'. For example:
#' \preformatted{
#' chrom <- c(3, 5, 10)
#' marker <- c(1, 2, 1000)
#' loci <- data.frame(chrom, marker)
#' }
#' @param num_loci Optioal. Number of loci (default is 1)
#' @export
#' @examples
#' # Assume 'info' is a valid BaseInfo object
#' add_trait_AD(info, "Trait1", mean=100.0, h2=0.6, H2=0.7, sd=10.0, num_loci = 2)
#' trait <- get_trait(info, 1)
#' summary(trait)
add_trait_AD <- function(info, name, mean, sd = NULL, h2 = NULL, H2 = NULL,
								a = NULL, d = NULL, loci = NULL, num_loci = 1) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	else if(!is.null(a) && !is.null(d)) {
		if(is.null(sd) && is.null(h2) && is.null(H2)) {
			message <- "Error: When 'a' and 'd' are provided, at least one of"
			message <- paste(message, "'sd', 'h2', or 'H2' must be provided.")
			stop(message)
		}
		else if(is.null(H2) && !is.null(h2)) {
			# H2 must be determined
			sd2 <- a%*%a / 2 / h2
			H2 <- h2 + d%*%d / 4 / sd2
		}
		else if(is.null(H2) && !is.null(sd)) {
			H2 <- (a%*%a / 2 + d%*%d / 4) / (sd * sd)
		}
	}
	else if(!is.null(a)) {
		message <- "Error: When 'a' is provided and 'd' is NULL"
		if(is.null(H2)) {
			message <- paste(message, "'H2' must be provided.")
			stop(message)
		}
		else if(is.null(sd) && is.null(h2)) {
      		message <- paste(message, "either 'sd' or 'h2' must be provided.")
      		stop(message)
		}
		else if(is.null(h2)) {
			h2 <- a%*%a / 2 / (sd * sd)
		}
	}
	else if(!is.null(d)) {
		message <- "Error: When 'd' is provided and 'a' is NULL"
		if(is.null(h2)) {
			message <- paste(message, "'h2' must be provided.")
			stop(message)
		}
		else if(is.null(sd) && is.null(H2)) {
      		message <- paste(message, "either 'sd' or 'H2' must be provided.")
      		stop(message)
		}
		else if(is.null(H2)) {
			H2 <- h2 + d%*%d / 4 / (sd * sd)
		}
	}
	else {
		if(is.null(h2) || is.null(sd) || is.null(H2)) {
			message <- "Error: When 'd' and 'a' are NULL"
			message <- paste(message, "'sd', 'h2', and 'H2' must be provided.")
			stop(message)
		}
	}
	
	if(!is.null(a) && !is.null(loci) && length(a) != dim(loci)[1]) {
		stop("Error: The length of 'a' and 'loci' must be the same.")
	}
	else if(!is.null(d) && !is.null(loci) && length(d) != dim(loci)[1]) {
		stop("Error: The length of 'd' and 'loci' must be the same.")
	}
	else if(!is.null(a) && !is.null(d) && length(a) != length(d)) {
		stop("Error: The length of 'a' and 'd' must be the same.")
	}
	else if(!is.null(a) && length(a) != num_loci) {
		if(num_loci == 1) {
			num_loci = length(a)
		}
		else {
			stop("Error: The length of 'a' and 'num_loci' must be the same.")
		}
	}
	else if(!is.null(d) && length(d) != num_loci) {
		if(num_loci == 1) {
			num_loci = length(d)
		}
		else {
			stop("Error: The length of 'd' and 'num_loci' must be the same.")
		}
	}
	else if(!is.null(loci) && dim(loci)[1] != num_loci) {
		if(num_loci == 1) {
			num_loci = dim(loci)[1]
		}
		else {
			stop("Error: The length of 'loci' and 'num_loci' must be the same.")
		}
	}
	
	if(H2 > 1) {
		stop("Error: H2 must be less than or equal to 1.")
	}
	add_Trait_AD_wrapper(info, name, mean, sd, h2, H2, a, d, loci, num_loci)
	num = getNumTraits(info)
	trait <- get_trait(info, num)
	summary_trait(trait)
}
