#' Create a BaseInfo object
#'
#' This function creates a BaseInfo object.
#'      If the seed is set to -1, a random seed is generated,
#'      resulting in different outcomes each time the function is called.
#'      If a specific seed is provided,
#'      the random number generation will be based on that seed,
#'      ensuring reproducible results.
#'
#' @param chrom_maps A list of data.frames, each representing a chromosome map.
#'        Each data.frame should have two columns: 'cM' for centiMorgans and
#'        'position' for base pair positions.
#'        The list should be named, with each name corresponding to
#'        a chromosome identifier (e.g., "chr1", "chr2", etc.).
#'        If chrom_maps is provided,
#'        the parameters num_chroms, num_markers, cM, and bp are ignored.
#' @param num_chroms An integer. Number of chromosomes.
#'        Ignored if chrom_maps is provided. Default is 10.
#' @param num_markers An integer. Number of markers per chromosome.
#'        Ignored if chrom_maps is provided. Default is 1000.
#' @param cM A numeric. Length of each chromosome in centiMorgans.
#'        Ignored if chrom_maps is provided. Default is 100.
#' @param bp An integer. Length of each chromosome in base pairs.
#'        Ignored if chrom_maps is provided. Default is 1000000.
#' @param seed An integer. A seed for random number generation.
#'        Default is -1, which generates a random seed.
#' @return An external pointer to a BaseInfo object.
#' @export
#' @examples
#' # Create a BaseInfo object with a random seed
#' base_info_random <- createBaseInfo()
#' getInfo(base_info_random)
#'
#' # Create a BaseInfo object with a specific seed for reproducible results
#' base_info_reproducible <- createBaseInfo(seed = 123)
#' getInfo(base_info_reproducible)
#'
#' # Create a chromosome map with 100 cM and 1 Mbp, containing 1000 markers
#' f <- function(x) { (x^3 / (1 + x^2) + 8/5) * 500 / 16 }
#' cM <- sapply(1:1000, function(i) f(i/250 - 2))
#' position <- sapply(1:1000, function(i) i * 1000)
#' chrom_map <- data.frame(cM, position)
#' chrom_maps <- replicate(10, chrom_map, simplify = FALSE)
#' names(chrom_maps) <- paste0("chr", 1:10)
#' info <- createBaseInfo(chrom_maps, seed = 123)
#' getInfo(info)
createBaseInfo <- function(chrom_maps=NULL, num_chroms=10, num_markers=1000,
													cM=100, bp=1000000, seed=-1) {
	if(!is.null(chrom_maps) && !is.list(chrom_maps)) {
		stop("Error: chrom_maps must be a list of data.frames.")
	}
	if(!is.numeric(num_chroms) || num_chroms <= 0) {
		stop("Error: num_chroms must be a positive integer.")
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
	if(is.null(chrom_maps)) {
		info <- .Call('_BitBreedingSim_createBaseInfoCpp', num_chroms,
													num_markers, cM, bp, seed)
	}
	else {
		info <- .Call('_BitBreedingSim_createBaseInfoWithMap', chrom_maps, seed)
	}
	class(info) <- "BaseInfo"
	return(info)
}

#' Get the values of a BaseInfo object
#'
#' @param info An external pointer to a BaseInfo object
#' @return The list of the values of a BaseInfo object
#' @export
getInfo <- function(info) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	num_chroms <- .Call('_BitBreedingSim_getNumChroms', info)
	num_traits <- .Call('_BitBreedingSim_getNumTraits', info)
	
	return(list(num_chroms = num_chroms, num_traits = num_traits))
}

#' Get a trait from a BaseInfo object
#'
#' @param info An external pointer to a BaseInfo object.
#' @param i An integer. The index of the trait to retrieve.
#' @return The trait at the specified index.
#' @export
getTrait <- function(info, i) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	num_traits <- .Call('_BitBreedingSim_getNumTraits', info)
	if(i < 1 || i > num_traits) {
		stop(paste("Index out of bounds: i should be between 1 and", num_traits, "but got", i))
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
#'   map <- getMap(info)
#'   print(map)
#' }
#' @export
getMap <- function(info) {
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
#' @param sd Optional Phenotype standard deviation
#' @param a Optional numeric vector of additive effects
#' @param loci Optional list of loci
#' @param num_loci Number of loci (default is 1)
#' @export
addTraitA <- function(info, name, mean, h2, sd = NULL, a = NULL,
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
	trait <- getTrait(info, num)
	summary.Trait(trait)
}

#' Add Trait with Additive and Dominance Effects to BaseInfo
#'
#' This function adds a trait with additive and dominance effects to the BaseInfo object.
#'
#' @param info External pointer to BaseInfo object
#' @param name Name of the trait
#' @param mean Phenotype mean
#' @param Optional sd Phenotype standard deviation
#' @param Optional h2 Narrow-sense heritability (proportion of variance due to additive genetic effects)
#' @param Optional H2 Broad-sense heritability (proportion of variance due to all genetic effects, can be NULL)
#' @param a Optional numeric vector of additive effects
#' @param d Optional numeric vector of dominance effects
#' @param loci Optional list of loci
#' @param num_loci Number of loci (default is 1)
#' @export
addTraitAD <- function(info, name, mean, sd = NULL, h2 = NULL, H2 = NULL,
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
	trait <- getTrait(info, num)
	summary.Trait(trait)
}
