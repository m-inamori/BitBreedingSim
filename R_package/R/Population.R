#' Create origins for a Population object
#'
#' @param num_inds An integer. The number of individuals.
#' @param info An external pointer to a BaseInfo object.
#' @param name_base A string. The base name for individuals.
#' @return An external pointer to a Population object.
#' @export
create_origins <- function(num_inds, info, name_base) {
	pop <- .Call('_BitBreedingSim_createOrigins', num_inds, info, name_base)
	class(pop) <- "Population"
	return(pop)
}

#' Create BaseInfo and Population from a VCF file
#'
#' This function takes a VCF object and a seed value, and returns both a Population object
#' and its associated BaseInfo object. It reads the input VCF and initializes the data
#' accordingly. The seed value is used to initialize the pseudo-random number generator
#' for the BaseInfo object. If the seed is set to -1, an appropriate value will be
#' automatically chosen.
#'
#' @param vcf An external pointer to a VCF object.
#' @param seed An integer. The seed value for initializing the BaseInfo object's
#'             pseudo-random number generator. Defaults to -1, which automatically
#'             selects a suitable seed.
#' @return A list containing two elements:
#' \describe{
#'   \item{info}{An external pointer to a BaseInfo object.}
#'   \item{pop}{An external pointer to a Population object.}
#' }
#' @export
#' @examples
#' # Assuming 'vcf_file' is a valid VCF file
#' vcf <- read_VCF(vcf_file)
#' result <- create_info_pop_from_VCF(vcf, seed = 42)
#' summary(result$info)
#' summary(result$pop)
create_info_pop_from_VCF <- function(vcf, seed=-1) {
	if(!inherits(vcf, "VCF")) {
		stop("Error: vcf must be a VCF object.")
	}
	if(!is.numeric(seed) || seed != as.integer(seed)) {
		stop("Error: seed must be an integer.")
	}
	
	result <- .Call('_BitBreedingSim_createInfoAndPopFromVCF', vcf, seed)
	class(result$info) <- "BaseInfo"
	class(result$pop) <- "Population"
	return(result)
}

#' Create a Population object from a HaploArray
#'
#' This function takes a 3-dimensional array (HaploArray) and a BaseInfo object to create a Population object.
#' The HaploArray should have dimensions representing individuals, markers, and alleles (Maternal and Paternal).
#' Allele values must be binary (0 or 1). A value of 0 indicates that the allele matches the reference genome, 
#' while a value of 1 indicates the alternative allele.
#' The function validates the input data and initializes the Population object accordingly.
#'
#' @param haploArray A 3-dimensional array where:
#' \describe{
#'   \item{dim[1]}{Number of individuals.}
#'   \item{dim[2]}{Number of markers.}
#'   \item{dim[3]}{Size must be 2, representing Maternal and Paternal alleles.}
#' }
#' Allele values must be binary (0 or 1), based on the reference genome.
#' @param info An external pointer to a BaseInfo object.
#' @return An external pointer to a Population object.
#' @export
#' @examples
#' # Create a HaploArray
#' haploArray <- array(
#'   data = rbinom(n = 30, size = 1, prob = 0.3),
#'   dim = c(3, 5, 2),
#'   dimnames = list(
#'     paste0("Ind_", 1:3),
#'     paste0("Mrk_", 1:5),
#'     c("Maternal", "Paternal")
#'   )
#' )
#' 
#' # Assuming 'info' is a valid BaseInfo object
#' pop <- create_pop_from_HaploArray(haploArray, info)
#' summary(pop)
create_pop_from_HaploArray <- function(haploArray, info) {
	if(!inherits(info, "BaseInfo")) {
		stop("Error: info must be a BaseInfo object.")
	}
	
	# Check if 'haploArray' is an array
	if(!is.array(haploArray)) {
		stop("Error: haploArray must be an array.")
	}
	# Verify the array has 3 dimensions
	if(length(dim(haploArray)) != 3) {
		stop("Error: haploArray must have exactly 3 dimensions.")
	}
	# Ensure the size of the innermost dimension (Allele) is 2
	if(dim(haploArray)[3] != 2) {
		stop("Error: The third dimension of haploArray must have a size of 2 (Maternal and Paternal).")
	}
	
	info_value <- get_info(info)
	if(info_value$num_markers != dim(haploArray)[2]) {
		stop(
			sprintf(
				"Error: The number of markers is inconsistent.\nhaploArray has %d markers, while info specifies %d markers.",
				dim(haploArray)[2],
				info_value$num_markers
			)
		)
	}
	
	pop <- .Call('_BitBreedingSim_createPopFromHaploArray', haploArray, info)
	class(pop) <- "Population"
	return(pop)
}

#' Create a 3-dimensional HaploArray from a Population object
#'
#' This function converts a Population object into a 3-dimensional array (HaploArray). 
#' The resulting HaploArray contains dimensions representing individuals, markers, 
#' and alleles (Maternal and Paternal).
#' The HaploArray has the following dimensions:
#' \describe{
#'   \item{dim[1]}{Number of individuals.}
#'   \item{dim[2]}{Number of markers.}
#'   \item{dim[3]}{Size is 2, representing Maternal and Paternal alleles.}
#' }
#' The array also includes dimension names such as individual names, marker names, 
#' and allele labels (Maternal and Paternal).
#'
#' Note: This function is the reverse of `create_pop_from_HaploArray`, which creates 
#' a Population object from a HaploArray.
#'
#' @param pop An external pointer to a Population object. This object must be properly initialized.
#' @return A 3-dimensional array (HaploArray) containing the following dimensions:
#' \describe{
#'   \item{Individuals}{Individual names as rows.}
#'   \item{Markers}{Marker names as columns.}
#'   \item{Alleles}{Maternal and Paternal alleles as the third dimension.}
#' }
#' @export
#' @examples
#' info <- create_base_info(num_chroms=2, num_markers=10, cM=100, bp=1e6, seed=2)
#' addTraitA(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
#' trait <- getTrait(info, 1)
#' 
#' haploArray <- array(
#'   data = rbinom(n = 120, size = 1, prob = 0.3),
#'   dim = c(3, 20, 2),
#'   dimnames = list(
#'     paste0("Ind_", 1:3),
#'     paste0("Mrk_", 1:20),
#'     c("Maternal", "Paternal")
#'   )
#' )
#' pop <- create_pop_from_HaploArray(haploArray, info)
#' prog <- cross_randomly(5, pop, pop, "prog_")
#' haploArrayProg <- create_HaploArray_from_pop(prog)
create_HaploArray_from_pop <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	return(.Call('_BitBreedingSim_createHaploArrayFromPop', pop))
}

#' Cross two Population randomly
#'
#' @param num_inds An integer. The number of individuals.
#' @param mat_pop An external pointer to a Population object representing mothers.
#' @param pat_pop An external pointer to a Population object representing fathers.
#' @param name_base A character string. The base name for individuals. If not provided, the 'names' parameter must be specified.
#' @param names A character vector. The specific names for individuals. If not provided, the 'name_base' parameter must be specified.
#' @param num_threads Optional. An integer. The number of threads to be used.
#'                             If not specified, the function will use
#'                             the maximum number of available threads.
#' @return An external pointer to a Population object.
#' @export
#' @examples
#' # Assuming 'mothers' and 'fathers' are valid Population objects
#' new_population <- cross_randomly(100, mothers, fathers, "prog_")
#' summary(new_population)
cross_randomly <- function(num_inds, mat_pop, pat_pop,
							name_base = NULL, names = NULL, num_threads = 0) {
	if(is.null(name_base) && is.null(names)) {
		stop("Either 'name_base' or 'names' must be specified.")
	}
	if(num_threads < 1) {
		num_threads <- parallel::detectCores()
	}
	cat("num_threads :", num_threads, "\n")
	
	# If names is NULL, pass name_base. Otherwise, pass names.
	if(is.null(names)) {
		# Passing an empty character vector if names are not specified
		names <- character(0)
	}
	pop <- .Call('_BitBreedingSim_crossPopsRandomly', num_inds,
								mat_pop, pat_pop, name_base, names, num_threads)
	class(pop) <- "Population"
	return(pop)
}

#' Check Parent Existence in Population
#'
#' This function checks whether the maternal and paternal names in the given cross table
#' are present in the specified maternal and paternal populations.
#'
#' @param df A data.frame representing the cross table containing 'mat' (maternal names) and 'pat' (paternal names) columns.
#' @param mat_pop An external pointer to the maternal Population object.
#' @param pat_pop An external pointer to the paternal Population object.
#' 
#' @return This function does not return a value. It outputs messages if any maternal or paternal names 
#'         in the cross table are not found in the respective populations.
#' 
#' @examples
#' # Assuming 'df', 'mat_pop', and 'pat_pop' are valid objects
#' check_parent_existance(df, mat_pop, pat_pop)
check_parent_existance <- function(df, mat_pop, pat_pop) {
	# Get name data from maternal and paternal Population objects
	mats <- getPopNames(mat_pop)
	pats <- getPopNames(pat_pop)
	
	# using env like a hash table
	name_env <- new.env(hash = TRUE, parent = emptyenv())
	# Store maternal names in the hash table
	for(name in mats$name) {
		assign(name, TRUE, envir = name_env)
	}
	# Store paternal names in the hash table
	for(name in pats$name) {
		assign(name, TRUE, envir = name_env)
	}
	
	# Check each row in the cross table for the existence of maternal
	# and paternal names
	for(i in 1:nrow(df)) {
		mat <- df$mats[i]
		pat <- df$pats[i]
		# If maternal name is not found, output an error message
		if(!exists(mat, envir=name_env, inherits=FALSE)) {
			message(paste(mat, "is not found in maternal population."))
		}
		# If paternal name is not found, output an error message
		if(!exists(pat, envir=name_env, inherits=FALSE)) {
			message(paste(pat, "is not found in paternal population."))
		}
	}
}

#' Cross populations according to a table
#'
#' @param df A data frame. Contains the crossing information with columns for mat, pat, and num.
#'           The 'mat' column represents the maternal population, 'pat' represents the paternal population,
#'           and 'num' represents the number of progenies resulting from the cross.
#' @param mat_pop An external pointer to a Population object representing mothers.
#' @param pat_pop An external pointer to a Population object representing fathers.
#' @param name_base A string. The base name for the new individuals.
#' @param num_threads Optional. An integer. The number of threads to be used.
#'                             If less than 1, the function will use the maximum number of available threads.
#' @return An external pointer to a Population object.
#' @export
#' @examples
#' # Assuming 'mat_pop' and 'pat_pop' are valid inputs
#' mats <- c("mat1", "mat2")
#' pats <- c("pat1", "pat2")
#' nums <- c(1, 2)
#' df <- data.frame(mats, pats, nums)
#' new_population <- cross_by_table(df, mat_pop, pat_pop, "prog_")
#' summary(new_population)
cross_by_table <- function(df, mat_pop, pat_pop, name_base, num_threads = 0) {
	if(num_threads < 1) {
		num_threads <- parallel::detectCores()
	}
	cat("num_threads :", num_threads, "\n")
	pop <- NULL		# Initialize pop to NULL in case of error
	tryCatch({
		pop <- .Call('_BitBreedingSim_crossPopsByTable', df,
									mat_pop, pat_pop, name_base, num_threads)
		class(pop) <- "Population"
	}, error = function(e) {
		message(e$message);
		check_parent_existance(df, mat_pop, pat_pop)
	})
	return(pop)
}

#' Write Population to VCF file
#'
#' This function writes a Population object to a VCF file.
#'
#' @param pop An external pointer to a Population object.
#' @param filename A string specifying the path to the output VCF file.
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object
#' write_VCF(pop, "output.vcf")
write_VCF <- function(pop, filename) {
    if (!inherits(pop, "Population")) {
        stop("Error: pop is not a Population object.")
    }
    tryCatch({
        .Call('_BitBreedingSim_writeVCF', pop, filename)
    }, error = function(e) {
        message("Error: Unable to write to file '", filename, "'. ", e$message)
    })
}

#' Get information for a Population object
#'
#' @param pop An external pointer to a BaseInfo object.
#' @return A list containing the number of individuals
#'         and the number of chromosomes in the population.
#' @export
get_pop_info <- function(pop) {
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
#' phenotypes <- get_phenotypes(pop, 1)
#' print(phenotypes)
get_phenotypes <- function(pop, i) {
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
get_genotypes <- function(pop) {
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
#' @noRd
get_genotypes_naive <- function(pop) {
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
get_phased_genotypes <- function(pop) {
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
get_phased_int_genotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info must be a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getPhasedIntGenotypes', pop)
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
#' selected_pop <- select_pop(pop, c(1, 2, 3))
#' print(selected_pop)
select_pop <- function(pop, indices) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
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

#' Join multiple Population objects
#'
#' @param ... External pointers to Population objects.
#' @return An external pointer to a new Population object
#'         containing the combined individuals.
#' @export
#' @examples
#' # Assuming 'pop1', 'pop2', and 'pop3' are valid Population objects
#' combined_pop <- join_pops(pop1, pop2, pop3)
#' print(combined_pop)
join_pops <- function(...) {
	pops <- list(...)
	if(length(pops) < 2) {
		stop("Error: At least two Population objects are required.")
	}
	
	combined_pop <- pops[[1]]
	
	for(i in 2:length(pops)) {
		combined_pop <- .Call('_BitBreedingSim_joinPop', combined_pop, pops[[i]])
	}
	
	class(combined_pop) <- "Population"
	return(combined_pop)
}

#' Get name data from a Population object
#'
#' @param pop An external pointer to a Population object.
#' @return A data.frame containing the following columns:
#' \describe{
#'   \item{name}{Names from the Population object}
#'   \item{mat}{Maternal names from the Population object}
#'   \item{pat}{Paternal names from the Population object}
#' }
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object
#' name_data <- get_pop_names(pop)
#' print(name_data)
get_pop_names <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	
	.Call('_BitBreedingSim_createNameDataFromPop', pop)
}

#' Summarize the details of a Population object
#'
#' This function extracts and displays key information from a Population object.
#' A Population object represents a collection of individuals, chromosomes, 
#' and genetic markers used in genetic simulations. The summary includes the 
#' number of individuals, chromosomes, and markers in the population, providing 
#' an overview of its structure and size.
#'
#' @param pop An external pointer to a Population object. This object must 
#'            be properly initialized and created using the appropriate 
#'            functions or constructors.
#' @param ... Additional arguments (not currently used).
#' @return Prints a formatted summary of the Population object directly 
#'         to the console.
#' @export
#' @examples
#' # Assuming 'pop' is a valid Population object
#' summary_population(pop)
summary_population <- function(pop, ...) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	
	pop_list <- .Call('_BitBreedingSim_getPopulationInfo', pop)
	cat("Population Summary:\n")
	cat("Num individuals: ", pop_list$num_inds, "\n")
	cat("Num chromomoses: ", pop_list$num_chroms, "\n")
	cat("Num markers: ", pop_list$num_markers, "\n")
}
