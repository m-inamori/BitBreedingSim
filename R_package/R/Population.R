#' Create origins for a Population object
#'
#' @param num_inds An integer. The number of individuals.
#' @param info An external pointer to a BaseInfo object.
#' @param genotype_ratio A numeric vector of length 3 representing the desired ratio of genotypes (0/0, 0/1, 1/1). 
#'   All elements must be non-negative, and the vector must not be all zeros. 
#'   The values will be normalized to sum to 1.
#'   Default is `c(0.25, 0.5, 0.25)`.
#' @param name_base A string. The base name for individuals.
#' @param names A character vector. The names of individuals. Either `names` or `name_base` must be provided, but not both.

#' @return An external pointer to a Population object.
#' @export
#' @examples
#' # Example 1: Specify individual names directly
#' pop <- create_origins(info, names = c("p1", "p2"))
#' get_individual_names(pop)
#'
#' # Example 2: Generate individual names using name_base and num_inds
#' pop2 <- create_origins(info, num_inds = 2, name_base = "q")
#' get_individual_names(pop2)
#'
#' # Example 3: Specify genotype_ratio explicitly
#' pop3 <- create_origins(info, genotype_ratio = c(1, 1, 1), num_inds = 3, name_base = "g")
#' get_individual_names(pop3)
create_origins <- function(info, genotype_ratio = NULL,
							names = NULL, num_inds = NULL, name_base = NULL) {
	if(is.null(genotype_ratio)) {
		genotype_ratio <- c(0.25, 0.5, 0.25)
	}
	else if(!is.numeric(genotype_ratio) || !is.vector(genotype_ratio) ||
												length(genotype_ratio) != 3) {
		stop("Error: genotype_ratio must be a numeric vector of length 3.")
	}
	else if(any(genotype_ratio < 0)) {
		stop("Error: genotype_ratio must be zero or positive.")
	}
	else if(all(genotype_ratio == 0)) {
		stop("Error: genotype_ratio must not be a zero vector.")
	}
	# normalize
	gratio <- genotype_ratio / sum(genotype_ratio)
	
	if(is.null(name_base) == is.null(names)) {
		stop("Error: Either 'name_base' or 'names' must be specified, but not both.")
	}
	
	if(is.null(names)) {
		names <- generate_names(name_base, num_inds)
	}
	# Check names length vs num_inds and handle mismatch
	if(!is.null(names) && !is.null(num_inds) && length(names) != num_inds) {
		warning(paste("Warning: Length mismatch detected.",
					  sprintf("'names' has %d elements, while 'num_inds' is %d.",
													length(names), num_inds_value),
					  "'num_inds' will be ignored."))
	}
	# If names is NULL, pass name_base. Otherwise, pass names.
	if(is.null(names)) {
		names <- paste0(name_base, 1:num_inds)
	}
	pop <- .Call('_BitBreedingSim_createOrigins', info, gratio, names)
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
#' add_trait_A(info, "Trait1", mean=100.0, h2=0.6, sd=10.0, num_loci = 2)
#' trait <- get_trait(info, 1)
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
#' prog <- cross_randomly(pop, pop, num_inds = 5, name_base = "prog_")
#' haploArrayProg <- create_HaploArray_from_pop(prog)
create_HaploArray_from_pop <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	return(.Call('_BitBreedingSim_createHaploArrayFromPop', pop))
}

#' Set individual names for a Population object
#'
#' This function assigns names to individuals in a Population object. 
#' The `names` argument must be a character vector, and `pop` must be an object 
#' of the "Population" class.
#'
#' @param names A character vector containing the names of individuals to assign.
#' @param pop A Population object where the individual names will be set.
#' @return None. The function modifies the Population object directly.
#' @examples
#' # Create a Population object and set individual names
#' info <- create_base_info()
#' pop <- create_origins(info, num_inds = 2, name_base = "ind")
#' set_individual_names(c("sample1", "sample2"), pop)
#'
#' # Access the updated individual names
#' get_individual_names(pop)
#' @export
set_individual_names <- function(names, pop) {
	if(!is.character(names) || !is.vector(names)) {
		stop("Error: names must be a character vector.")
	}
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	ret <- get_pop_info(pop)
	if(ret$num_individuals != length(names)) {
		stop(paste("Error: length of names (", length(names), 
			") must match the number of individuals in the Population (",
			num, ")."))
	}
	.Call('_BitBreedingSim_setSampleNames', names, pop)
}

#' Generate names for individuals based on a base name and number of individuals
#'
#' This function creates a vector of names for individuals when `names` is not provided. 
#' It uses a base name (`name_base`) and a specified number of individuals (`num_inds`) 
#' to generate names in the format `name_base1`, `name_base2`, ..., up to the total number.
#'
#' @param name_base A character string used as the base for generating names. 
#'                  For example, `name_base = "p"` will generate names like `"p1", "p2", ...`.
#' @param num_inds A positive integer indicating the number of names to generate.
#' @return A character vector of generated names.
#' @examples
#' # Generate 3 names with base "p"
#' generated_names <- generate_names("p", 3)
#' print(generated_names) # Output: "p1", "p2", "p3"
generate_names <- function(name_base, num_inds) {
	# Check inputs
	if(is.null(name_base) || !is.character(name_base)) {
		stop("Error: 'name_base' must be a single character string.")
	}
	if(is.null(num_inds)) {
		stop("Error: When 'name_base' is specified, 'num_inds' must also be specified.")
	}
	if(!is.numeric(num_inds) || num_inds <= 0 || num_inds %% 1 != 0) {
		stop("Error: 'num_inds' must be a positive integer.")
	}
	
    # Generate names
    return(paste0(name_base, 1:num_inds))
}

#' Cross Population randomly
#'
#' This function performs random crossing of individuals from two Population objects 
#' (maternal and paternal parents) to generate a new Population object.
#' For each offspring, a parent is randomly selected from each Population 
#' (one from the maternal Population and one from the paternal Population).
#' You can specify either `name_base` and `num_inds` or provide 
#' a complete `names` vector for the offspring.
#'
#' @param mat_pop A Population object representing the maternal parents.
#'                This should be an object created using the relevant Population functions in the package.
#' @param pat_pop A Population object representing the paternal parents.
#'                This should be an object created using the relevant Population functions in the package.
#' @param names Optional. A character vector containing specific names for individuals. 
#'              If provided, the length of `names` determines the number of offspring.
#' @param num_inds Optional. A positive integer representing the number of offspring. 
#'                 Required if `name_base` is specified.
#' @param name_base Optional. A character string used as the base name for generating individual names. 
#'                  For example, if `name_base = "p"` and `num_inds = 3`, the generated names will be 
#'                  `c("p1", "p2", "p3")`. This parameter must be used in combination with `num_inds`.
#' @param num_threads Optional. A positive integer representing the number of threads to use.
#'                    If not specified or set to 0, the maximum number of available threads is used.
#' @return A Population object with offspring individuals.
#' @examples
#' # Example using name_base and num_inds
#' new_population <- cross_randomly(mothers, fathers, num_inds = 100, name_base = "offspring_")
#'
#' # Example using predefined names
#' names_vector <- c("child1", "child2", "child3")
#' new_population <- cross_randomly(mothers, fathers, names = names_vector)
#'
#' # Summary of the new Population object
#' summary(new_population)
#' @export
cross_randomly <- function(mat_pop, pat_pop, names = NULL,
							num_inds = NULL, name_base = NULL, num_threads = 0) {
	if(is.null(name_base) == is.null(names)) {
		stop("Error: Either 'name_base' or 'names' must be specified, but not both.")
	}
	if(is.null(names)) {
		names <- generate_names(name_base, num_inds)
	}
	# Check names length vs num_inds and handle mismatch
	if(!is.null(names) && !is.null(num_inds) && length(names) != num_inds) {
		warning(paste("Warning: Length mismatch detected.",
					  sprintf("'names' has %d elements, while 'num_inds' is %d.",
													length(names), num_inds_value),
					  "'num_inds' will be ignored."))
	}
	if(num_threads < 1) {
		num_threads <- parallel::detectCores()
	}
	cat("num_threads :", num_threads, "\n")
	
	# If names is NULL, pass name_base. Otherwise, pass names.
	if(is.null(names)) {
		names <- paste0(name_base, 1:num_inds)
	}
	pop <- .Call('_BitBreedingSim_crossPopsRandomly',
									mat_pop, pat_pop, names, num_threads)
	class(pop) <- "Population"
	return(pop)
}

#' Validate the structure and content of the data frame used for crossing information
#'
#' This function checks if the input data frame `df` meets the required 
#' structure for the `cross_by_table` function. It ensures that `df` is a valid 
#' `data.frame`, contains the required columns (`mat`, `pat`, `num`), 
#' and verifies that the `num` column has positive numeric values.
#' If any of these checks fail, the function stops execution and displays 
#' an error message.
#'
#' @param df A data frame. Must include the following columns:
#'           - `mat`: Names of the maternal parents for each cross.
#'           - `pat`: Names of the paternal parents for each cross.
#'           - `num`: Number of offspring to generate for each cross. Must contain positive numeric values.
#' @return No return value. If the validation is successful, the function proceeds silently.
#'         If the validation fails, the function stops execution and displays an error message.
#' @examples
#' # Valid input example
#' df <- data.frame(mat = c("mat1", "mat2"), pat = c("pat1", "pat2"), num = c(1, 2))
#' check_dataframe(df) # Passes silently
#'
#' # Invalid input example: Missing 'num' column
#' df_invalid <- data.frame(mat = c("mat1", "mat2"), pat = c("pat1", "pat2"))
#' check_dataframe(df_invalid) # Throws an error
check_dataframe <- function(df) {
	# Check if df is a data.frame
	if(!is.data.frame(df)) {
		stop("Error: 'df' must be a data.frame.")
	}
	
	# Required columns
	required_columns <- c("mat", "pat", "num")
	
	# Check if all required columns are in df
	missing_columns <- setdiff(required_columns, colnames(df))
	if(length(missing_columns) > 0) {
		stop(sprintf("Error: Missing required column(s): %s",
						paste(missing_columns, collapse = ", ")))
	}
	
	# Check if 'num' column contains valid numbers
	if(!is.numeric(df$num) || any(df$num <= 0)) {
		stop("Error: 'num' column must contain only positive numeric values.")
	}
}

#' Check Parent Existence in Population
#'
#' This function checks whether the maternal and paternal names in the given
#' cross table are present in the specified maternal and paternal populations.
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
	mats <- get_individual_names(mat_pop)
	pats <- get_individual_names(pat_pop)
	
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
		mat <- df$mat[i]
		pat <- df$pat[i]
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
#' This function performs crossing of individuals from two Population objects (maternal and paternal parents)
#' based on specified crossing information provided in a table. Each row of the table defines a specific
#' cross, including the number of offspring to generate for that cross.
#' You can specify either `name_base` or provide a complete `names` vector for the offspring.
#'
#' @param df A data frame containing crossing information. It should include the following columns:
#' \describe{
#'   \item{`mat`}{Names of the maternal parents for each cross.}
#'   \item{`pat`}{Names of the paternal parents for each cross.}
#'   \item{`num`}{Number of offspring to generate for each cross.}
#' }
#' @param mat_pop A Population object representing the maternal parents.
#'                This should be an object created using the relevant Population functions in the package.
#' @param pat_pop A Population object representing the paternal parents.
#'                This should be an object created using the relevant Population functions in the package.
#' @param names Optional. A character vector containing specific names for individuals.
#'              If provided, the length of `names` should match the total number of offspring specified in `df`.
#' @param name_base Optional. A character string used as the base name for generating individual names.
#'                  For example, if `name_base = "prog_"`, the generated names will follow the format 
#'                  `c("prog_1", "prog_2", ...)`.
#' @param num_threads Optional. A positive integer representing the number of threads to use.
#'                    If not specified or set to 0, the maximum number of available threads is used.
#' @return A Population object with offspring individuals.
#' @examples
#' # Example using name_base
#' df <- data.frame(mat = c("mat1", "mat2"), pat = c("pat1", "pat2"), num = c(1, 2))
#' new_population <- cross_by_table(df, mat_pop, pat_pop, name_base = "prog_")
#'
#' # Example using predefined names
#' names_vector <- c("child1", "child2", "child3")
#' new_population <- cross_by_table(df, mat_pop, pat_pop, names = names_vector)
#'
#' # Summary of the new Population object
#' summary(new_population)
#' @export
cross_by_table <- function(df, mat_pop, pat_pop, names = NULL,
									name_base = NULL, num_threads = 0) {
    if(!inherits(mat_pop, "Population")) {
        stop("Error: mat_pop is not a Population object.")
	}
    if(!inherits(pat_pop, "Population")) {
        stop("Error: pat_pop is not a Population object.")
	}
	
	check_dataframe(df)
	num_inds <- sum(df[, "num"])
	
	if(is.null(name_base) == is.null(names)) {
		stop("Error: Either 'name_base' or 'names' must be specified, but not both.")
	}
	if(is.null(names)) {
		names <- generate_names(name_base, num_inds)
	}
	if(!is.null(names) && length(names) != num_inds) {
		stop(sprintf(
			"Error: Length mismatch detected. 'names' has %d elements, while 'num_inds' calculated from the 'num' column in 'df' is %d.",
			length(names), num_inds
		))
	}
	
	if(num_threads < 1) {
		num_threads <- parallel::detectCores()
	}
	cat("num_threads :", num_threads, "\n")
	
	pop <- NULL		# Initialize pop to NULL in case of error
	tryCatch({
		pop <- .Call('_BitBreedingSim_crossPopsByTable', df,
										mat_pop, pat_pop, names, num_threads)
		class(pop) <- "Population"
	}, error = function(e) {
		message(e$message);
		check_parent_existance(df, mat_pop, pat_pop)
	})
	return(pop)
}

#' Get information for a Population object
#'
#' This function retrieves detailed information about a Population object, including the number
#' of individuals, the number of chromosomes, and the number of markers in the population.
#'
#' @param pop A Population object. This should be an object created using the relevant Population functions in the package.
#' @return A list containing:
#'         - `num_individuals`: The number of individuals in the population.
#'         - `num_chromosomes`: The number of chromosomes in the population.
#'         - `num_markers`: The number of markers in the population.
#' @export
#' @examples
#' # Create base information
#' info <- create_base_info()
#' # Create a population with 2 individuals
#' pop <- create_origins(info, num_inds = 2, names = c("p1", "p2"))
#' # Retrieve population information
#' get_pop_info(pop)
get_pop_info <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	num_inds <- .Call('_BitBreedingSim_getNumInds', pop)
	num_chroms <- .Call('_BitBreedingSim_getNumChromsPop', pop)
	num_markers <- .Call('_BitBreedingSim_getNumMarkersPop', pop)
	
	return(list(num_individuals = num_inds,
				num_chromosomes = num_chroms,
				num_markers = num_markers))
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
#' The genotypes are represented in a matrix format where rows correspond to
#' individuals and columns correspond to markers.
#'
#' Genotype encoding:\cr
#' - 0/0 is encoded as -1\cr
#' - 0/1 is encoded as 0\cr
#' - 1/1 is encoded as 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are individuals and columns are markers.
#' @export
#' @examples
#' info <- create_base_info()
#' pop <- create_origins(info, num_inds = 2, names = c("p1", "p2"))
#' geno <- get_genotypes(pop)
#' geno[, 1:5]
get_genotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getGenotypes', pop)
}

#' Get genotypes from a Population object (slow version)
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to
#' individuals and columns correspond to markers.
#'
#' Genotype encoding:\cr
#' - 0/0 is encoded as -1\cr
#' - 0/1 is encoded as 0\cr
#' - 1/1 is encoded as 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are individuals and columns are markers.
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
#' and columns correspond to individuals.
#'
#' Genotype is 0|0, 0|1, 1|0, or 1|1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are individuals and columns are markers.
#' @export
#' @examples
#' info <- create_base_info()
#' pop <- create_origins(info, num_inds = 2, names = c("p1", "p2"))
#' geno <- get_phased_genotypes(pop)
#' geno[, 1:5]
get_phased_genotypes <- function(pop) {
	if(!inherits(pop, "Population")) {
		stop("Error: info is not a BaseInfo object.")
	}
	
	.Call('_BitBreedingSim_getPhasedGenotypes', pop)
}

#' Get phased integer genotypes from a Population object
#'
#' This function retrieves the genotypes from a given Population object.
#' The genotypes are represented in a matrix format where rows correspond to
#' individuals and columns correspond to markers. Each individual has two rows:
#' the first row represents the maternal allele and the second row represents
#' the paternal allele.
#'
#' Genotype is represented as integers: 0 or 1
#'
#' @param pop An external pointer to a Population object.
#' @return A matrix of genotypes where rows are individuals and columns are markers.
#'         Each individual has two rows: the first row is the maternal allele and the
#'         second row is the paternal allele.
#' @export
#' @examples
#' info <- create_base_info()
#' pop <- create_origins(info, num_inds = 2, names = c("p1", "p2"))
#' geno <- get_phased_int_genotypes(pop)
#' geno[, 1:5]
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
#' name_data <- get_individual_names(pop)
#' print(name_data)
get_individual_names <- function(pop) {
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
