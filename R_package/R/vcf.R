#' Read a VCF file and return a VCF object
#'
#' This function reads a VCF file from the specified filename and returns
#' a VCF object as an external pointer. If the input is not a valid character
#' string, a message is displayed and NULL is returned.
#'
#' @param filename A character string specifying the path to the VCF file.
#' @return An external pointer to a VCF object, or NULL if an error occurs.
#' @export
#' @examples
#' # Assuming 'example.vcf' is a valid VCF file
#' vcf <- read_VCF("example.vcf")
#' if (is.null(vcf)) {
#'   cat("Failed to read the VCF file.\n")
#' }
read_VCF <- function(filename) {
	if (!is.character(filename) || length(filename) != 1) {
		message("Error: 'filename' must be a single character string.")
		return(NULL)
	}
	
	vcf <- .Call('_BitBreedingSim_readVCF', filename)
	class(vcf) <- "VCF"
	return(vcf)
}

#' Get information from a VCF object
#'
#' This function retrieves detailed information about a VCF object, including the number
#' of individuals and the number of markers contained within the VCF.
#'
#' @param vcf A VCF object. This should be an object created using the relevant VCF functions in the package.
#' @return A list containing:
#'         - `num_individuals`: The number of individuals in the VCF object.
#'         - `num_markers`: The number of markers in the VCF object.
#' @export
#' @examples
#' # Assume 'vcf' is a valid VCF object
#' vcf_info <- get_VCF_info(vcf)
#' print(vcf_info) # Displays the number of individuals and markers
get_VCF_info <- function(vcf) {
	if(!inherits(vcf, "VCF")) {
		stop("Error: vcf must be a VCF object.")
	}
	
	num_inds <- .Call('_BitBreedingSim_getNumIndsVCF', vcf)
	num_markers <- .Call('_BitBreedingSim_getNumMarkersVCF', vcf)
	
	return(list(num_individuals = num_inds,
				num_markers = num_markers))
}

#' Add a Population object to a VCF object
#'
#' This function adds a Population object to an existing VCF object and creates a new VCF object
#' that includes the population information.
#'
#' @param vcf A VCF object. This should be an object created using the relevant VCF functions in the package.
#' @param pop A Population object. This should be an object created using the relevant Population functions in the package.
#' @return A new VCF object that includes the added Population information.
#' @export
#' @examples
#' # Assume 'vcf' is a valid VCF object and 'pop' is a valid Population object
#' new_vcf <- add_pop_to_VCF(vcf, pop)
#' summary(new_vcf) # Inspect the new VCF object
add_pop_to_VCF <- function(vcf, pop) {
	if(!inherits(vcf, "VCF")) {
		stop("Error: vcf must be a VCF object.")
	}
	if(!inherits(pop, "Population")) {
		stop("Error: pop must be a Population object.")
	}
	
	vcf_info <- get_VCF_info(vcf)
	pop_info <- get_pop_info(pop)
	if(vcf_info$num_markers != pop_info$num_markers) {
		stop(sprintf(pasete(
			"Error: Marker count mismatch detected. The VCF object contains",
			"%d markers, but the Population object contains %d markers.",
			"Please ensure both objects have the same",
			"number of markers before proceeding."),
								vcf_info$num_markers, pop_info$num_markers)
		)
	}
	
	new_vcf <- .Call('_BitBreedingSim_addPopToVCF', vcf, pop)
	class(new_vcf) <- "VCF"
	return(new_vcf)
}

#' Write VCF or Population object to a VCF file
#'
#' This function writes either a VCF object or a Population object to a VCF file. If the input is a
#' VCF object, it directly writes the VCF data to the specified file. If the input is a Population
#' object, it generates VCF data from the Population and writes it to the file.
#'
#' @param obj An object of class VCF or Population. The data to write to the VCF file.
#' @param filename A string specifying the path to the output VCF file.
#' @export
#' @examples
#' # Example 1: Write a Population object to a VCF file
#' pop <- create_origins(info, num_inds = 2, names = c("p1", "p2"))
#' write_VCF(pop, "population_output.vcf")
#'
#' # Example 2: Write a VCF object to a VCF file
#' vcf <- read_VCF("input.vcf")
#' write_VCF(vcf, "vcf_output.vcf")
write_VCF <- function(obj, filename) {
	tryCatch({
		if(inherits(obj, "VCF")) {
			.Call('_BitBreedingSim_writeVCF', obj, filename)
		}
		else if(inherits(obj, "Population")) {
			.Call('_BitBreedingSim_writePopToVCF', obj, filename)
		}
		else {
			stop("Error: The object must be either a VCF or a Population object.")
		}
	}, error = function(e) {
		message("Error: Unable to write to file '", filename, "'. ", e$message)
	})
}
