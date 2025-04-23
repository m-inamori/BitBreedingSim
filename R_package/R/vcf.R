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
	return(vcf)
}
