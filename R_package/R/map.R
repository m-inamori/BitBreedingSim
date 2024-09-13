#' @export
summary.Map <- function(map, ...) {
	pop_list <- .Call('_BitBreedingSim_getPopulationInfo', pop)
	cat("Population Summary:\n")
	cat("Num individuals: ", pop_list$num_inds, "\n")
	cat("Num chromomoses: ", pop_list$num_chroms, "\n")
	cat("Num markers: ", pop_list$num_markers, "\n")
}
