.onLoad <- function(libname, pkgname) {
	assign("summary", function(object, ...) {
		if(inherits(object, "Trait")) {
			summary_trait(object, ...)
		}
		else if(inherits(object, "Population")) {
			summary_population(object, ...)
		}
		else {
			base::summary(object, ...)
		}
	}, envir = .GlobalEnv)
	
	assign("print", function(object, ...) {
		if(inherits(object, "Trait")) {
			print_trait(object, ...)
		}
		else {
			base::print(object, ...)
		}
	}, envir = .GlobalEnv)
}
