.onLoad <- function(libname, pkgname) {
	assign("summary", function(object, ...) {
		if(inherits(object, "Trait")) {
			summary.Trait(object, ...)
		}
		else if(inherits(object, "Population")) {
			summary.Population(object, ...)
		}
		else {
			base::summary(object, ...)
		}
	}, envir = .GlobalEnv)
	
	assign("print", function(object, ...) {
		if(inherits(object, "Trait")) {
			print.Trait(object, ...)
		}
		else {
			base::print(object, ...)
		}
	}, envir = .GlobalEnv)
}
