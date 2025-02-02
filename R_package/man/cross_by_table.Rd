% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Population.R
\name{cross_by_table}
\alias{cross_by_table}
\title{Cross populations according to a table}
\usage{
cross_by_table(df, mat_pop, pat_pop, name_base, num_threads = 0)
}
\arguments{
\item{df}{A data frame. Contains the crossing information with columns for mat, pat, and num.
The 'mat' column represents the maternal population, 'pat' represents the paternal population,
and 'num' represents the number of progenies resulting from the cross.}

\item{mat_pop}{An external pointer to a Population object representing mothers.}

\item{pat_pop}{An external pointer to a Population object representing fathers.}

\item{name_base}{A string. The base name for the new individuals.}

\item{num_threads}{Optional. An integer. The number of threads to be used.
If less than 1, the function will use the maximum number of available threads.}
}
\value{
An external pointer to a Population object.
}
\description{
Cross populations according to a table
}
\examples{
# Assuming 'mat_pop' and 'pat_pop' are valid inputs
mats <- c("mat1", "mat2")
pats <- c("pat1", "pat2")
nums <- c(1, 2)
df <- data.frame(mats, pats, nums)
new_population <- cross_by_table(df, mat_pop, pat_pop, "prog_")
summary(new_population)
}
