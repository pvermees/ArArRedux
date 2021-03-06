% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/age.R
\name{getages}
\alias{getages}
\title{Calculate 40Ar/39Ar ages}
\usage{
getages(RJ)
}
\arguments{
\item{RJ}{an object of class \code{Redux}
containing the amalgamated $^{40}$Ar$^*$/$^{39}$Ar$_K$-ratios and
J-factors with their covariance matrix}
}
\value{
an object of class \code{results} containing the
ages and their covariance matrix
}
\description{
Calculate 40Ar/39Ar ages from a vector of 40Ar/39Ar-ratios and
J-factors
}
\examples{
data(Melbourne)
R <- get4039(Melbourne$X,Melbourne$irr)
J <- getJfactors(R)
ages <- getages(J)
plotcorr(ages)
}
