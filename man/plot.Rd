% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.timeresolved}
\alias{plot.timeresolved}
\alias{plot.PHdata}
\title{Plot a time resolved mass spectrometry signal}
\usage{
\method{plot}{timeresolved}(x, label, mass, ...)

\method{plot}{PHdata}(x, label, mass, ...)
}
\arguments{
\item{x}{an object of class \code{\link{timeresolved}} or
\code{\link{PHdata}}}

\item{label}{a string with the name of the run}

\item{mass}{a string indicating the isotope of interest}

\item{...}{optional parameters}
}
\description{
Plots the raw signal of a given isotope against time.
}
\examples{
samplefile <- system.file("Samples.csv",package="ArArRedux")
masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
mMC <- loaddata(samplefile,masses)
plot(mMC,"MD2-1a","Ar40")
mPH <- loaddata(samplefile,masses,PH=TRUE)
plot(mPH,"MD2-1a","Ar40")
}
