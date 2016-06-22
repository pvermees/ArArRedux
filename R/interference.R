#' define the interference corrections
#'
#' create a new object of class \code{logratios} containing
#' the interferences from neutron reactions on Ca and K
#'
#' @param intercepts a vector with logratios
#' @param covmat the covariance matrix of the logratios
#' @param num a vector of strings marking the numerator isotopes of
#' \code{intercepts}
#' @param den a vector of strings marking the denominator isotopes of
#' \code{intercepts}
#' @param irr an object of class \code{irradiations}
#' @param label a string with a name which can be used to identify the
#' interference data in subsequent calculations
#' @return an object of class \code{logratios}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' irrfile <- system.file("irradiations.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' X <- read(samplefile,masses,blabel="EXB#",Jpos=c(3,15))
#' irr <- loadirradiations(irrfile)
#'# assume log(36Ar/37Ar) = log(39Ar/37Ar) = 1 in co-irradiate Ca-salt
#'# with variances of 0.0001 and zero covariances
#' ca <- interference(intercepts=c(1,1),
#'                    covmat=matrix(c(0.001,0,0,0.001),nrow=2),
#'                    num=c("Ar39","Ar36"),den=c("Ar37","Ar37"),
#'                    irr=X$irr[1],label="Ca-salt")
#'# assume log(39Ar/40Ar) = 4.637788 in co-irradiate K-glass
#'# with variance 7.9817e-4
#' k <- interference(intercepts=4.637788,covmat=7.9817e-4,
#'                   num="Ar39",den="Ar40",irr=X$irr[1],
#'                   label="K-glass")
#' ages <- process(X,irr,ca=ca,k=k)
#' summary(ages)
#' @export
interference <- function(intercepts,covmat,num,den,irr,label){
    out <- list()
    class(out) <- "logratios"
    out$irr <- irr
    out$pos <- 0
    out$labels <- label
    out$thedate <- as.numeric(as.Date("1970-01-01 00:00:00"))
    out$num <- num
    out$den <- den
    out$nlr <- length(intercepts)
    out$intercepts <- intercepts
    out$covmat <- covmat
    return(out)
}
