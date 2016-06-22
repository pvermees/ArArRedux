#' Calculate 40Ar/39Ar ages
#'
#' Calculate 40Ar/39Ar ages from a vector of 40Ar/39Ar-ratios and
#' J-factors
#' 
#' @param RJ an object of class \code{Redux}
#' containing the amalgamated $^{40}$Ar$^*$/$^{39}$Ar$_K$-ratios and
#' J-factors with their covariance matrix
#' @return an object of class \code{results} containing the
#' ages and their covariance matrix
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' J <- getJfactors(R)
#' ages <- getages(J)
#' plotcorr(ages)
#' @export
getages <- function(RJ){
    out <- list()
    class(out) <- "results"
    ns <- (length(RJ$intercepts)-1)/2
    lambda <- utils::tail(RJ$intercepts,n=1)
    out$thedate <- RJ$thedate[1:ns]
    out$labels <- RJ$labels[1:ns]
    out$ages <- log(1+RJ$intercepts[1:ns]*
                    RJ$intercepts[(ns+1):(2*ns)])/lambda
    J <- matrix(0,nrow=ns,ncol=2*ns+1)
    for (i in 1:ns){
        r <- RJ$intercepts[i]
        j <- RJ$intercepts[ns+i]
        J[i,i] <- j/(lambda*(1+j*r))
        J[i,ns+i] <- r/(lambda*(1+j*r))
        J[i,2*ns+1] <- - log(1+j*r)/(lambda^2)
    }
    out$covmat <- J %*% RJ$covmat %*% t(J)
    return(out)
}

#' Summary table
#'
#' Plots the ages and their standard errors
#' 
#' @param object an objct of class \code{\link{results}}
#' @param ... no other arguments
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' summary(ages)[1:5,]
#' @rdname summary
#' @method summary results
#' @export
summary.results <- function(object,...){
    tab <- cbind(object$labels,object$ages,sqrt(diag(object$covmat)))
    colnames(tab) <- c("name","age[Ma]","s.e.[Ma]")
    print(tab)
    return(tab)
}
