#' Correct for radioactive decay occurred since irradiation
#'
#' Correct for radioactive decay of neutron-induced 37Ar and 39Ar
#' occurred since irradiation
#' 
#' @param X an objects of class \code{redux}
#' @param irr the irradiation schedule
#' @param isotope a string denoting the isotope that needs correcting
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' A <- massfractionation(C,Melbourne$fract)
#' D9 <- decaycorrection(A,Melbourne$irr,"Ar39")
#' plotcorr(D9)
#' @export
decaycorrection <- function(X,irr,isotope){
    out <- X
    D <- getDmatrix(X,irr,isotope)
    J <- Jdecay(X,isotope)
    intercepts <- c(X$intercepts,D$intercepts)
    covmat <- mergematrices(X$covmat,D$covmat)
    out$intercepts <- J %*% intercepts
    out$covmat <- J %*% covmat %*% t(J)
    return(out)
}

getD <- function(P,T,dt,lambda){
    num <- sum(P*dt)
    den <- sum(P*(exp(-lambda*T)-exp(-lambda*(T+dt)))/lambda)
    return(log(num)-log(den))
}

getdDdL <- function(P,T,dt,lambda){
    num <- sum(P*((lambda*T-lambda*dt+1)*exp(-lambda*(T+dt)) -
                  (lambda*T          +1)*exp(-lambda*T)))
    den <- sum(P*(exp(-lambda*(T+dt))-exp(-lambda*T)))
    return(num/den)
}

Jdecay <- function(X,isotope){
    ni <- length(X$intercepts)
    ns <- nruns(X)
    J <- matrix(0,nrow=ni,ncol=ni+ns)
    ii <- 0
    for (i in 1:ns){
        for (j in 1:X$nlr[i]){
            ii <- ii + 1
            J[ii,ii] <- 1
            if (!is.na(X$num[ii]) & X$num[ii]==isotope) J[ii,ni+i] <- 1
            if (!is.na(X$den[ii]) & X$den[ii]==isotope) J[ii,ni+i] <- -1
        }
    }
    return(J)
}

getTdt <- function(irr,thedate){
    dt <- (irr$tout-irr$tin)/(365*24*3600)
    T <- (thedate-irr$tout)/(365*24*3600)
    return(list(T=T,dt=dt))
}

# computes decay correction
getDmatrix <- function(X,irradiations,isotope){
    ns <- nruns(X)
    D <- rep(0,ns)
    dDdL <- rep(0,ns)
    if (isotope=="Ar37"){
        lambda <- X$param$l7
        vlambda <- X$param$sl7^2
    }
    if (isotope=="Ar39"){
        lambda <- X$param$l9
        vlambda <- X$param$sl9^2
    }
    for (i in 1:ns){ # loop through the samples
        irr <- irradiations[[X$irr[i]]]
        Tdt <- getTdt(irr,X$thedate[i])
        D[i] <- getD(irr$P,Tdt$T,Tdt$dt,lambda)
        dDdL[i] <- getdDdL(irr$P,Tdt$T,Tdt$dt,lambda)
    }
    covmatD <- (dDdL * vlambda) %*% t(dDdL)
    return(list(intercepts=D,covmat=covmatD))
}
