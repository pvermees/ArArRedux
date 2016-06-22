#' Cl-interference correction
#'
#' Apply the interference correction for the Cl-decay products
#' 
#' @param X an object of class \code{redux}
#' @param irr the irradiation schedule
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' Cl <- clcorrection(Melbourne$X,Melbourne$irr)
#' plotcorr(Cl)
#' @export
clcorrection <- function(X,irr){
    E <- getEmatrix(X,irr)
    return(concat(list(X,E)))
}

getE <- function(P,T,dt,lambda){
    if (any(T<0)) return(0)
    num <- sum(P*(exp(-lambda*(T+dt))-exp(-lambda*T)))
    den <- lambda*sum(P*dt)
    return(log(1+num/den))
}

getdEdL <- function(P,T,dt,lambda){
    if (any(T<0)) return(0)
    num <- sum(P*((lambda*T          +1)*exp(-lambda* T    ) -
                  (lambda*T+lambda*dt+1)*exp(-lambda*(T+dt))))
    den <- lambda*sum(P*(lambda*dt + exp(-lambda*(T+dt)) - exp(-lambda*T)))
    return(num/den)
}

# computes Cl correction
getEmatrix <- function(X,irradiations){
    out <- X
    ns <- nruns(X)
    lpcl <- log(X$param$pcl)
    slpcl <- X$param$spcl/X$param$pcl
    out$intercepts <- rep(lpcl,ns)
    out$num <- rep("Ar36",ns) 
    out$den <- rep("Ar38",ns)
    out$nlr <- rep(1,ns)
    dEdL <- rep(0,ns)
    lambda <- X$param$l6
    for (i in 1:ns){ # loop through the samples
        irr <- irradiations[[X$irr[i]]]
        Tdt <- getTdt(irr,X$thedate[i])
        out$intercepts[i] <- out$intercepts[i] +
                             getE(irr$P,Tdt$T,Tdt$dt,lambda)
        dEdL[i] <- getdEdL(irr$P,Tdt$T,Tdt$dt,lambda)
        out$labels[i] <- paste("Cl:",X$labels[i],sep='')
    }
    J <- cbind(rep(1,ns),dEdL)
    covmatE <- matrix(c(slpcl^2,0,0,X$param$sl6^2),nrow=2)
    out$covmat <- J %*% covmatE %*% t(J)
    return(out)
}

expired <- function(irr,thedate,l7){
    Tdt <- getTdt(irr,thedate)
    return(any(Tdt$T>5*log(2)/l7))
}
