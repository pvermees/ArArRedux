#' Correct for radioactive decay occurred since irradiation
#'
#' Correct for radioactive decay of neutron-induced 37Ar and 39Ar
#' occurred since irradiation
#' 
#' @param x an objects of class \code{redux} or \code{WiscAr}
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
decaycorrection <- function(x,irr,isotope){
    out <- x
    D <- getDmatrix(x,irr,isotope)
    J <- Jdecay(x,isotope)
    intercepts <- c(x$intercepts,D$intercepts)
    covmat <- mergematrices(x$covmat,D$covmat)
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

Jdecay <- function(x,isotope){
    ni <- length(x$intercepts)
    ns <- nruns(x)
    J <- matrix(0,nrow=ni,ncol=ni+ns)
    ii <- 0
    for (i in 1:ns){
        for (j in 1:x$nlr[i]){
            ii <- ii + 1
            J[ii,ii] <- 1
            if (!is.na(x$num[ii]) & x$num[ii]==isotope) J[ii,ni+i] <- 1
            if (!is.na(x$den[ii]) & x$den[ii]==isotope) J[ii,ni+i] <- -1
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
getDmatrix <- function(x,...){ UseMethod("getDmatrix",x) }
getDmatrix.default <- function(x,irradiations,isotope){
    ns <- nruns(x)
    D <- rep(0,ns)
    dDdL <- rep(0,ns)
    if (isotope=="Ar37"){
        lambda <- x$param$l7
        vlambda <- x$param$sl7^2
    }
    if (isotope=="Ar39"){
        lambda <- x$param$l9
        vlambda <- x$param$sl9^2
    }
    for (i in 1:ns){ # loop through the samples
        irr <- irradiations[[x$irr[i]]]
        Tdt <- getTdt(irr,x$thedate[i])
        D[i] <- getD(irr$P,Tdt$T,Tdt$dt,lambda)
        dDdL[i] <- getdDdL(irr$P,Tdt$T,Tdt$dt,lambda)
    }
    covmatD <- (dDdL * vlambda) %*% t(dDdL)
    return(list(intercepts=D,covmat=covmatD))
}
getDmatrix.WiscAr <- function(x,irradiations,clabel){
    w <- Wiscalgas()
    icalgas <- findrunindices(x,prefixes=clabel,invert=FALSE)
    iothers <- findrunindices(x,prefixes=clabel,invert=TRUE)
    ngas <- length(icalgas)/2 # only use hop 101
    noth <- length(iothers)/2 # only use hop 101
    icalgas <- icalgas[1:ngas]
    iothers <- iothers[1:noth]
    nr <- noth+ngas
    nD <- 2*noth+3*ngas
    out <- list()
    D <- rep(0,nD)
    J <- matrix(0,nD,5)
    l7 <- x$param$l7
    l9 <- x$param$l9
    out <- list()
    out$num <- rep(NA,nD)
    out$den <- rep('NA',nD)
    out$nlr <- rep(NA,nr)
    out$thedate <- rep(NA,nr)
    out$irr <- rep(NA,nr)
    out$pos <- rep(NA,nr)
    out$hop <- rep('NA',nr)
    out$labels <- rep(NA,nr)
    j <- 1
    for (i in 1:nr){
        if (i %in% iothers){ # the samples may contain both Ar37 and Ar39
            irr <- irradiations[[x$irr[i]]]
            Tdt <- getTdt(irr,x$thedate[i])
            D[j] <- getD(irr$P,Tdt$T,Tdt$dt,l7)        # Ar37 decay correction
            D[j+1] <- getD(irr$P,Tdt$T,Tdt$dt,l9)      # Ar39 decay correction
            J[j,1] <- getdDdL(irr$P,Tdt$T,Tdt$dt,l7)   # dD7/dl7
            J[j+1,2] <- getdDdL(irr$P,Tdt$T,Tdt$dt,l9) # dD9/dl9
            out$num[j:(j+1)] <- c('Ar37','Ar39')
            out$nlr[i] <- 2
            j <- j + 2
        } else if (i %in% icalgas){ # the calibration gas does not contain Ar37
            dt <- (x$thedate[i]-w$thedate)/(365.25*24*3600)
            D[j] <- w$intercepts[1] - l9*dt   # restored l[0/9] of gas shot i
            D[j+1] <- w$intercepts[2] - l9*dt # restored l[8/9] of gas shot i
            D[j+2] <- w$intercepts[3] - l9*dt # restored l[6/9] of gas shot i
            J[j:(j+2),2] <- -dt   # d[corr]/dl9
            J[j,3] <- 1           # d[09]corr/d09
            J[j+1,4] <- 1         # d[89]corr/d89
            J[j+2,5] <- 1         # d[69]corr/d69
            out$num[j:(j+2)] <- c('Ar40','Ar38','Ar36')
            out$den[j:(j+2)] <- 'Ar39'
            out$nlr[i] <- 3
            j <- j + 3
        }
    }
    E <- matrix(0,5,5)
    E[1,1] <- x$param$sl7^2
    E[2,2] <- x$param$sl9^2
    E[3:5,3:5] <- w$covmat
    covmat <- J %*% E %*% t(J)
    out$labels <- unlist(lapply('DCAL:',paste,x$labels[1:nr]))
    out$intercepts <- D
    out$covmat <- covmat
    return(out)
}
