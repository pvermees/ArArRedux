#' Calculate the irradiation parameter ('J factor')
#'
#' Interpolate the irradiation parameters for the samples
#' given the 40Ar*/39ArK ratios of the samples and fluence monitors
#' 
#' @param R a vector of 40Ar*/39ArK ratios
#' @return an object of class \code{redux} containing, as
#' \code{intercepts}, the 40Ar*/39ArK ratios of the samples, the
#' interpolated J-factors, and the 40K decay constant; and as
#' \code{covmat}: the covariance matrix. All other class properties
#' are inherited from \code{R}.
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' J <- getJfactors(R)
#' plotcorr(J)
#' @export
getJfactors <- function(R){
    RS <- interpolateRJ(R)
    lambda <- R$param$l0
    ts <- R$param$ts
    covmat <- mergematrices(RS$covmat,
              matrix(c(R$param$sl0^2,0,0,R$param$sts^2),nrow=2))
    ns <- length(RS$labels)/2
    out <- RS
    out$intercepts[(ns+1):(2*ns)] <-
        (exp(lambda*ts)-1)/RS$intercepts[(ns+1):(2*ns)]
    out$intercepts[2*ns+1] <- lambda
    out$labels[2*ns+1] <- 'lambda'
    out$thedate[2*ns+1] <- NA
    J <- getJJ(RS$intercepts,ns,lambda,ts)
    out$covmat <- J %*% covmat %*% t(J)
    return(out)
}

interpolateRJ <- function(R){ 
    S <- averagebypos(R,R$Jpos,newlabel="J")
    ji <- getindices(S,pos=R$Jpos)
    si <- getindices(S,pos=R$Jpos,invert=TRUE)
    ns <- length(si)
    J <- matrix(0,nrow=2*ns,ncol=length(S$intercepts))
    for (i in 1:length(si)){ # loop through samples
        J[i,si[i]] <- 1
        px <- S$pos[si[i]]
        dx <- px-R$Jpos
        pm <- R$Jpos[which(dx==max(dx[dx<0]))]
        pM <- R$Jpos[which(dx==min(dx[dx>0]))]
        mi <- getindices(S,pos=pm)
        Mi <- getindices(S,pos=pM)
        J[ns+i,mi] <- (px-pm)/(pM-pm)
        J[ns+i,Mi] <- (pM-px)/(pM-pm)
    }
    out <- R
    samps <- subset(S,si)
    out$thedate <- rep(samps$thedate,2)
    out$irr <- rep(samps$irr,2)
    out$pos <- rep(samps$pos,2)
    out$num <- rep(samps$num,2)
    out$den <- rep(samps$den,2)
    out$nlr <- rep(samps$nlr,2)
    out$labels <- c(samps$labels,paste('J:',samps$labels,sep=''))
    out$intercepts <- J %*% S$intercepts
    out$covmat <- J %*% S$covmat %*% t(J)    
    return(out)
}

getJJ <- function(RS,ns,lambda,ts){
    J <- matrix(0,nrow=2*ns+1,ncol=2*ns+2)
    J[1:ns,1:ns] <- diag(ns)
    J[2*ns+1,2*ns+1] <- 1
    for (i in (ns+1):(2*ns)){
        J[i,i] <- (1-exp(lambda*ts))/(RS[i]^2)
        J[i,2*ns+1] <- ts*exp(lambda*ts)/RS[i]
        J[i,2*ns+2] <- lambda*exp(lambda*ts)/RS[i]
    }
    return(J)
}
