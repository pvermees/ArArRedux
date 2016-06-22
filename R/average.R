#' Calculate the arithmetic mean
#'
#' Calculate the arithmetic mean of some logratio data
#' 
#' @param x an object of class \code{redux} or \code{logratios}
#' @param i (optional) vector of sample indices
#' @param newlabel (optional) string with the new label to be assigned
#' to the averaged values
#' @return an object of the same class as \code{x}
#' @examples
#' data(Melbourne)
#' K <- average(Melbourne$X,grep("K:",Melbourne$X$labels),newlabel="K-glass")
#' plotcorr(K)
#' @export
average <- function(x,i=NULL,newlabel=NULL){
    if (length(i)==0) return(x)
    if (is.null(i)) i <- 1:nruns(x)
    if (is.null(newlabel)) newlabel <- as.character(x$labels[i[1]])
    j <- which(!((1:nruns(x)) %in% i))
    out <- subset(x,c(i[1],j))
    out$labels[1] <- newlabel
    J <- Jmean(x$nlr,i)
    out$intercepts <- J %*% x$intercepts
    out$covmat <- J %*% x$covmat %*% t(J)
    return(out)
}

#' Average all the data collected on the same day.
#'
#' This function is useful for grouping a number of replicate air
#' shots or calibration experiments
#' 
#' @param x an object of class \code{timeresolved}, \code{logratios},
#' \code{PHdata} or \code{redux}
#' @param newlabel a string with the new label that should be given to
#' the average
#' @return an object of the same class as x
#' @examples
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' dlabels <- c("H1","AX","L1","L2")
#' md <- loaddata(dfile,dlabels,PH=TRUE)
#' ld <- fitlogratios(blankcorr(md))
#' d <- averagebyday(ld,"DCAL")
#' plotcorr(d)
#' @export
averagebyday <- function(x,newlabel){
    out <- x
    thedays <- rle(theday(x$thedate))$values
    for (i in 1:length(thedays)){
        j <- which(theday(out$thedate)==thedays[i])
        thelabel <- paste(newlabel,i,sep="-")
        out <- average(out,j,thelabel)
    }
    return(out)
}

averagebypos <- function(X,pos,newlabel){
    out <- X
    for (i in 1:length(pos)){
        j <- which(out$pos==pos[i])
        thelabel <- paste(newlabel,i,sep="-")
        out <- average(out,j,thelabel)
    }
    return(out)
}

#' Calculate the weighted mean age
#'
#' Computes the error weighted mean and MSWD of some samples taking
#' into covariances.
#' 
#' @param ages an object of class \code{results}
#' @param prefix is either a string with the prefix of the
#' samples that need to be averaged, or a vector of sample names.
#' @return a list with items:
#'  
#' \code{avgt}: the weighted mean age\cr
#' \code{err}: the standard error of \code{avgt}\cr
#' \code{MSWD}: the Mean Square of the Weighted Deviates
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' weightedmean(ages,"MD2-1")
#' @export
weightedmean <- function(ages,prefix=NULL){
    if (is.null(prefix)){
        slabs <- ages$labels
    } else if (length(prefix)==1){
        slabs <- ages$labels[grep(prefix,ages$labels)]
    } else {
        slabs <- prefix
    }
    subs <- subset(ages,labels=slabs)
    n <- nruns(subs)
    W <- rep(1,n) # design matrix
    invSigma <- solve(subs$covmat)
    vart <- 1/(W %*% invSigma %*% W)
    avgt <- vart*(W %*% invSigma %*% subs$ages)
    mswd <- ((subs$ages - avgt) %*% solve(subs$covmat) %*%
             (subs$ages - avgt))/(n-1)    
    return(list(avgt=avgt,err=sqrt(vart),MSWD=mswd))
}

# x = object of class "logratios"
# ireplace = vector of indices with runs that need replacing
Jmean <- function(nlr,ireplace){
    nruns <- length(nlr)
    idontreplace <- which(!((1:nruns) %in% ireplace))
    nlrmean <- nlr[ireplace[1]]
    nlrout <- c(nlrmean,nlr[-ireplace])
    n <- length(ireplace)
    J <- matrix(0,nrow=sum(nlrout),ncol=sum(nlr))
    for (i in ireplace){ # first calculate the mean
        j <- 1:nlrmean
        k <- getindices.logratios(list(nlr=nlr),i)
        J[j,k] <- diag(nlrmean)/n
    }
    if (length(idontreplace)==0) return(J)
    for (i in 1:length(idontreplace)){
        j <- getindices.logratios(list(nlr=nlrout),i+1)
        k <- getindices.logratios(list(nlr=nlr),idontreplace[i])
        J[j,k] <- diag(nlrout[i+1])
    }
    return(J)
}
