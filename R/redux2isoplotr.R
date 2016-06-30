#' Export \code{ArArRedux} data to \code{IsoplotR}
#'
#' Creates a data object compatible with the \code{IsoplotR} package
#'
#' @param x an object of class \code{\link{redux}}
#' @param irr the irradiation schedule
#' @param fract list with air shot data (one item per denominator
#'     isotope)
#' @param ca an object of class \code{\link{logratios}} with
#'     Ca-interference data (not necessary if interferences are
#'     included in X)
#' @param k an object of class \code{\link{logratios}} with
#'     K-interference data (not necessary if interferences are
#'     included in X)
#' @param format input format for \code{IsoplotR}. I.e. one of
#'
#' \code{1}: 39/40, s[39/40], 36/40, s[36/40], 39/36, s[39/36]
#' 
#' (other formats will be added later)
#' @param file optional (\code{.csv}) file name to write the output
#'     to.
#' @return an object of class \code{ArAr}, i.e. a table with the
#'     following columns: \code{'Ar4036'}, \code{'errAr4036'},
#'     \code{'Ar3936'}, \code{'errAr3936'}, \code{'Ar4039'}, and
#'     \code{'errAr4039'}
#' @examples
#' data(Melbourne)
#' print(redux2isoplotr(Melbourne$X,Melbourne$irr))
#' @export
redux2isoplotr <- function(x,irr,fract=NULL,ca=NULL,
                           k=NULL,format=1,file=NULL){
    X <- redux2ArAr(x,irr,fract=fract,ca=ca,k=k,format=format,file=file)
    if (format==1) out <- format1(X)
    else out <- format2(X)
    if (is.null(file)) return(out)
    else utils::write.table(out,file=file,col.names=FALSE,
                            row.names=FALSE,sep=',')
}

# convert 'redux' data object to 'ArAr' data object
redux2ArAr <- function(x,irr,fract=NULL,ca=NULL,
                       k=NULL,format=1,file=NULL){
    Cl <- corrections(x,irr,fract=fract,ca=ca,k=k)
    Y <- getABCDEF(Cl)
    ni <- length(Y$intercepts)
    out <- list()
    class(out) <- "ArAr"
    R <- get4039(Cl,irr)
    JJ <- getJfactors(R)
    Jfact <- subset(JJ,labels="J:")
    ns <- nruns(Jfact)
    if (ns==1) {
        out$J <- c(Jfact$intercepts,sqrt(Jfact$covmat))
    } else if (diff(range(Jfact$pos))==0) {
        out$J <- c(Jfact$intercepts[1],sqrt(Jfact$covmat[1,1]))
    } else {
        stop("The dataset contains more than one J-factor.")
    }
    labels <- rep(c('Ar36Ar40','Ar39Ar40','Ar39Ar36','Ar40Ar36'),ns)
    out$x <- rep(0,4*ns)
    out$covmat <- matrix(0,4*ns,4*ns)
    J <- matrix(0,nrow=4*ns,ncol=ni)
    hasKglass <- "K-glass" %in% Cl$labels
    hasCasalt <- "Ca-salt" %in% Cl$labels
    for (i in 1:ns){
        ri <- (i-1)*4 # row index (36/40,39/40,39/36,40/36)
        ci <- (i-1)*6 # column index (A,B,C,D,E,F)
        label <- Y$labels[i]
        AA <- Y$intercepts[getindices(Y,label,num='A')]
        CC <- Y$intercepts[getindices(Y,label,num='C')]
        EE <- Y$intercepts[getindices(Y,label,num='E')]
        if (hasKglass) {
            DD <- Y$intercepts[getindices(Y,label,num='D')]
        } else {
            DD <- 0
        }
        if (!hasCasalt | expired(irr[[Y$irr[i]]],Y$thedate[i],Y$param$l7)) {
            BB <- 0
            FF <- 0
        } else {
            BB <- Y$intercepts[getindices(Y,label,num='B')]
            FF <- Y$intercepts[getindices(Y,label,num='F')]
        }
        out$x[ri+1] <- (AA-BB-CC)/(1-DD)
        out$x[ri+2] <- (EE-FF)/(1-DD)
        out$x[ri+3] <- (EE-FF)/(AA-BB-CC)
        out$x[ri+4] <- (1-DD)/(AA-BB-CC)
        J[ri+1,ci+1] <- 1/(1-DD)             # dX1dA
        J[ri+1,ci+2] <- 1/(DD-1)             # dX1dB
        J[ri+1,ci+3] <- 1/(DD-1)             # dX1dC
        J[ri+1,ci+4] <- (AA-BB-CC)/(1-DD)^2  # dX1dD
        J[ri+2,ci+4] <- (EE-FF)/(1-DD)^2     # dY1dD
        J[ri+2,ci+5] <- 1/(1-DD)             # dY1dE
        J[ri+2,ci+6] <- 1/(DD-1)             # dY1dF
        J[ri+3,ci+1] <- (FF-EE)/(AA-BB-CC)^2 # dX2dA
        J[ri+3,ci+2] <- (EE-FF)/(AA-BB-CC)^2 # dX2dB
        J[ri+3,ci+3] <- (EE-FF)/(AA-BB-CC)^2 # dX2dC
        J[ri+3,ci+4] <- 1/(BB+CC-AA)         # dX2dD
        J[ri+4,ci+1] <- (DD-1)/(AA-BB-CC)^2  # dY2dA
        J[ri+4,ci+2] <- (1-DD)/(AA-BB-CC)^2  # dY2dB
        J[ri+4,ci+3] <- (1-DD)/(AA-BB-CC)^2  # dY2dC
        J[ri+4,ci+5] <- 1/(AA-BB-CC)         # dY2dE
        J[ri+4,ci+5] <- 1/(BB+CC-AA)         # dY2dF
    }
    out$covmat <- J %*% Y$covmat %*% t(J)
    names(out$x) <- labels
    rownames(out$covmat) <- labels
    colnames(out$covmat) <- labels
    out
}

format1 <- function(x){
    ns <- length(x$x)/4
    out <- matrix('',ns+3,6)
    out[1,1:2] <- c("J","err[J]")
    out[2,1:2] <- x$J
    out[3,] <- c("Ar39Ar40","errAr39Ar40",
                 "Ar36Ar40","errAr36Ar40",
                 "Ar39Ar36","errAr39Ar36")
    i90 <- findmatches(labels=names(x$x),prefixes=c("Ar39Ar40"))
    i60 <- findmatches(labels=names(x$x),prefixes=c("Ar36Ar40"))
    i96 <- findmatches(labels=names(x$x),prefixes=c("Ar39Ar36"))
    for (i in 1:ns){
        out[i+3,1] <- x$x[i90[i]]
        out[i+3,2] <- sqrt(x$covmat[i90[i],i90[i]])
        out[i+3,3] <- x$x[i60[i]]
        out[i+3,4] <- sqrt(x$covmat[i60[i],i60[i]])
        out[i+3,5] <- x$x[i96[i]]
        out[i+3,6] <- sqrt(x$covmat[i96[i],i96[i]])
    }
    out
}

format2 <- function(x){
    x
}

getJABCDEF <- function(Z,Slabels,nl){
    J <- matrix(0,nrow=nl,ncol=length(Z$intercepts))
    i67ca <- getindices(Z,"Ca-salt","Ar36","Ar37")
    i97ca <- getindices(Z,"Ca-salt","Ar39","Ar37")
    i09k <- getindices(Z,"K-glass","Ar40","Ar39")
    for (i in 1:length(Slabels)){
        j <- (i-1)*6
        label <- Slabels[i]
        i60 <- getindices(Z,label,"Ar36","Ar40")
        i70 <- getindices(Z,label,"Ar37","Ar40")
        i80 <- getindices(Z,label,"Ar38","Ar40")
        i90 <- getindices(Z,label,"Ar39","Ar40")
        i68cl <- getindices(Z,paste("Cl:",label,sep=""),"Ar36","Ar38")
        J[j+1,i60]    <- 1         # A
        J[j+2,c(i67ca,i70)] <- 1   # B
        J[j+3,c(i68cl,i80)] <- 1   # C
        J[j+4,c(i09k,i90)] <- 1    # D
        J[j+5,i90] <- 1            # E
        J[j+6,c(i70,i97ca)] <- 1   # F
    }
    return(J)
}

getABCDEF <- function(Z){
    si <- findrunindices(Z,c("Ca-salt","K-glass","Cl:"),invert=TRUE)
    ns <- length(si)
    out <- subset(Z,si)
    out$num <- c(rep(c("A","B","C","D","E","F"),ns))
    out$den <- rep(NA,6*ns)
    out$nlr <- rep(6,ns)
    Jv <- getJABCDEF(Z,out$labels,6*ns)
    out$intercepts <- exp(Jv %*% Z$intercepts)
    Jw <- apply(Jv,2,"*",out$intercepts)
    out$covmat <- Jw %*% Z$covmat %*% t(Jw)
    return(out)
}
