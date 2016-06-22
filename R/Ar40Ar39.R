#' Calculate the 40Ar*/39ArK-ratios
#'
#' Calculate the 40Ar*/39ArK-ratios of interference corrected logratio
#' intercept data
#' 
#' @param X an object of class \code{redux} containing some
#' interference corrected logratio intercept data
#' @param irr the irradiation schedule
#' @return an object of class \code{link{redux}} containing the
#' 40Ar*/39ArK-ratios as \code{intercepts} and its covariance matrix
#' as \code{covmat}
#' @examples
#' data(Melbourne)
#' R <- get4039(Melbourne$X,Melbourne$irr)
#' plotcorr(R)
#' @export
get4039 <- function(X,irr){
    Y <- getabcdef(X)
    ns <- nruns(Y)-1
    ni <- length(Y$intercepts)
    out <- subset(Y,1:ns)
    out$intercepts <- rep(0,ns)
    out$covmat <- matrix(0,ns,ns)
    out$num <- rep(NA,ns)
    out$den <- rep(NA,ns)
    out$nlr <- rep(1,ns)
    J <- matrix(0,nrow=ns,ncol=ni)
    hasKglass <- "K-glass" %in% X$labels
    hasCasalt <- "Ca-salt" %in% X$labels
    for (i in 1:ns){
        j <- (i-1)*6
        label <- Y$labels[i]
        aa <- Y$intercepts[getindices(Y,label,num='a')]
        bb <- Y$intercepts[getindices(Y,label,num='b')]
        cc <- Y$intercepts[getindices(Y,label,num='c')]
        dd <- Y$intercepts[getindices(Y,label,num='d')]
        ee <- Y$intercepts[getindices(Y,label,num='e')]
        if (hasKglass){
            ff <- Y$intercepts[getindices(Y,label,num='f')]
            J[i,j+6] <- -1
        } else {
            ff <- 0
        }
        if (!hasCasalt | expired(irr[[Y$irr[i]]],Y$thedate[i],Y$param$l7)){
            out$intercepts[i] <- (1-aa+cc)/dd-ff
            J[i,j+1] <- -1/dd           # dR/da
            J[i,j+2] <-  0              # dR/db
            J[i,j+3] <-  1/dd           # dR/dc
            J[i,j+4] <- -(1-aa+cc)/dd^2 # dR/dd
            J[i,j+5] <-  0              # dR/de
        } else {
            out$intercepts[i] <- (1-aa+bb+cc)/(dd-ee)-ff
            J[i,j+1] <- -1/(dd-ee)
            J[i,j+2] <-  1/(dd-ee)
            J[i,j+3] <-  1/(dd-ee)
            J[i,j+4] <- -(1-aa+bb+cc)/(dd-ee)^2
            J[i,j+5] <-  (1-aa+bb+cc)/(dd-ee)^2
        }
    }
    out$covmat <- J %*% Y$covmat %*% t(J)
    return(out)
}

getJabcdef <- function(Z,Slabels,nl){
    J <- matrix(0,nrow=nl,ncol=length(Z$intercepts))
    i67ca <- getindices(Z,"Ca-salt","Ar36","Ar37")
    i97ca <- getindices(Z,"Ca-salt","Ar39","Ar37")
    i09k <- getindices(Z,"K-glass","Ar40","Ar39")
    for (i in 1:length(Slabels)){
        j <- (i-1)*6
        label <- Slabels[i]
        iair <- getindices(Z,"air-ratio","Ar40","Ar36")
        i60 <- getindices(Z,label,"Ar36","Ar40")
        i70 <- getindices(Z,label,"Ar37","Ar40")
        i80 <- getindices(Z,label,"Ar38","Ar40")
        i90 <- getindices(Z,label,"Ar39","Ar40")
        i68cl <- getindices(Z,paste("Cl:",label,sep=""),"Ar36","Ar38")
        J[j+1,c(iair,i60)] <- 1
        J[j+2,c(iair,i67ca,i70)] <- 1
        J[j+3,c(iair,i68cl,i80)] <- 1
        J[j+4,i90] <- 1
        J[j+5,c(i97ca,i70)] <- 1
        J[j+6,i09k] <- 1
    }
    return(J)
}

getabcdef <- function(Cl){
    Z <- concat(list(Cl,air(Cl))) # matrix with everything
    si <- getrunindices(Z,c("Ca-salt","K-glass","Cl:","air-ratio"),
                        invert=TRUE)
    ns <- length(si)
    theS <- subset(Z,si) # contains only samples
    theK <- subset(Z,labels="K-glass") # contains only K-glass
    out <- Z
    out$irr <- c(theS$irr,theK$irr)
    out$pos <- c(theS$pos,theK$pos)
    out$labels <- c(theS$labels,theK$labels)
    out$num <- c(rep(c("a","b","c","d","e","f"),ns))
    out$den <- rep(NA,6*ns)
    out$nlr <- rep(6,ns)
    out$thedate <- c(theS$thedate,theK$thedate)
    Jv <- getJabcdef(Z,theS$labels,length(out$num))
    out$intercepts <- exp(Jv %*% Z$intercepts)
    Jw <- apply(Jv,2,"*",out$intercepts)
    out$covmat <- Jw %*% Z$covmat %*% t(Jw)
    return(out)
}
