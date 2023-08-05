#' Extrapolation to 'time zero'
#'
#' This function extrapolates time resolved mass spectrometer data to
#' t=0. When fed with multicollector data, it forms the ratios of the
#' raw signals, forms their logs and performs linear regression to t=0
#' When fed with single collector data, the function first takes their
#' logs and extrapolates them to t=0 before taking ratios, unless
#' \code{denmass}=NULL, in which case the logs of the raw signals are
#' extrapolated.
#' 
#' @param x an object of class \code{timeresolved} or \code{PHdata}
#' @param ... further arguments (see below)
#' @return an object of class \code{logratios}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' blanklabel <- "EXB#"
#' l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
#' plotcorr(l)
#' @export
fitlogratios <- function(x,...){ UseMethod("fitlogratios",x) }
#' @rdname fitlogratios
#' @export
fitlogratios.default <- function(x,...){ stop() }
#' @param denmass a string denoting the denominator isotope
#' @rdname fitlogratios
#' @export
fitlogratios.timeresolved <- function(x,denmass="Ar40",tmin=0,...){
    r <- takeratios(x,denmass)
    l <- takelogs(r)
    f <- timezero(l,tmin=tmin)
    return(cast(f,"logratios"))
}
#' @rdname fitlogratios
#' @export
fitlogratios.PHdata <- function(x,denmass=NULL,tmin=0,...){
    for (i in 1:nmasses(x)){
        z <- newfit(x)
        l <- takelogs(x$signals[[i]])
        z <- setmasses(x=z,mass=z$masses[i],value=timezero(l,tmin=tmin))
    }
    if (is.null(denmass)) {
        f <- z
        f$denmass <- NA
    } else {
        f <- takeratios(z,denmass)
    }
    return(cast(f,"logratios"))
}
#' @rdname fitlogratios
#' @export
fitlogratios.WiscAr <- function(x,denmass='Ar39',...){
    out <- x
    for (hop in names(x)){
        out[[hop]] <- fitlogratios(x[[hop]],denmass=denmass,...)
    }
    out
}

Jtakeratios <- function(nruns,inum,iden){
    nlr <- length(inum)
    nmasses <- nlr+1
    J <- matrix(0,nrow=nlr,ncol=nmasses) # elementary matrix
    for (i in 1:length(inum)){
        J[i,inum[i]] <- 1
        J[,iden] <- -1
    }
    out <- matrix(0,nrow=nlr*nruns,ncol=nmasses*nruns)
    for (i in 1:nruns){
        irow <- ((i-1)*nlr+1):(i*nlr)
        icol <- ((i-1)*nmasses+1):(i*nmasses)
        out[irow,icol] <- J
    }
    return(out)
}

takeratios <- function(x,...){ UseMethod("takeratios",x) }
takeratios.default <- function(x,...){ stop() }
takeratios.timeresolved <- function(x,denmass,...){
    den <- getmasses(x,denmass)
    num <- getmasses(x,denmass,invert=TRUE)
    out <- num
    for (mass in num$masses){ # loop through the isotopes
        ratio <- getmasses(num,mass)$d/den$d
        out <- setmasses(out,mass,ratio)
    }
    out$denmass <- denmass
    class(out) <- append(class(out),"ratio")
    return(out)
}
takeratios.fit <- function(x,denmass,...){
    out <- x
    iden <- which(x$masses == denmass)
    inum <- which(x$masses != denmass)
    J <- Jtakeratios(nruns(x),inum,iden)
    out$masses <- x$masses[inum]
    out$intercepts <- J %*% x$intercepts
    out$covmat <- J %*% x$covmat %*% t(J)
    out$denmass <- denmass
    return(out)
}

takelogs <- function(x){
    out <- x
    out <- replacenegatives(x)
    out$d <- log(out$d)
    class(out) <- append(class(out),"logged")
    return(out)
}

newfit <- function(x,...){ UseMethod("newfit",x) }
newfit.default <- function(x,nmasses=NULL,nruns=NULL,...){
    out <- x
    class(out) <- "fit"
    if (is.null(nmasses)) nmasses <- nmasses(x)
    if (is.null(nruns)) nruns <- nruns(x)
    out$d <- NULL
    out$thetime <- NULL
    out$intercepts <- rep(0,nmasses*nruns)
    out$covmat <- matrix(0,nrow=nmasses*nruns,ncol=nmasses*nruns)
    return(out)
}
newfit.PHdata <- function(x,...){
    x1 <- x$signals[[1]]
    out <- newfit(x1,nmasses(x),nruns(x1))
    out$masses <- x$masses
    out$signals <- NULL
    return(out)
}

timezero <- function(x,blankindices=NA,tmin=0){
    if (all(is.na(blankindices)))
        blankindices <- x$blankindices
    nmasses <- nmasses(x)
    out <- newfit(x,nmasses)
    bi <- rle(as.vector(blankindices))$values # blank indices
    irunsout <- 0
    for (i in bi){ # loop through the groups
        irunsx <- which(i==blankindices)
        irunsout <- utils::tail(irunsout,n=1) + (1:length(irunsx))
        g <- subset(x,irunsx) # extract group
        out <- setfit(out,fit(g,tmin=tmin),nmasses,irunsout)
    }
    return(out)
}

# x and f are both objects of class "fit"
# nmasses = number of masses or isotope ratios per sample
# iruns = indicates which samples need replacing
setfit <- function(x,f,nmasses,iruns){
    out <- x
    ii <- getindices(nmasses=nmasses,iruns=iruns)
    out$intercepts[ii] <- f$intercepts
    out$covmat[ii,ii] <- f$covmat
    return(out)
}

fit <- function(x,tmin=0,mass=NA,returnfit=FALSE){
    out <- newfit(x)
    if (nruns(x)>1) { # average time for all samples in group
        thetime <- apply(x$thetime,1,"mean")
    } else {
        thetime <- x$thetime
    }
    i <- which(thetime>tmin)
    if (is.na(mass)){
        dat <- x$d[i,]
    } else {
        m <- which(x$masses %in% mass)
        dat <- x$d[i,m]
    }
    tt <- thetime[i]
    f <- stats::lm(dat ~ tt)# + I(tt^2))
    if (nruns(x) == 1 & (nmasses(x)==1 | !is.na(mass))) { # only one intercept
        out$intercepts <- stats::coef(f)["(Intercept)"]
        covcolname <- "(Intercept)"
    } else { # an entire row of intercepts
        out$intercepts <- stats::coef(f)["(Intercept)",]
        covcolname <- ":(Intercept)"
    }
    myvcov <- stats::vcov(f)
    j <- which(colnames(myvcov)==covcolname)
    out$covmat <- myvcov[j,j]
    if (returnfit) out$fit <- f
    return(out)
}
