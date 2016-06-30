#' Apply the mass fractionation correction
#'
#' Applies the fractionation obtained from air shot data by
#' \code{\link{fractionation}} to the denominator detector in order to
#' correct it for the mass difference between the numerator and
#' denominator isotopes.
#' 
#' @param X an object of class \code{redux}
#' @param fract a list with fractionation data for the detectors used
#' to measure the denominator isotopes
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' A <- massfractionation(C,Melbourne$fract)
#' plotcorr(A)
#' @export
massfractionation <- function(X,fract){
    fdet <- X$detectors[c("Ar37","Ar39","Ar40")]
    # add the air shot data
    errorlessair <- air(X)
    errorlessair$covmat <- 0
    if (methods::is(fract,"logratios")){
        Y <- concat(list(X,fract,errorlessair))
    } else {
        Y <- concat(c(list(X),fract,list(errorlessair)))
    }
    # apply the fractionation correction
    J <- Jair(Y,fdet)
    j <- findrunindices(Y,c(unlist(fdet),"air-ratio"),invert=TRUE)
    out <- subset(Y,j)
    out$intercepts <- J %*% Y$intercepts
    out$covmat <- J %*% Y$covmat %*% t(J)
    return(out)
}

#' Compute the mass fractionation correction
#'
#' Compares the measured 40Ar/36Ar ratio of an air shot on a given
#' detector with the atmospheric ratio.
#' 
#' @param fname a .csv file with the air shot data
#' @param detector the name of the ion detector
#' @param MS the type of mass spectrometer
#' @param PH TRUE if the data were recorded in 'peak hopping' mode,
#' FALSE if recorded in multicollector mode.
#' @examples
#' data(Melbourne)
#' fd37file <- system.file("AirL2.csv",package="ArArRedux")
#' fd40file <- system.file("AirH1.csv",package="ArArRedux")
#' fract <- list(fractionation(fd37file,"L2",PH=TRUE),
#'               fractionation(fd40file,"H1",PH=FALSE))
#' if (isTRUE(all.equal(Melbourne$fract,fract))){
#'   print("We just re-created the fractionation correction for the Melbourne dataset")
#' }
#' @return an object of class \code{\link{logratios}}
#' @export
fractionation <- function(fname,detector,MS="ARGUS-VI",PH=FALSE){
    mf <- loaddata(fname,c("Ar40","Ar36"),MS,PH)
    lf <- fitlogratios(blankcorr(mf),"Ar40")
    f <- averagebyday(lf,detector)
    return(f)
}

Jair <- function(X,detectors){
    ns <- nruns(X)
    di <- findrunindices(X,detectors)
    ai <- findrunindices(X,"air-ratio")
    dai <- c(di,ai)
    si <- (1:ns)[-dai] # sample indices
    ndai <- sum(X$nlr[dai])
    nsi <- sum(X$nlr[si])
    J <- matrix(0,nsi,nsi+ndai)
    for (i in si){ # loop through the samples
        # intercept indices of current sample
        sj <- getindices(X,prefix=X$labels[i])
        for (js in sj){ # loop through all the masses
            num <- X$num[js]
            den <- X$den[js]
            a <- as.integer(substr(num,start=3,stop=6))
            b <- as.integer(substr(den,start=3,stop=6))
            # indices of the nearest calibration data
            dj <- array(grep(detectors[[den]],X$labels))
            id <- di[nearest(X$thedate[i],X$thedate[dj])]
            jd <- getindices(X,prefix=X$labels[id])
            ja <- getindices(X,prefix="air-ratio")
            J[js,js] <- 1
            J[js,jd] <- (log(a)-log(b))/(log(40)-log(36))
            J[js,ja] <- (log(a)-log(b))/(log(40)-log(36))
        }
    }
    return(J)
}

air <- function(X){
    out <- list(
        num = "Ar40",
        den = "Ar36",
        intercepts = log(X$param$air), 
        covmat = (X$param$sair/X$param$air)^2, # variance of the air ratio
        irr = NULL,
        pos = NULL,
        labels = "air-ratio",
        thedate = as.numeric(as.Date("1970-01-01 00:00:00")),
        nlr = 1
    )
    return(out)
}
