# X = object of class "redux"
Jcal <- function(X,clabel,detectors){
    ci <- array(grep(clabel, X$labels)) # calibration indices
    si <- array(grep(clabel, X$labels, invert=TRUE))
    nc <- length(ci) # number of calibrations
    ns <- length(si)
    nci <- sum(X$nlr[ci]) # number of calibration intercepts
    nsi <- sum(X$nlr[-ci]) # number of sample intercepts
    J <- matrix(0,nrow=nsi,ncol=nsi+nci)
    for (i in si){ # loop through all the samples
        # indices of the nearest calibration data
        ic <- ci[nearest(X$thedate[i],X$thedate[ci])]
        # intercept indices of current sample
        sj <- getindices(X,prefix=X$labels[i])
        for (js in sj){ # loop through all the masses
            num <- X$num[js]
            den <- X$den[js]
            # get the detectors for the numerator and denominator masses
            ndet <- detectors[[num]]
            ddet <- detectors[[den]]
            if (ndet != ddet) {
                jn <- getindices(X,prefix=X$labels[ic],num=ndet)
                jd <- getindices(X,prefix=X$labels[ic],num=ddet)
                J[js,jn] <- -1
                J[js,jd] <- 1
            }
            J[js,js] <- 1
        }
    }
    return(J)
}

#' Detector calibration
#'
#' Apply the detector calibration for multicollector data
#' 
#' @param X an object of class \code{redux}
#' @param clabel the label of the detector calibration data
#' @return an object of class \code{redux}
#' @examples
#' data(Melbourne)
#' C <- calibration(Melbourne$X,"DCAL")
#' plotcorr(C)
#' @export
calibration <- function(X,clabel){
    j <- grep(clabel,X$labels,invert=TRUE)
    if (length(j)==length(X$labels)) return(X)
    J <- Jcal(X,clabel,X$detectors)
    out <- subset(X,j)
    out$intercepts <- J %*% X$intercepts
    out$covmat <- J %*% X$covmat %*% t(J)
    return(out)
}
