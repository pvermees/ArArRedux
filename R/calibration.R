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
calibration <- function(x,...){ UseMethod("calibration",x) }
#' @rdname calibration
#' @export
calibration.default <- function(x,clabel,...){
    j <- grep(clabel,x$labels,invert=TRUE)
    if (length(j)==length(x$labels)) return(X)
    J <- Jcal(x,clabel,x$detectors)
    out <- subset(x,j)
    out$intercepts <- J %*% x$intercepts
    out$covmat <- J %*% x$covmat %*% t(J)
    return(out)
}
#' @rdname calibration
#' @export
calibration.WiscAr <- function(x,irradiations,clabel){
    D <- getDmatrix(x,irradiations,clabel)
    X <- concat(list(x,D))
    icalgas <- findrunindices(X,prefixes=clabel)
    iothers <- findrunindices(X,prefixes=clabel,invert=TRUE)
    ng <- length(icalgas)/3
    ns <- length(iothers)/3
    icalgas <- icalgas[1:ng]
    iothers <- iothers[1:ns]
    nc <- length(X$intercepts)
    J <- matrix(0,4*ns,nc)
    # standard gas indices
    i09g101 <- getindices(X,prefix=clabel,num='Ar40',den='Ar39',hop='101')
    i79g101 <- getindices(X,prefix=clabel,num='Ar37',den='Ar39',hop='101')
    i69g101 <- getindices(X,prefix=clabel,num='Ar36',den='Ar39',hop='101')
    i89g102 <- getindices(X,prefix=clabel,num='Ar38',den='Ar39',hop='102')
    i69g102 <- getindices(X,prefix=clabel,num='Ar36',den='Ar39',hop='102')
    # all measurement indices
    i09a101 <- getindices(X,num='Ar40',den='Ar39',hop='101')
    i79a101 <- getindices(X,num='Ar37',den='Ar39',hop='101')
    i69a101 <- getindices(X,num='Ar36',den='Ar39',hop='101')
    i89a102 <- getindices(X,num='Ar38',den='Ar39',hop='102')
    i69a102 <- getindices(X,num='Ar36',den='Ar39',hop='102')
    # sample indices are difference between 'g' and 'a'
    i09s101 <- i09a101[!(i09a101 %in% i09g101)]
    i79s101 <- i79a101[!(i79a101 %in% i79g101)]
    i69s101 <- i69a101[!(i69a101 %in% i69g101)]
    i89s102 <- i89a102[!(i89a102 %in% i89g102)]
    i69s102 <- i69a102[!(i69a102 %in% i69g102)]
    # sample decay correction indices
    iD7 <- getindices(X,prefix='DCAL',num='Ar37',den='NA')
    iD9 <- getindices(X,prefix='DCAL',num='Ar39',den='NA')
    # true standard gas composition indices
    i09w <- getindices(X,prefix='DCAL',num='Ar40',den='Ar39')
    i89w <- getindices(X,prefix='DCAL',num='Ar38',den='Ar39')
    i69w <- getindices(X,prefix='DCAL',num='Ar36',den='Ar39')
    # match the samples with the nearest standard gas
    inearestgas <- nearest(X$thedate[iothers],X$thedate[icalgas])
    for (i in 1:ns){
        j <- inearestgas[i]
        # s09corr
        J[(4*i-3),i09g101[j]] <- -1
        J[(4*i-3),i09s101[i]] <-  1
        J[(4*i-3),iD9[i]]     <- -1
        J[(4*i-3),i09w[j]]    <-  1
        # s89corr
        J[(4*i-2),i89g102[j]] <- -1
        J[(4*i-2),i89s102[i]] <-  1
        J[(4*i-2),iD9[i]]     <- -1
        J[(4*i-2),i89w[j]]    <-  1
        # s79corr
        J[(4*i-1),i79s101[i]] <-  1
        J[(4*i-1),i89g102[j]] <- -1
        J[(4*i-1),i69g102[j]] <-  1
        J[(4*i-1),iD7[i]]     <-  1
        J[(4*i-1),iD9[i]]     <- -1
        J[(4*i-1),i89w[j]]    <-  1
        J[(4*i-1),i69w[j]]    <- -1
        # s69corr
        J[4*i,i69g101[j]] <- -1/2
        J[4*i,i69s101[i]] <-  1/2
        J[4*i,i69g102[j]] <- -1/2
        J[4*i,i69s102[i]] <-  1/2
        J[4*i,iD9[i]]     <- -1
        J[4*i,i69w[j]]    <-  1
    }
    isamp <- findrunindices(x,prefixes=clabel,invert=TRUE)[1:ns]
    out <- subset(X,i=isamp)
    out$intercepts <- J %*% X$intercepts
    out$covmat <- J %*% X$covmat %*% t(J)
    out$hop <- NULL
    out$num <- rep(c("Ar40","Ar38","Ar37","Ar36"),ns)
    out$den <- rep("Ar39",4*ns)
    out$nlr <- rep(4,ns)
    out
}

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
