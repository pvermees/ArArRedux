#' Create a new \code{\link{redux}} object
#'
#' Initialises a new \code{\link{redux}} object by packing a
#' \code{\link{logratios}} dataset together with all the parameters
#' needed for age calculation
#' 
#' @param X an object of class \code{\link{logratios}}
#' @param Jpos a vector of integers denoting the positions of the
#' fluence monitors in the irradiation stack
#' @param detectors a list of strings denoting the detectors for each
#' argon isotope
#' @return an object of class \code{\link{redux}}
#' @export
newredux <- function(X,Jpos,detectors=
            list(Ar36="H1",Ar37="L2",Ar38="L1",Ar39="AX",Ar40="H1")){
    out <- X
    out$Jpos <- Jpos
    out$detectors <- detectors
    out$param$l0 = 5.5492e-4 # 40K decay constant [Ma-1]
    out$param$sl0 = 0.0047e-4 # Renne et al. (2010)
    out$param$l7 = 365.25*log(2)/34.95 # 37Ar decay constant [yr-1]
    out$param$sl7 = out$param$l7*0.04/34.95 # Renne and Norman (2001)
    out$param$l9 = log(2)/269 # 39Ar decay constant [a-1]
    out$param$sl9 = out$param$l9*1.5/269 # Stoenner et al. (1965)
    out$param$l6 = log(2)/301200 # 36Cl decay constant [a-1]
    out$param$sl6 = out$param$l6*1000/301200 # uncertainty
    out$param$pcl = 252.7 # P(36Cl/38Cl) for OSTR reactor
    out$param$spcl = 1.8 # Renne et al. (2008)
    out$param$ts = 28.201 # age of the Fish Canyon Tuff
    out$param$sts = 0.023 # 1-sigma age uncertainty
    out$param$air = 298.56 # atmospheric 40/36-ratio
    out$param$sair = 0.155 # Lee et al. (2006)
    class(out) <- "redux"
    return(out)
}

cast <- function(x,...){ UseMethod("cast",x) }
cast.default <- function(x,...){stop()}
cast <- function(from,to){
    out <- from
    class(out) <- to
    if (to == "logratios"){
        nruns <- nruns(out)
        num <- from$masses
        den <- rep(from$denmass,nmasses(from))
        out$num <- rep(num,nruns)
        out$den <- rep(den,nruns)
        out$masses <- NULL
        out$denmass <- NULL
        out$blankindices <- NULL
        out$nlr <- rep(nmasses(from),nruns)
    }
    return(out)
}

readthedate <- function(x){
    thedate <- strptime(x,"%d-%b-%Y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d/%m/%y %H:%M")
    if (any(is.na(thedate))) thedate <- strptime(x,"%d-%b-%Y")
    return(as.numeric(thedate))
}

# masses = vector of strings
# cirr = column with irradiation can label
# cpos = column with position in irradiation stack
# clabel = column containing the sample labels
# cdate = column containing the date
# ci = matrix with column indices of the masses of interest
newtimeresolved <- function(thetable,masses,cirr,cpos,clabel,cdate,ci){
    out <- list()
    class(out) <- "timeresolved"
    out$masses <- masses
    out$irr <- as.character(thetable[,cirr])
    out$pos <- as.numeric(thetable[,cpos])
    out$labels <- as.character(thetable[,clabel])
    out$thedate <- readthedate(thetable[,cdate])
    nmasses <- nmasses(out)
    nruns <- nruns(out)
    ncycles <- dim(ci)[1]
    out$d <- matrix(0,nrow=ncycles,ncol=nmasses*nruns)
    out$thetime <- t(thetable[,ci[,1]+1])
    for (i in 1:ncycles){
        for (j in 1:nmasses){
            k <- seq(from=j,to=nmasses*nruns,by=nmasses)
            out$d[i,k] <- thetable[,ci[i,j]]
        }
    }
    return(out)
}

newPHdata <- function(thetable,masses,cirr,cpos,clabel,cdate,ci){
    out <- list()
    class(out) <- "PHdata"
    out$masses <- masses
    for (i in 1:nmasses(out)){
        out$signals[[masses[i]]] <-
            newtimeresolved(thetable,masses[i],
                cirr,cpos,clabel,cdate,as.matrix(ci[,i]))
    }
    return(out)
}

#' Load mass spectrometer data
#'
#' Loads a .csv file with raw mass spectrometer data
#' 
#' @param fname the file name, must end with .csv
#' @param masses a vector of strings denoting the order of the
#' isotopes listed in the table
#' @param MS the type of mass spectrometer
#' @param PH a boolean indicating whether the data are to be treated
#' as multicollector (PH=FALSE) or 'peak hopping' (PH=TRUE) data. The
#' default is PH=FALSE.
#' @return if PH=FALSE: an object of class \code{timeresolved}\cr
#' if PH=TRUE: an object of class \code{PHdata}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' plot(m,"MD2-1a","Ar40")
#' @export
loaddata <- function(fname,masses,MS="ARGUS-VI",PH=FALSE){
    thetable <- utils::read.csv(file=fname,header=FALSE,skip=3)
    nrows <- dim(thetable)[1] # number of MS runs
    ncols <- dim(thetable)[2] # number of columns
    nmass <- length(masses) # number of isotopes
    ncycles <- (ncols-3)/(2*nmass)
    cirr <- 1
    cpos <- 2
    clabel <- 3
    cdate <- 4 # column with the date
    ci <- matrix(NA,nrow=ncycles,ncol=nmass) # column indices
    for (i in 1:nmass){
        ci[,i] <- seq(from=3+2*i,to=ncols,by=2*nmass)
    }
    if (PH) {
        out <- newPHdata(thetable,masses,cirr,cpos,clabel,cdate,ci)
    } else {
        out <- newtimeresolved(thetable,masses,cirr,cpos,clabel,cdate,ci)
    }
    return(out)
}

#' Load the irradiation schedule
#'
#' Loads a .csv file with the schedule of a multi-stage neutron
#' irradiation
#' 
#' @param fname file name (in .csv format)
#' @return a list of irradiations, where each irradiation is a named
#' list containing:
#' 
#' \code{tin}: vector with the start times of irradiations \cr
#' \code{tout}: vector with the end times of irradiations \cr
#' \code{P}: vector with the power of the irradiations
#' @examples
#' irrfile <- system.file("irradiations.csv",package="ArArRedux")
#' irr <- loadirradiations(irrfile)
#' str(irr)
#' @export
loadirradiations <- function(fname){
    f <- file(fname)
    open(f)
    out <- list()
    while (TRUE) { # read the file line by line
        line <- readLines(f, n=1, warn=FALSE)
        if (length(line) <= 0){ # EOF
            close(f)
            return(out)
        }
        l <- unlist(strsplit(line, split=','))
        if (l[1]=='In') {
            out[[irrname]]$tin <- c(out[[irrname]]$tin,readthedate(l[2]))
            out[[irrname]]$P <- c(out[[irrname]]$P,as.numeric(l[3]))
        } else if (l[1]=='Out'){
            out[[irrname]]$tout <- c(out[[irrname]]$tout,readthedate(l[2]))
        } else {
            irrname <- l[1]
            out[[irrname]] <- list(P=c(),tin=c(),tout=c())
        }
    }
}

#' Read mass spectrometer data
#'
#' Reads raw mass spectrometer data and parses it into a
#' \code{\link{redux}} format for further processing.
#' 
#' @param xfile a .csv file with samples and fluence monitor data
#' @param masses a list which specifies the order in which the isotopes
#' are recorded by the mass spectrometer
#' @param blabel a prefix string denoting the blanks
#' @param Jpos a vector of integers denoting the positions of the
#' fluence monitors in the irradiation stack
#' @param kfile a .csv file with the K-interference monitor data
#' (optional)
#' @param cafile a .csv file with the Ca-interference monitor data
#' (optional)
#' @param dfile a .csv file with the detector calibration data
#' (optional)
#' @param dlabels a list which specifies the names of the detectors
#' and the order in which they were recorded by the mass spectrometer
#' @param MS a string denoting the type of mass spectrometer. At the
#' moment only the ARGUS-IV is listed. Please contact the author to
#' add other file formats to Ar-Ar_Redux.
#' @return an object of class \code{\link{redux}}.
#' @examples
#' samplefile <-  system.file("Samples.csv",package="ArArRedux")
#' kfile <- system.file("K-glass.csv",package="ArArRedux")
#' cafile <- system.file("Ca-salt.csv",package="ArArRedux")
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' dlabels <- c("H1","AX","L1","L2")
#' X <- read(samplefile,masses,"EXB#",c(3,15),kfile,cafile,dfile,dlabels)
#' plotcorr(X)
#' @export
read <- function(xfile,masses,blabel,Jpos,kfile=NULL,cafile=NULL,
                  dfile=NULL,dlabels=NULL,MS="ARGUS-VI"){
    # load the .csv files
    m <- loaddata(xfile,masses,MS) # samples and J-standards
    x <- fitlogratios(blankcorr(m,blabel),"Ar40")    
    if (!is.null(kfile)){ # K-interference data
        # subset of the relevant isotopes
        mk <- loaddata(kfile,masses,MS)
        lk <- fitlogratios(blankcorr(mk,blabel,"K:"),"Ar39")
        k <- getmasses(lk,"Ar40","Ar39")
        x <- concat(list(x,k))
    }
    if (!is.null(cafile)){ # Ca interference data
        mca <- loaddata(cafile,masses,MS)
        lca <- fitlogratios(blankcorr(mca,blabel,"Ca:"),"Ar37")
        ca <- getmasses(lca,c("Ar36","Ar39"),c("Ar37","Ar37"))
        x <- concat(list(x,ca))
    }
    if (!is.null(dfile)){
        md <- loaddata(dfile,dlabels,MS,PH=TRUE)
        ld <- fitlogratios(blankcorr(md))
        d <- averagebyday(ld,"DCAL")
        x <- concat(list(x,d))
    }
    return(newredux(x,Jpos))
}

#' Set or get Ar-Ar_Redux parameters
#'
#' This function is used to query and modify the half lives, standard
#' ages etc. associated with an object of class \code{\link{redux}}
#'
#' \code{\link{param}} grants access to the following parameters:
#'  
#' \code{l0}: 40K decay constant (default value = 5.5492e-4 Ma-1,
#' Renne et al. [2010])\cr
#' \code{sl0}: standard error of the 40K decay constant (default value
#' = 0.0047e-4 Ma-1)\cr
#' \code{l7}: 37Ar decay constant (default value = 7.2438 yr-1, Renne
#' and Norman [2001])\cr
#' \code{sl7}: standard error of the 37Ar decay constant (default
#' value = 0.0083 yr-1)\cr
#' \code{l9}: 39Ar decay constant (0.002577 yr-1 Stoenner et
#' al. [1965])\cr
#' \code{sl9}: standard error of the 39Ar decay constant (0.000014
#' yr-1)\cr
#' \code{l6}: 36Cl decay constant (default value = 2301.3e-9 yr-1)\cr
#' \code{sl6}: standard error of the 36Cl decay constant (default
#' value = 7.6e-9 yr-1\cr
#' \code{pcl}: (36Cl/38Cl)-production rate (default value = 252.7 for
#' OSTR reactor, Renne et al. [2008])\cr
#' \code{spcl}: standard error of the (36Cl/38Cl)-production rate
#' (default value = 1.8)\cr
#' \code{ts}: age of the fluence monitor (default = 28.201 Myr for the
#' Fish Canyon Tuff, Kuiper et al. [2008])\cr
#' \code{sts}: standard error of the fluence monitor age (default
#' value = 0.023 Myr)\cr
#' \code{air}: atmospheric 40Ar/36Ar ratio (default value = 298.56,
#' Lee et al. [2006])\cr
#' \code{sair}: standard error of the atmospheric 40Ar/36Ar ratio
#' (default value = 0.155)
#' 
#' @param X an object of class \code{\link{redux}}
#' @param ... any combination of the parameters given below
#' @return returns the modified \code{\link{redux}} object OR the
#' current parameter values if no optional arguments are supplied.
#' @examples
#' data(Melbourne)
#' param(Melbourne$X)$air
#' Y <- param(Melbourne$X,air=295.5)
#' param(Y)$air
#' @export
param <- function(X,...){
    arguments <- list(...)
    if (length(arguments)>0){
        X$param[names(arguments)] <- arguments
        return(X)
    } else {
        return(X$param)
    }
}

#' Export \code{ArArRedux} data to \code{IsoplotR}
#'
#' Creates a data object compatible with the \code{IsoplotR} package
#'
#' @param x an object of class \code{\link{redux}}
#' @return an object of class \code{ArAr}, i.e. a table with the
#'     following columns: \code{'Ar4036'}, \code{'errAr4036'},
#'     \code{'Ar3936'}, \code{'errAr3936'}, \code{'Ar4039'}, and
#'     \code{'errAr4039'}
#' @examples
#' data(Melbourne)
#' print(redux2isoplotr(Melbourne$X))
redux2isoplotr <- function(X,irr,fract=NULL,ca=NULL,k=NULL,format=1){
    Cl <- corrections(X,irr,fract=fract,ca=ca,k=k)
    # calculate the 40Ar*/39ArK-ratios 
    R <- get4039(Cl,irr)
    # calculate J factors
    J <- getJfactors(R)
    i <- findmatches(Cl$pos,prefixes=Cl$Jpos,invert=TRUE)
    Y <- subset(Cl,i)
    if (format==1){
        out <- redux2isoplotr1(Y,J)
    }
    out
}

redux2isoplotr1 <- function(x,J){
    ss <- getmasses(x,c('Ar36','Ar39'),c('Ar40')) # subset
    snames <- ss$labels
    nn <- length(snames)
    out <- list()
    out$format <- 1
    out$x <- matrix(0,nn,6)
    colnames(out$x) <- c('Ar39Ar40','errAr39Ar40','Ar36Ar40',
                       'errAr36Ar40','Ar39Ar36','errAr39Ar36')
    rownames(out$x) <- snames
    Jisoplotr <- matrix(0,3*nn,2*nn)
    for (ii in 1:nn){
        if (ss$nlr[ii]==2){
            l3640 <- getsignal(ss,snames[ii],num='Ar36')[1]
            l3940 <- getsignal(ss,snames[ii],num='Ar39')[1]
            Ar3940 <- exp(l3940)
            Ar3640 <- exp(l3640)
            Ar3936 <- exp(l3940-l3640)
            dAr3940dl3940 <-  Ar3940
            dAr3640dl3640 <-  Ar3640
            dAr3936dl3940 <-  Ar3936
            dAr3936dl3640 <- -Ar3936
            j <- 3*(ii-1)+1
            k <- 2*(ii-1)+1
            out$x[ii,'Ar39Ar40'] <- Ar3940
            out$x[ii,'Ar36Ar40'] <- Ar3640
            out$x[ii,'Ar39Ar36'] <- Ar3936
            Jisoplotr[j,k]     <- dAr3940dl3940
            Jisoplotr[j+1,k+1] <- dAr3640dl3640
            Jisoplotr[j+2,k]   <- dAr3936dl3940
            Jisoplotr[j+2,k+1] <- dAr3936dl3640
        }
    }
    covmat <- Jisoplotr %*% ss$covmat %*% t(Jisoplotr)
    err <- sqrt(diag(covmat))
    l <- seq(from=1,to=3*nn,by=3)
    out$x[,'errAr39Ar40'] <- err[l]
    out$x[,'errAr36Ar40'] <- err[l+1]
    out$x[,'errAr39Ar36'] <- err[l+2]
    JJ <- subset(J,label="J:")
    if (length(JJ$intercepts)==1)
        out$J <- c(JJ$intercepts,sqrt(JJ$covmat))
    else if (abs(max(JJ$intercepts)-min(JJ$intercepts))<1e-10)
        out$J <- c(JJ$intercepts[1],sqrt(JJ$covmat[1,1]))
    else
        stop('All samples must have the same J-factor')
    out
}
