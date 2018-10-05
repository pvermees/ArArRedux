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
loaddata <- function(fname,MS="ARGUS-VI"){
    if (identical(MS,'ARGUS-VI'))
        out <- loadArgusData(fname)
    else if (identical(MS,'WiscAr'))
        out <- loadWiscAr(fname)
    out
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
#' @param x a \code{.csv} file with samples and fluence monitor data
#'     (if \code{MS = 'ARGUS-VI'}) OR a directory with \code{.RUN}
#'     files (if \code{MS = 'WiscAr'})
#' @param blabel a prefix string denoting the blanks
#' @param Jpos a vector of integers denoting the positions of the
#'     fluence monitors in the irradiation stack
#' @param kdat a \code{.csv} file with the K-interference monitor data
#'     (if \code{MS = 'ARGUS-VI'}) OR a directory with \code{.RUN}
#'     files of K-interference monitor data (if \code{MS = 'WiscAr'})
#'     (optional)
#' @param cadat a \code{.csv} file with the Ca-interference monitor
#'     data (if \code{MS = 'ARGUS-VI'}) OR a directory with
#'     \code{.RUN} files of Ca-interference monitor data (if \code{MS
#'     = 'WiscAr'}) (optional)
#' @param ddat a \code{.csv} file with the detector calibration data
#'     (if \code{MS = 'ARGUS-VI'}) OR the prefix of the standard gas
#'     measurements (if \code{MS = 'WiscAr'}) (optional)
#' @param MS a string denoting the type of mass spectrometer. At the
#'     moment there are two options: \code{'ARGUS-VI'}, for data that
#'     are acquired on the eponymous instrument following the
#'     measurement protocols at Melbourne University, and
#'     \code{'WiscAr'} for the 5-detector Noblesse instrument at the
#'     University of Wisconsin. Please contact the author to add other
#'     file formats to Ar-Ar_Redux.
#' @return an object of class \code{\link{redux}}.
#' @examples
#' samplefile <-  system.file("Samples.csv",package="ArArRedux")
#' kfile <- system.file("K-glass.csv",package="ArArRedux")
#' cafile <- system.file("Ca-salt.csv",package="ArArRedux")
#' dfile <- system.file("Calibration.csv",package="ArArRedux")
#' X <- read(samplefile,blabel="EXB#",Jpos=c(3,15),
#'           kdat=kfile,cadat=cafile,ddat=dfile)
#' plotcorr(X)
#' @export
read <- function(x,blabel,Jpos,kdat=NULL,cadat=NULL,ddat=NULL,MS="ARGUS-VI"){
    # load the .csv files
    m <- loaddata(x,MS) # samples and J-standards
    x <- fitlogratios(blankcorr(m,blabel),"Ar40")
    if (!is.null(kdat)){ # K-interference data
        # subset of the relevant isotopes
        mk <- loaddata(kdat,MS)
        lk <- fitlogratios(blankcorr(mk,blabel,"K:"),"Ar39")
        k <- getmasses(lk,"Ar40","Ar39")
        x <- concat(list(x,k))
    }
    if (!is.null(cadat)){ # Ca interference data
        mca <- loaddata(cadat,MS)
        lca <- fitlogratios(blankcorr(mca,blabel,"Ca:"),"Ar37")
        ca <- getmasses(lca,c("Ar36","Ar39"),c("Ar37","Ar37"))
        x <- concat(list(x,ca))
    }
    if (!is.null(ddat)){
        md <- loaddata(ddat,MS)
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
#'     current parameter values if no optional arguments are supplied.
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
