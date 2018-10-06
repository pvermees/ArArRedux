loadArgusData <- function(fname){
    header <- utils::read.csv(file=fname,header=FALSE,skip=1,nrows=1)
    dat <- utils::read.csv(file=fname,header=FALSE,skip=3)
    PHtags <- parseArgusHeader(header) # extract isotopes and detectors from header
    PH <- PHtags$PH
    tags <- PHtags$tags
    nrows <- dim(dat)[1] # number of MS runs
    ncols <- dim(dat)[2] # number of columns
    ntags <- length(tags) # number of isotopes or detectors
    ncycles <- (ncols-3)/(2*ntags)
    cirr <- 1
    cpos <- 2
    clabel <- 3
    cdate <- 4 # column with the date
    ci <- matrix(NA,nrow=ncycles,ncol=ntags) # column indices
    for (i in 1:ntags){
        ci[,i] <- seq(from=3+2*i,to=ncols,by=2*ntags)
    }
    if (PH) {
        out <- newPHdata(dat,tags,cirr,cpos,clabel,cdate,ci)
    } else {
        out <- newtimeresolved(dat,tags,cirr,cpos,clabel,cdate,ci)
    }
    class(out) <- append(class(out),"ArgusVI")
    return(out)
}

parseArgusHeader <- function(header){
    isotopes <- NULL
    detectors <- NULL
    for (tag in header){
        dat <- unlist(strsplit(as.character(tag),':'))
        if (length(dat)==3){
            detector <- dat[2]
            isotope <- dat[3]
            if (!(detector %in% detectors))
                detectors <- c(detectors,detector)
            if (!(isotope %in% isotopes))
                isotopes <- c(isotopes,isotope)
        }
    }
    if ((length(detectors)>1) && (length(isotopes)==1)){
        PH <- TRUE
        tags <- detectors
    }
    if ((length(detectors)>1) && (length(isotopes)>1)){
        PH <- FALSE
        tags <- paste0('Ar',isotopes)
    }
    if ((length(detectors)==1) | (length(isotopes)==1)){
        PH <- TRUE
        tags <- paste0('Ar',isotopes)
    }
    list(PH=PH,tags=tags)
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
    out <- list(PH=TRUE)
    class(out) <- "PHdata"
    out$masses <- masses
    for (i in 1:nmasses(out)){
        out$signals[[masses[i]]] <-
            newtimeresolved(thetable,masses[i],
                cirr,cpos,clabel,cdate,as.matrix(ci[,i]))
    }
    return(out)
}

timeresolvedplot.ArgusVI <- function(x,mass,label=NULL,run=1,...){
    if (!is.null(label))
        run <- which(x$labels==label)-1
    if (length(run)!=1){
        print('invalid input into plot function')
        return(NA)
    }
    k <- which(x$masses==mass)
    i <- run*nmasses(x)+k
    graphics::plot(x$thetime[,i],x$d[,i],type='p',
                   xlab='time',ylab=mass)
}
