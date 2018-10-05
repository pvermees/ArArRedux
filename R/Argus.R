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
