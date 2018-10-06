loadWiscArData <- function(dname){ 
    out <- list()
    filenames <- Sys.glob(paste0(dname,"/*.RUN"))
    nf <- length(filenames)
    labels=list('1'=c('IC0','IC1','Far','IC2','IC3'),
                '101'=c('Ar40','Ar39','Ar38','Ar37','Ar36'),
                '102'=c('Ar39','Ar38','Ar37','Ar36','35'))
    for (i in 1:nf){
        con = file(filenames[i], "r")
        while (TRUE){
            line = readLines(con,n=1)
            if (grepl("Start Of Run",line)){
                txt <- unlist(strsplit(line,'#'))[2]
                thedate <- readthedate(txt)
            } else if (grepl(".mdf",line)){
                line = readLines(con,n=1) # the next line contains the sample name
                sname <- unlist(strsplit(line,'"'))[2]
                naliquots <- length(names(out) %in% sname)
                if (naliquots > 0) sname <- paste(sname,naliquots,sep=' ')
            } else if (identical(line,'""')){ # start of signal
                d <- readWiscArSignal(con,masses=labels)
                out[[sname]] <- list(thedate=thedate,d=d)
                break
            }
        }
        close(con)
    }
    class(out) <- c('WiscAr','timeresolved')
    out
}

readWiscArSignal <- function(con,masses){
    out <- list()
    while (TRUE){
        line = readLines(con,n=1)
        if (length(line)==0) break
        txt <- unlist(strsplit(line,','))
        tag <- txt[[7]]
        if (tag %in% names(masses)){
            out[[tag]] <- rbind(out[[tag]],
                                as.numeric(txt[c(6,1,2,3,4,5)]))
        }
    }
    for (tag in names(masses)){
        colnames(out[[tag]]) <- c('time',masses[[tag]])
    }
    out
}

timeresolvedplot.WiscAr <- function(x,mass='Ar40',label=NULL,run=1,...){
    if (!is.null(label))
        run <- which(label %in% names(x))
    samp <- x[[run]]
    dat <- samp$d
    X <- NULL
    Y <- NULL
    hops <- NULL
    for (hop in names(dat)){
        signal <- dat[[hop]]
        if (mass %in% colnames(signal)){
            hops <- c(hops,mass)
            X <- cbind(X,signal[,'time'])
            Y <- cbind(Y,signal[,mass])
        }
    }
    graphics::plot(range(X),range(Y),type='n',
                   xlab='time',ylab=mass)
    title(names(x)[[run]])
    cols <- c('red','blue','green')
    for (i in 1:ncol(X)){
        graphics::points(X[,i],Y[,i],col=cols[i],...)
        graphics::lines(X[,i],Y[,i],col=cols[i],...)
    }
}
