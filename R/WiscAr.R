loadWiscArData <- function(dname){ 
    out <- list()
    filenames <- Sys.glob(paste0(dname,"/*.RUN"))
    nf <- length(filenames)
    out$thedate <- rep(NA,nf)
    out$labels <- rep(NA,nf)
    out$irr <- rep(NA,nf)
    out$pos <- rep(NA,nf)
    out$hops=list('1'=c('IC0','IC1','Far','IC2','IC3'),
                  '101'=c('Ar40','Ar39','Ar38','Ar37','Ar36'),
                  '102'=c('Ar39','Ar38','Ar37','Ar36','35'))
    out[['101']] <- list(masses=out$hops[['101']],d=NULL,thetime=NULL)
    out[['102']] <- list(masses=out$hops[['102']],d=NULL,thetime=NULL)
    for (i in 1:nf){
        con = file(filenames[i], "r")
        while (TRUE){
            line = readLines(con,n=1)
            if (grepl("Start Of Run",line)){
                txt <- unlist(strsplit(line,'#'))[2]
                out$thedate[i] <- readthedate(txt)
            } else if (grepl(".mdf",line,ignore.case=TRUE)){
                line = readLines(con,n=1) # the next line contains the sample name
                sname <- unlist(strsplit(line,'"'))[2]
                naliquots <- length(names(out) %in% sname)
                if (naliquots > 0) sname <- paste(sname,naliquots,sep=' ')
                out$labels[i] <- sname
            } else if (identical(line,'""')){ # start of signal
                out <- addWiscArSignal(out,con)
                break
            }
        }
        close(con)
    }
    class(out) <- c('timeresolved','WiscAr')
    out
}

addWiscArSignal <- function(dat,con){
    out <- dat
    thetime <- list()
    d <- list()
    while (TRUE){
        line = readLines(con,n=1)
        if (length(line)==0) break
        txt <- unlist(strsplit(line,','))
        hop <- txt[[7]]
        if (hop %in% names(dat$hops)){
            thetime[[hop]] <- append(thetime[[hop]],as.numeric(txt[6]))
            d[[hop]] <- rbind(d[[hop]],as.numeric(txt[1:5]))
        }
    }
    for (hop in names(dat$hops)){
        out[[hop]]$thetime <- cbind(out[[hop]]$thetime,thetime[[hop]])
        out[[hop]]$d <- cbind(out[[hop]]$d,d[[hop]])
    }
    out
}

plottimeresolved.WiscAr <- function(x,mass='Ar40',label=NULL,run=1,...){
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
