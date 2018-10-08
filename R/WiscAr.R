loadWiscArData <- function(dname){ 
    out <- list()
    filenames <- Sys.glob(paste0(dname,"/*.RUN"))
    nf <- length(filenames)
    thedate <- rep(NA,nf)
    labels <- rep(NA,nf)
    irr <- rep(NA,nf)
    pos <- rep(NA,nf)
    hops=list('101'=c('Ar40','Ar39','Ar38','Ar37','Ar36'),
              '102'=c('Ar39','Ar38','Ar37','Ar36','35'))
            # '1'=c('IC0','IC1','Far','IC2','IC3') # not recorded
    for (i in 1:nf){
        con = file(filenames[i], "r")
        while (TRUE){
            line = readLines(con,n=1)
            if (grepl("Start Of Run",line)){
                txt <- unlist(strsplit(line,'#'))[2]
                thedate[i] <- readthedate(txt)
            } else if (grepl(".mdf",line,ignore.case=TRUE)){
                line = readLines(con,n=1) # the next line contains the sample name
                sname <- unlist(strsplit(line,'"'))[2]
                naliquots <- length(names(out) %in% sname)
                if (naliquots > 0) sname <- paste(sname,naliquots,sep=' ')
                labels[i] <- sname
            } else if (identical(line,'""')){ # start of signal
                out <- addWiscArSignal(out,con,hops,thedate=thedate,
                                       thelabel=thelabel,theirr=NA,thepos=NA)
                break
            }
        }
        close(con)
    }
    for (hop in names(hops)){
        out[[hop]]$thedate <- thedate
        out[[hop]]$labels <- labels
        out[[hop]]$irr <- irr
        out[[hop]]$pos <- pos
        class(out[[hop]]) <- 'timeresolved'
    }
    class(out) <- 'WiscAr'
    out
}

addWiscArSignal <- function(dat,con,hops,thedate=NA,thelabel=NA,theirr=NA,thepos=NA){
    out <- dat
    thetime <- list()
    d <- list()
    while (TRUE){
        line = readLines(con,n=1)
        if (length(line)==0) break
        txt <- unlist(strsplit(line,','))
        hop <- txt[[7]]
        if (hop %in% names(hops)){
            thetime[[hop]] <- append(thetime[[hop]],as.numeric(txt[6]))
            d[[hop]] <- rbind(d[[hop]],as.numeric(txt[1:5]))
        }
    }
    for (hop in names(hops)){
        out[[hop]]$masses <- hops[[hop]]
        out[[hop]]$thetime <- cbind(out[[hop]]$thetime,thetime[[hop]])
        out[[hop]]$d <- cbind(out[[hop]]$d,d[[hop]])
    }
    out
}
