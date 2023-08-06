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
    keep <- list('101'=c(TRUE,TRUE,FALSE,TRUE,TRUE),
                 '102'=c(TRUE,TRUE,FALSE,TRUE,FALSE))
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
                naliquots <- length(which(grepl(sname,labels)))
                if (naliquots > 0) sname <- paste(sname,naliquots,sep=' ')
                labels[i] <- sname
            } else if (identical(line,'""')){ # start of signal
                out <- addWiscArSignal(out,con,hops,keep,thedate=thedate)
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

addWiscArSignal <- function(dat,con,hops,keep,thedate=NA,
                            thelabel=NA,theirr=NA,thepos=NA){
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
            sweep <- txt[1:5]
            d[[hop]] <- rbind(d[[hop]],
                              as.numeric(sweep[keep[[hop]]]))
        }
    }
    for (hop in names(hops)){
        i <- keep[[hop]]
        out[[hop]]$masses <- hops[[hop]][i]
        out[[hop]]$thetime <- cbind(out[[hop]]$thetime,thetime[[hop]])
        out[[hop]]$d <- cbind(out[[hop]]$d,d[[hop]])
    }
    out
}

Wiscalgas <- function(){
    # standard gas composition from WiscAr spreadsheet
    meas <- log(c(4.2055,69.948,0.9714))
    names(meas) <- c('09','96','68')
    covmat <- matrix(0,3,3)
    diag(covmat) <- (c(0.0027/4.2055,
                       0.076/69.948,
                       0.0015/0.9714)/2)^2
    out <- list()
    out$irr <- rep(NA,3)
    out$pos <- rep(NA,3)
    out$labels <- 'CAL'
    out$thedate <- readthedate("3-12-2016 12:00")
    out$num <- c('Ar40','Ar38','Ar36')
    out$den <- 'Ar39'
    out$hop <- NA
    out$nlr <- 3
    out$intercepts <- rep(0,3)
    out$intercepts[1] <- meas['09']
    out$intercepts[2] <- - meas['96']
    out$intercepts[3] <- - meas['68'] - meas['96']
    J <- matrix(0,3,3)
    J[1,1] <- 1
    J[2,2] <- -1
    J[3,2] <- -1
    J[3,3] <- -1
    out$covmat <- J %*% covmat %*% t(J)
    out
}
renormalise.WiscAr <- function(x){
    i0 <- getindices(x,num='Ar40')
    i8 <- getindices(x,num='Ar38')
    i7 <- getindices(x,num='Ar37')
    i6 <- getindices(x,num='Ar36')
    ns <- length(x$labels)
    J <- matrix(0,4*ns,4*ns)
    for (i in 1:ns){
        J[4*i-3,i6[i]] <- 1
        J[4*i-2,i7[i]] <- 1
        J[4*i-1,i8[i]] <- 1
        J[4*i-(0:3),i0[i]] <- -1
    }
    out <- x
    out$intercepts <- J %*% x$intercept
    out$covmat <- J %*% x$covmat %*% t(J)
    out$num <- rep(c("Ar36","Ar37","Ar38","Ar39"),ns)
    out$den <- rep("Ar40",4*ns)
    out
}
