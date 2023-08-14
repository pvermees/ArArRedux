setwd("/home/pvermees/Documents/Programming/R/ArArRedux/R")

rm(list=ls())

source("age.R")
source("Ar40Ar39.R")
source("ArArRedux.R")
source("Argus.R")
source("average.R")
source("blank.R")
source("calibration.R")
source("Cl.R")
source("decay.R")
source("documentation.R")
source("fractionation.R")
source("interference.R")
source("io.R")
source("isoratios.R")
source("J.R")
source("plot.R")
source("timezero.R")
source("toolbox.R")
source("WiscAr.R")
source("/home/pvermees/Documents/Programming/R/IsoplotR/R/york.R")
graphics.off()

setwd("/home/pvermees/Documents/Programming/R/WiscAr")

options(warn=2)

plotsignals <- function(dat,hop,tmin=0){
    nr <- nruns(dat[[hop]])
    for (i in 1:nr){
        plot(dat,run=i,hop=hop,tmin=tmin)
        readline(prompt="Press [enter] to continue")
    }
}

plotellipses <- function(x){
    tt <- x$ages*1000
    covmat <- x$covmat*1e6
    ns <- length(tt)
    par(mfrow=rep(ns,2),mar=rep(1.5,4),oma=c(3,3,0,0))
    minage <- min(tt - 2.5*sqrt(diag(covmat)))
    maxage <- max(tt + 2.5*sqrt(diag(covmat)))
    for (j in 1:ns){
        for (i in 1:ns){
            ell <- IsoplotR::ellipse(tt[i],tt[j],
                                     covmat=covmat[c(i,j),c(i,j)])
            if (i==1) ylab <- j
            else ylab <- ''
            if (j==ns) xlab <- i
            else xlab <- ''
            plot(ell,type='l',xlim=c(minage,maxage),
                 ylim=c(minage,maxage),xlab=xlab,ylab=ylab,xpd=NA)
            points(tt[i],tt[j],pch=20)
        }
    }
}

WiscArPreProcess <- function(dat,irr,irrlab,pos,Jpos,tmin=0,
                             perfectcal=FALSE,dodecaycorr=TRUE,
                             propDecayConstErr=TRUE){
    dat[['101']]$irr <- irrlab
    dat[['102']]$irr <- irrlab
    b <- blankcorr(dat,blanklabel="Blank")
    x <- fitlogratios(b,"Ar39",tmin=tmin)
    j <- newredux(x,Jpos=Jpos)
    if (!dodecaycorr){
        j$param$l7 <- 0
        j$param$l9 <- 0
        j$param$sl7 <- 0
        j$param$sl9 <- 0
    }
    if (!propDecayConstErr){
        j$param$sl7 <- 0
        j$param$sl9 <- 0
        j$param$sl0 <- 0
    }
    j$param$pcl <- 263
    j$param$spcl <- 2.0
    CC <- calibration_test(j,irr=irr,clabel='Air',
                           perfectcal=perfectcal,
                           dodecaycorr=dodecaycorr)
    CC$pos <- pos
    return(CC)
}

knockoutblanks <- function(dat,prefix="Blank"){
    out <- dat
    for (hop in names(dat)){
        i <- findrunindices(dat[[hop]],prefixes=prefix)
        j <- getindices(nmasses=length(dat[[hop]]$masses),iruns=i)
        out[[hop]]$d[,j] <- 0
    }
    out
}

calibration_test <- function(x,irradiations,clabel,perfectcal=FALSE,
                             dodecaycorr=TRUE){
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
        J[(4*i-3),i09s101[i]] <-  1
        if (dodecaycorr){
            J[(4*i-3),iD9[i]] <- -1
        }
        if (!perfectcal){
            J[(4*i-3),i09g101[j]] <- -1
            J[(4*i-3),i09w[j]]    <-  1
        }
        # s89corr
        J[(4*i-2),i89s102[i]] <-  1
        if (dodecaycorr){
            J[(4*i-2),iD9[i]] <- -1
        }
        if (!perfectcal){
            J[(4*i-2),i89g102[j]] <- -1
            J[(4*i-2),i89w[j]]    <-  1
        }
        # s79corr
        J[(4*i-1),i79s101[i]] <-  1
        if (dodecaycorr){
            J[(4*i-1),iD7[i]] <-  1
            J[(4*i-1),iD9[i]] <- -1
        }
        if (!perfectcal){
            J[(4*i-1),i89g102[j]] <- -1
            J[(4*i-1),i69g102[j]] <-  1
            J[(4*i-1),i89w[j]]    <-  1
            J[(4*i-1),i69w[j]]    <- -1
        }
        # s69corr
        J[4*i,i69s101[i]] <-  1/2
        J[4*i,i69s102[i]] <-  1/2
        if (dodecaycorr){
            J[4*i,iD9[i]] <- -1
        }
        if (!perfectcal){
            J[4*i,i69g101[j]] <- -1/2
            J[4*i,i69g102[j]] <- -1/2
            J[4*i,i69w[j]]    <-  1
        }
    }
    isamp <- findrunindices(x,prefixes=clabel,invert=TRUE)[1:ns]
    out <- subset(X,i=isamp)
    out$covmat <- J %*% X$covmat %*% t(J)
    out$intercepts <- J %*% X$intercepts
    out$hop <- NULL
    out$num <- rep(c("Ar40","Ar38","Ar37","Ar36"),ns)
    out$den <- rep("Ar39",4*ns)
    out$nlr <- rep(4,ns)
    out
}

get4039_test <- function(X,irr,doCl=TRUE,doAir=TRUE){
    Y <- getabcdef(X)
    ns <- nruns(Y)
    ni <- length(Y$intercepts)
    out <- Y
    out$intercepts <- rep(0,ns)
    out$covmat <- matrix(0,ns,ns)
    out$num <- rep(NA,ns)
    out$den <- rep(NA,ns)
    out$nlr <- rep(1,ns)
    J <- matrix(0,nrow=ns,ncol=ni)
    hasKglass <- "K-glass" %in% X$labels
    hasCasalt <- "Ca-salt" %in% X$labels
    for (i in 1:ns){
        j <- (i-1)*6
        label <- Y$labels[i]
        aa <- Y$intercepts[getindices(Y,label,num='a')]
        bb <- Y$intercepts[getindices(Y,label,num='b')]
        cc <- Y$intercepts[getindices(Y,label,num='c')]
        dd <- Y$intercepts[getindices(Y,label,num='d')]
        ee <- Y$intercepts[getindices(Y,label,num='e')]
        if (hasKglass){
            ff <- Y$intercepts[getindices(Y,label,num='f')]
        } else {
            ff <- 0
        }
        if (!hasCasalt){
            bb <- 0
            ee <- 0
        }
        if (!doAir){
            aa <- 0
            bb <- 0
            cc <- 0
        }
        if (!doCl){
            cc <- 0
        }
        J[i,j+1] <- -1/(dd-ee)
        J[i,j+2] <-  1/(dd-ee)
        J[i,j+3] <-  1/(dd-ee)
        J[i,j+4] <- -(1-aa+bb+cc)/(dd-ee)^2
        J[i,j+5] <-  (1-aa+bb+cc)/(dd-ee)^2
        if (hasKglass){
            J[i,j+6] <- -1
        }
        if (!hasCasalt){
            J[i,j+2] <- 0  # dR/db
            J[i,j+5] <- 0  # dR/de
        }
        if (!doAir){
            J[i,j+1] <- 0 # dR/da
            J[i,j+2] <- 0 # dR/db
            J[i,j+3] <- 0 # dR/dc
        }
        if (!doCl){
            J[i,j+3] <- 0 # dR/dc
        }
        out$intercepts[i] <- (1-aa+bb+cc)/(dd-ee)-ff
    }
    out$covmat <- J %*% Y$covmat %*% t(J)
    return(out)
}

crunch <- function(doblank=TRUE,docalibration=TRUE,dodecaycorr=TRUE,
                   doCaInterference=TRUE,doClInterference=TRUE,
                   doKinterference=TRUE,doAirCorr=TRUE,
                   propDecayConstErr=TRUE,propTSerr=TRUE){
    tmin <- 60
    setwd("/home/pvermees/Dropbox/Programming/R")
    irr <- loadirradiations('WiscAr/irradiations.csv')
    dat1 <- loaddata("WiscAr/NAH",MS='WiscAr')
    if (!doblank) dat1 <- knockoutblanks(dat1,prefix="Blank")
    irrlab <- c(NA,'CAL','UW137','UW137',NA,'CAL',
                'UW137',NA,'UW137',NA,'UW137')
#    irrlab <- c(NA,'CAL','UW137',NA,'CAL',
#                'UW137',NA,'CAL')
    Jpos <- c(1,3)
    CCsamp <- WiscArPreProcess(dat1,irr,irrlab,pos=rep(2,4),
                               Jpos=Jpos,tmin=tmin,
                               perfectcal=!docalibration,
                               dodecaycorr=dodecaycorr,
                               propDecayConstErr=propDecayConstErr)
    dat2 <- loaddata("WiscAr/NAG",MS='WiscAr')
    if (!doblank) dat2 <- knockoutblanks(dat2,prefix="Blank")
    irrlab <- c(NA,'UW137',NA,'CAL','UW137',NA,'UW137',NA,
                'CAL','UW137',NA,'UW137',NA,'CAL','UW137',NA,'UW137',
                NA,'CAL','UW137',NA,'UW137',NA,'CAL','UW137',NA,
                'UW137',NA,'CAL','UW137')
    CCstand <- WiscArPreProcess(dat2,irr,irrlab,
                                pos=c(rep(1,6),rep(3,6)),
                                Jpos=Jpos,tmin=tmin,
                                perfectcal=!docalibration,
                                dodecaycorr=dodecaycorr,
                                propDecayConstErr=propDecayConstErr)
    X <- concat(list(CCsamp,CCstand))
    Y <- renormalise.WiscAr(X)
    if (FALSE){ # get Ar40Ar36 ratios of irradiation standards
        prefix <- 'S-228'
        iJ <- getindices(Y,prefix=prefix,num='Ar36')
        Ar40Ar36 <- exp(-Y$intercepts[iJ])
        riJ <- findrunindices(Y,prefixes=prefix)
        cbind(Y$labels[riJ],Ar40Ar36)
    }
    Y$param$ts <- 1.1864
    if (propTSerr) Y$param$sts <- 0.0006
    else Y$param$sts <- 0
    if (doKinterference){
        k <- interference(intercepts=log(c(0.0121,0.00054)),
                          covmat=matrix(c(0.00000361,0,0,0.06718464),nrow=2),
                          num=c("Ar38","Ar40"),
                          den=c("Ar39","Ar39"),
                          irr="UW137",label="K-glass")
        Y <- concat(list(Y,k))
    }
    if (doCaInterference){
        ca <- interference(intercepts=log(c(0.000265,0.000695)),
                           num=c("Ar39","Ar36"),
                           den=c("Ar37","Ar37"),
                           irr="UW137",label="Ca-salt",
                           covmat=matrix(c(6.892132e-05,0,0,0.0001676932),
                                         nrow=2))
        Y <- concat(list(Y,ca))
    }
    Z <- clcorrection(Y,irr)
    IR <- isoratios(Y,irr=irr,inverse=TRUE)
    yd <- data2york(IR)
    R <- get4039_test(Z,irr,
                      doCl=doClInterference,
                      doAir=doAirCorr)
    if (TRUE){ # only use low Ar36 J-standards
        R <- subset(R,i=c(1,2,3,4,7,11)) # 7, 11
    }
    J <- getJfactors(R)
    ages <- getages(J)
    rbind(ages$ages,sqrt(diag(ages$covmat)))
}

collate <- function(i,tests){
    nr <- length(tests)
    ns <- ncol(tests[[1]])
    out <- matrix(NA,nr,ns)
    rownames(out) <- names(tests)
    for (j in 1:nr){
        out[j,] <- tests[[j]][i,]
    }
    out
}

accuracyVSprecision <- function(tt,st,sn=1){
    accuracy <- tt[-1,]
    precision <- st[-1,]
    ns <- nrow(tt)
    for (i in 2:ns){
        accuracy[i-1,] <- 100*(tt[i,]-tt[1,])/tt[1,]
        precision[i-1,] <- 100*(st[i,]-st[1,])/tt[1,]
    }
    colour <- rep('blue',ns)
    neg <- which(accuracy<0)
    colour[neg] <- 'red'
    if (FALSE){
        plot(accuracy[,sn],precision[,sn],
             xlab='% change in age (accuracy)',
             ylab='% change in precision',pch='n')
        text(accuracy[,sn],precision[,sn],labels=rownames(accuracy))
    } else {
        plot(accuracy[,sn],precision[,sn],
             xlab='% change in age (accuracy)',
             ylab='% change in precision',
             type='p',pch=19,bg='black',cex=1)   
        text(accuracy[,sn],precision[,sn],
             labels=rownames(accuracy),pos=1)
    }
}
ageVSprecision <- function(tt,st,sn=1){
    nr <- nrow(tt)
    plot(1000*tt[,sn],1000*st[,sn],
         type='n',xlab='age (ka)',ylab='error (1s)')
    points(1000*tt[1,sn],1000*st[1,sn],
           pch=22,bg='white',cex=1)
    points(1000*tt[2:nr,sn],1000*st[2:nr,sn],
         pch=19,bg='black',cex=1)
    text(1000*tt[,sn],1000*st[,sn],labels=rownames(tt),pos=1)
}

sensitivityTest <- function(sn=1,doall=TRUE){
    tests <- list()
    tests$accurate <- crunch()
    if (doall){
        tests$calibration <- crunch(docalibration=FALSE)
        tests$blank <- crunch(doblank=FALSE)
        tests$air <- crunch(doAirCorr=FALSE)
    }
    tests$decay <- crunch(dodecaycorr=FALSE)
    tests$Cl <- crunch(doClInterference=FALSE)
    tests$K <- crunch(doKinterference=FALSE)
    tests$Ca <- crunch(doCaInterference=FALSE)
    tests$lambda <- crunch(propDecayConstErr=FALSE)
    tests$sts <- crunch(propTSerr=FALSE)
    tt <- collate(1,tests)
    st <- collate(2,tests)
    list(tt=tt,st=st)
}

option <- 4

if (option==1){
    tmin <- 60
    dat2 <- loaddata("NAH",MS='WiscAr')
    plotsignals(dat2,hop='101',tmin=tmin)
} else if (option==2){
    plotellipses(x)
} else if (option==3){
    sn <- 2
    X11(width=14,height=7)
    par(mfrow=c(1,2))
    tst <- sensitivityTest(sn=sn,doall=TRUE)
    #ageVSprecision(tst$tt,tst$st,sn=sn)
    pie(abs(tst$st[2:nrow(tst$tt),sn]-tst$st[1,sn]))
    tst <- sensitivityTest(sn=sn,doall=FALSE)
    #ageVSprecision(tst$tt,tst$st,sn=sn)
    pie(abs(tst$st[2:nrow(tst$tt),sn]-tst$st[1,sn]))
} else if (option==4){
    out <- crunch()
}
