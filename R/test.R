test <- function(option='full'){
    
    samplefile <- "../inst/Samples.csv"
    kfile <- "../inst/K-glass.csv"
    cafile <- "../inst/Ca-salt.csv"
    fd37file <- "../inst/AirL2.csv"
    fd39file <- "../inst/AirAX.csv"
    fd40file <- "../inst/AirH1.csv"
    irrfile <- "../inst/irradiations.csv"
    dfile <- "../inst/Calibration.csv"
    masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
    blanklabel <- "EXB#"
    dlabels <- c("H1","AX","L1","L2")
    Jpos <- c(3,15)
    
    X <- read(samplefile,masses,blanklabel,Jpos,
              kfile,cafile,dfile,dlabels)
    irr <- loadirradiations(irrfile)
    fract <- list(fractionation(fd37file,"L2",PH=TRUE),
                  fractionation(fd39file,"AX",PH=TRUE),
                  fractionation(fd40file,"H1",PH=FALSE))
    if (identical(option,'full')){ # full propagation
        ages <- process(X,irr,fract)
        return(ages)
    } else if (identical(option,'simple')){
        mMC <- loaddata(samplefile,masses)
        graphics::plot(mMC,"MD2-1a","Ar37")
    } else if (identical(option,'builddata')){
        Melbourne <- list(X=X,irr=irr,fract=fract)
        save(Melbourne,file="../data/Melbourne.rda")
    } else if (identical(option, 'subset')){
        ages <- process(X,irr,fract)
        out <- subset(ages,labels=c("MD2-1","MD2-2","MD2-3","MD2-4","MD2-5"))
    } else if (identical(option, 'isoplotr')){
        X <- read(samplefile,masses,blanklabel,Jpos)
        Y <- subset(X,labels=c("MD2-2h","MD2-2i","MD2-2j","MD2-2k","MD2-2l",
                               "MD2-2m","MD2-2n","MD2-2o","MD2-2p","MD2-2q",
                               "MD2-2r","MD2-2s","MD2-2t","MD2-2u"),include.J=TRUE)
        out <- redux2isoplotr(Y,irr,format=1,
                              file='/home/pvermees/Dropbox/Programming/R/IsoplotR/inst/ArAr3.csv')
        #out <- redux2isoplotr(Y,irr,fract=fract,format=2)
#        ArAr <- out; save(ArAr,file="/home/pvermees/Dropbox/Programming/R/IsoplotR/data/ArAr.rda")
    }
    out
}
