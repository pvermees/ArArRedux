test <- function(builddata=FALSE,option=TRUE){
    
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
    
    if (option){ # full propagation
        X <- read(samplefile,masses,blanklabel,Jpos,
                  kfile,cafile,dfile,dlabels)
        irr <- loadirradiations(irrfile)
        fract <- list(fractionation(fd37file,"L2",PH=TRUE),
                      fractionation(fd39file,"AX",PH=TRUE),
                      fractionation(fd40file,"H1",PH=FALSE))
        ages <- process(X,irr,fract)
        return(ages)
    } else {
        mMC <- loaddata(samplefile,masses)
        graphics::plot(mMC,"MD2-1a","Ar37")
    }

    if (builddata){
        Melbourne <- list(X=X,irr=irr,fract=fract)
        save(Melbourne,file="../data/Melbourne.rda")
    }

}
