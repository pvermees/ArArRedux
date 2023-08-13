#' Isochron ratios
#'
#' Calculate the Calculate the inverse isochron ratios of interference
#' corrected logratio intercept data
#' 
#' @param X an object of class \code{redux} containing some
#'     interference corrected logratio intercept data
#' @param irr the irradiation schedule
#' @param fract list with air shot data (one item per denominator
#'     isotope)
#' @param ca an object of class \code{\link{logratios}} with
#'     Ca-interference data (not necessary if interferences are
#'     included in X)
#' @param k an object of class \code{\link{logratios}} with
#'     K-interference data (not necessary if interferences are
#'     included in X)
#' @param inverse logical. If \code{TRUE}, returns inverse isochron
#'     ratios (\eqn{^{36}}Ar/\eqn{^{40}}Ar
#'     vs. \eqn{^{39}}Ar/\eqn{^{40}}Ar). Otherwise returns normal
#'     isochron ratios (\eqn{^{40}}Ar/\eqn{^{36}}Ar
#'     vs. \eqn{^{39}}Ar/\eqn{^{36}}Ar).
#' @return an object of class \code{link{isoratios}} containing the
#'     isochron ratios and their covariance matrix.
#' @examples
#' data(Melbourne)
#' IR <- isoratios(Melbourne$X,irr=Melbourne$irr,fract=Melbourne$fract)
#' @export
isoratios <- function(X,irr,fract=NULL,ca=NULL,k=NULL,inverse=TRUE){
    Cl <- corrections(X=X,irr=irr,fract=fract,ca=ca,k=k)
    Y <- getABCDEFGHI(Cl)
    ns <- nruns(Y)
    ni <- length(Y$intercepts)
    out <- Y
    out$inverse <- inverse
    out$intercepts <- rep(0,2*ns)
    out$covmat <- matrix(0,2*ns,2*ns)
    if (inverse){
        out$num <- rep(c('Ar39','Ar36'),ns)
        out$den <- rep('Ar40',2*ns)
    } else {
        out$num <- rep(c('Ar39','Ar40'),ns)
        out$den <- rep('Ar36',2*ns)
    }
    out$nlr <- rep(2,ns)
    J <- matrix(0,nrow=2*ns,ncol=ni)
    hasKglass <- "K-glass" %in% Cl$labels
    hasCasalt <- "Ca-salt" %in% Cl$labels
    for (i in 1:ns){
        j <- (i-1)*9
        k <- (i-1)*2
        label <- Y$labels[i]
        aa <- Y$intercepts[getindices(Y,label,num='A')]
        bb <- Y$intercepts[getindices(Y,label,num='B')]
        cc <- Y$intercepts[getindices(Y,label,num='C')]
        dd <- Y$intercepts[getindices(Y,label,num='D')]
        ee <- Y$intercepts[getindices(Y,label,num='E')]
        ff <- Y$intercepts[getindices(Y,label,num='F')]
        gg <- Y$intercepts[getindices(Y,label,num='G')]
        hh <- Y$intercepts[getindices(Y,label,num='H')]
        ii <- Y$intercepts[getindices(Y,label,num='I')]
        if (!hasKglass){
            dd <- ee <- ff <- gg <- 0
        }
        if (!hasCasalt | expired(irr[[Y$irr[i]]],Y$thedate[i],Y$param$l7)){
            bb <- ee <- gg <- ii <- 0
        }
        if (inverse){
            out$intercepts[k+1] <- (hh-ii)/(1-ff+gg)          # 90
            J[k+1,j+6] <- (hh-ii)/(1-ff+gg)^2                 # d90/df
            J[k+1,j+7] <- -(hh-ii)/(1-ff+gg)^2                # d90/dg
            J[k+1,j+8] <- 1/(1-ff+gg)                         # d90/dh
            J[k+1,j+9] <- -1/(1-ff+gg)                        # d90/di
            out$intercepts[k+2] <- (aa-bb-cc+dd-ee)/(1-ff+gg) # 60
            J[k+2,j+1] <- 1/(1-ff+gg)                         # d60/da
            J[k+2,j+2] <- -1/(1-ff+gg)                        # d60/db
            J[k+2,j+3] <- -1/(1-ff+gg)                        # d60/dc
            J[k+2,j+4] <- 1/(1-ff+gg)                         # d60/dd
            J[k+2,j+5] <- -1/(1-ff+gg)                        # d60/de
            J[k+2,j+6] <- (aa-bb-cc+dd-ee)/(1-ff+gg)^2        # d60/df
            J[k+2,j+7] <- -(aa-bb-cc+dd-ee)/(1-ff+gg)^2       # d60/dg
        } else {
            out$intercepts[k+1] <- (hh-ii)/(aa-bb-cc+dd-ee)   # 96
            J[k+1,j+1] <- -(hh-ii)/(aa-bb-cc+dd-ee)^2         # d96/da
            J[k+1,j+2] <- (hh-ii)/(aa-bb-cc+dd-ee)^2          # d96/db
            J[k+1,j+3] <- (hh-ii)/(aa-bb-cc+dd-ee)^2          # d96/dc
            J[k+1,j+4] <- -(hh-ii)/(aa-bb-cc+dd-ee)^2         # d96/dd
            J[k+1,j+5] <- (hh-ii)/(aa-bb-cc+dd-ee)^2          # d96/de
            J[k+1,j+8] <- 1/(aa-bb-cc+dd-ee)                  # d96/dh
            J[k+1,j+9] <- -1/(aa-bb-cc+dd-ee)                 # d96/di
            out$intercepts[k+2] <- (1-ff+gg)/(aa-bb-cc+dd-ee) # 06
            J[k+2,j+1] <- -(1-ff+gg)/(aa-bb-cc+dd-ee)^2       # d06/da
            J[k+2,j+2] <- (1-ff+gg)/(aa-bb-cc+dd-ee)^2        # d06/db
            J[k+2,j+3] <- (1-ff+gg)/(aa-bb-cc+dd-ee)^2        # d06/dc
            J[k+2,j+4] <- -(1-ff+gg)/(aa-bb-cc+dd-ee)^2       # d06/dd
            J[k+2,j+5] <- (1-ff+gg)/(aa-bb-cc+dd-ee)^2        # d06/de
            J[k+2,j+6] <- -1/(aa-bb-cc+dd-ee)                 # d06/df
            J[k+2,j+7] <- 1/(aa-bb-cc+dd-ee)                  # d06/dg
        }
    }
    out$covmat <- J %*% Y$covmat %*% t(J)
    class(out) <- append(class(out),"isoratios")
    return(out)
}

getJABCDEFGHI <- function(Z,Slabels,nl){
    J <- matrix(0,nrow=nl,ncol=length(Z$intercepts))
    i67ca <- getindices(Z,"Ca-salt","Ar36","Ar37")
    i97ca <- getindices(Z,"Ca-salt","Ar39","Ar37")
    i09k <- getindices(Z,"K-glass","Ar40","Ar39")
    i89k <- getindices(Z,"K-glass","Ar38","Ar39")
    for (i in 1:length(Slabels)){
        j <- (i-1)*9
        label <- Slabels[i]
        i60 <- getindices(Z,label,"Ar36","Ar40")
        i70 <- getindices(Z,label,"Ar37","Ar40")
        i80 <- getindices(Z,label,"Ar38","Ar40")
        i90 <- getindices(Z,label,"Ar39","Ar40")
        i68cl <- getindices(Z,paste0("Cl:",label),"Ar36","Ar38")
        J[j+1,i60] <- 1
        J[j+2,c(i67ca,i70)] <- 1
        J[j+3,c(i68cl,i80)] <- 1
        J[j+4,c(i90,i89k,i68cl)] <- 1
        J[j+5,c(i70,i97ca,i89k,i68cl)] <- 1
        J[j+6,c(i90,i09k)] <- 1
        J[j+7,c(i70,i97ca,i09k)] <- 1
        J[j+8,i90] <- 1
        J[j+9,c(i70,i97ca)] <- 1
    }
    return(J)
}

getABCDEFGHI <- function(Cl){
    Z <- concat(list(Cl,air(Cl))) # matrix with everything
    i <- findrunindices(Z,c("Ca-salt","K-glass","Cl:"),invert=TRUE)
    j <- findmatches(Z$pos,NA,invert=TRUE)
    ii <- 1:nruns(Z)
    si <- which((ii %in% i) & (ii %in% j))
    ns <- length(si)
    out <- subset(Z,si)
    out$num <- c(rep(c("A","B","C","D","E","F","G","H","I"),ns))
    out$den <- rep(NA,9*ns)
    out$nlr <- rep(9,ns)
    Jv <- getJABCDEFGHI(Z,Z$labels[si],length(out$num))
    out$intercepts <- exp(Jv %*% Z$intercepts)
    Jw <- apply(Jv,2,"*",out$intercepts)
    out$covmat <- Jw %*% Z$covmat %*% t(Jw)
    return(out)
}
