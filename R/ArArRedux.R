#' Process logratio data and calculate 40Ar/39Ar ages
#'
#' Performs detector calibration, mass fractionation correction, decay
#' corrections, interference corrections, interpolates J-factors and
#' calculates ages.
#'
#' @param X an object of class \code{\link{redux}}
#' @param irr the irradiation schedule
#' @param fract list with air shot data (one item per denominator isotope)
#' @param ca an object of class \code{\link{logratios}} with
#' Ca-interference data (not necessary if interferences are included in X)
#' @param k an object of class \code{\link{logratios}} with
#' K-interference data (not necessary if interferences are included in X)
#' @examples
#' data(Melbourne)
#' ages <- process(Melbourne$X,Melbourne$irr,Melbourne$fract)
#' summary(ages)
#' @export
process <- function(X,irr,fract=NULL,ca=NULL,k=NULL){
    Cl <- corrections(X,irr,fract=fract,ca=ca,k=k)
    # calculate the 40Ar*/39ArK-ratios 
    R <- get4039(Cl,irr)
    # calculate J factors
    J <- getJfactors(R)
    # calculate ages
    ages <- getages(J)
    return(ages)
}

# apply calibration, fractionation, decay and interference corrections
corrections <- function(X,irr,fract=NULL,ca=NULL,k=NULL){
    # apply the detector calibration (this won't affect the Ar40/Ar36 ratio)
    C <- calibration(X,"DCAL")
    # apply the mass fractionation correction
    if (is.null(fract)){
        A <- C
    } else {
        A <- massfractionation(C,fract)
    }
    # decay corrections
    D9 <- decaycorrection(A,irr,"Ar39")
    D7 <- decaycorrection(D9,irr,"Ar37")
    if (is.null(k)){
    # interference corrections
        K <- average(D7,grep("K:",A$labels),newlabel="K-glass")
    } else {
        K <- concat(list(D7,k))
    }
    if (is.null(ca)){
        Ca <- average(K,grep("Ca:",K$labels),newlabel="Ca-salt")
    } else {
        Ca <- concat(list(K,ca))
    }
    clcorrection(Ca,irr)
}
