#' Apply a blank correction
#'
#' Applies a blank correction to some time-resolved mass spectrometer data
#' 
#' @param x an object of class \code{\link{timeresolved}} or
#' \code{\link{PHdata}}
#' @param ... other arguments
#' @param blanklabel as string denoting the prefix of the blanks
#' @param prefix a string to be prepended to the non-blanks
#' @return an object of class \code{\link{blankcorrected}}
#' @examples
#' samplefile <- system.file("Samples.csv",package="ArArRedux")
#' masses <- c("Ar37","Ar38","Ar39","Ar40","Ar36")
#' m <- loaddata(samplefile,masses) # samples and J-standards
#' blanklabel <- "EXB#"
#' l <- fitlogratios(blankcorr(m,blanklabel),"Ar40")
#' plotcorr(l)
#' @export
blankcorr <- function(x,...){ UseMethod("blankcorr",x) }
#' @rdname blankcorr
#' @export
blankcorr.default <- function(x,...){stop()}
#' @rdname blankcorr
#' @export
blankcorr.timeresolved <- function(x,blanklabel=NULL,prefix='',...){
    if (is.null(blanklabel)){
        out <- x
        out$blankindices <- 1:nruns(x)
    } else {
        # find indices of the blanks and non-blanks
        iblanks <- array(grep(blanklabel,x$labels))
        iothers <- array(grep(blanklabel,x$labels,invert=TRUE))
        blanks <- subset(x,iblanks)
        others <- subset(x,iothers)
        out <- others
        out$labels <- unlist(lapply(prefix,paste0,others$labels))
        inearestblanks <- nearest(others$thedate,blanks$thedate)
        out$d <- others$d - getruns(blanks,inearestblanks)
        out$blankindices <- as.vector(inearestblanks)
    }
    class(out) <- append(class(out),"blankcorrected")
    return(out)
}
#' @rdname blankcorr
#' @export
blankcorr.PHdata <- function(x,blanklabel=NULL,prefix='',...){
    out <- x
    for (mass in out$masses){
        out$signals[[mass]] <-
            blankcorr.timeresolved(out$signals[[mass]],blanklabel,prefix)
    }
    return(out)
}
#' @rdname blankcorr
#' @export
blankcorr.WiscAr <- function(x,blanklabel=NULL,prefix='',...){
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@10@"]]));##:ess-bp-end:##
    iblanks <- array(grep(blanklabel,x$labels))
}
